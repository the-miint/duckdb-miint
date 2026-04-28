#!/usr/bin/env node
// Introspect the miint extension by diffing duckdb_functions() before and
// after the extension is loaded. Writes site/build/functions.json.
//
// By default the docs build introspects the *released* extension from the
// public miint repository at https://ftp.microbio.me/pub/miint, so the docs
// reflect what users actually install. To introspect a local development
// build instead, set MIINT_EXT to the path of a built .duckdb_extension.
//
// Inputs (env vars override defaults):
//   DUCKDB_BIN   stock duckdb client (NOT the repo's build/release/duckdb,
//                which auto-loads miint and so cannot snapshot a clean
//                baseline). Default: 'duckdb' (resolved via PATH).
//   MIINT_REPO   Custom miint extension repository to install from.
//                Default: 'https://ftp.microbio.me/pub/miint'.
//   MIINT_EXT    Path to a local .duckdb_extension to introspect instead of
//                the released artifact. When set, MIINT_REPO is ignored and
//                duckdb is invoked with -unsigned. Useful for previewing docs
//                while iterating on a registration locally.
//
// Output:
//   site/build/functions.json
//
// Why a diff and not WHERE internal=false? The miint extension registers
// functions with internal=true (same as core extensions), so a filter on
// internal would return zero miint functions. The diff isolates exactly the
// catalog entries miint adds.

import { spawnSync } from 'node:child_process';
import { mkdirSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = dirname(fileURLToPath(import.meta.url));

const DUCKDB_BIN = process.env.DUCKDB_BIN ?? 'duckdb';
const MIINT_REPO = process.env.MIINT_REPO ?? 'https://ftp.microbio.me/pub/miint';
const MIINT_EXT = process.env.MIINT_EXT ?? null; // null => use MIINT_REPO
const OUT = resolve(__dirname, '..', 'build', 'functions.json');

const useLocalBuild = MIINT_EXT !== null;
// FORCE INSTALL: silently overrides any existing community-installed miint
// pointing at a different repository, so the same tooling works whether the
// developer's local DuckDB has miint pre-installed from community or not.
// httpfs is loaded first because the released repo is HTTPS-served and
// httpfs does not autoload reliably in every DuckDB packaging.
const loadStmt = useLocalBuild
  ? `LOAD '${MIINT_EXT}'`
  : `INSTALL httpfs; LOAD httpfs; FORCE INSTALL miint FROM '${MIINT_REPO}'; LOAD miint`;

function runDuckDB(sql) {
  // -unsigned is needed for both modes: the local-build path is unsigned,
  // and the released artifact at MIINT_REPO is not signed by a key the
  // stock duckdb client recognizes for that URL. (Community-extensions
  // installs ARE signed; the dedicated miint repo is not yet.)
  const args = ['-unsigned', '-json', '-c', sql];
  const res = spawnSync(DUCKDB_BIN, args, {
    encoding: 'utf8',
    maxBuffer: 64 * 1024 * 1024,
  });
  if (res.status !== 0) {
    process.stderr.write(`duckdb exited ${res.status}\nSTDERR:\n${res.stderr}\n`);
    process.exit(1);
  }
  const out = res.stdout.trim();
  if (!out) return [];
  return JSON.parse(out);
}

const FIELDS = [
  'function_name', 'function_type', 'alias_of', 'description',
  'parameters', 'parameter_types', 'return_type', 'varargs',
  'has_side_effects', 'examples', 'categories',
];

function snapshot(load) {
  const select = `SELECT ${FIELDS.join(', ')} FROM duckdb_functions();`;
  return runDuckDB(load ? `${load}; ${select}` : select);
}

console.error(`introspect: baseline (no extension) via ${DUCKDB_BIN}`);
const before = snapshot('');
const beforeKey = new Set(
  before.map((r) => `${r.function_name} ${r.function_type} ${(r.parameter_types ?? []).join(',')}`)
);

if (useLocalBuild) {
  console.error(`introspect: after LOAD '${MIINT_EXT}' (local-build override)`);
} else {
  console.error(`introspect: after INSTALL miint FROM '${MIINT_REPO}' (released)`);
}
const after = snapshot(loadStmt);

const newFns = after.filter(
  (r) => !beforeKey.has(`${r.function_name} ${r.function_type} ${(r.parameter_types ?? []).join(',')}`)
);

// Group overloads under a single name so the generator emits one page per name.
const grouped = new Map();
for (const fn of newFns) {
  const key = `${fn.function_type}:${fn.function_name}`;
  if (!grouped.has(key)) {
    grouped.set(key, {
      name: fn.function_name,
      type: fn.function_type,
      alias_of: fn.alias_of || null,
      variants: [],
    });
  }
  grouped.get(key).variants.push({
    parameters: fn.parameters ?? [],
    parameter_types: fn.parameter_types ?? [],
    return_type: fn.return_type ?? null,
    description: fn.description ?? null,
    examples: fn.examples ?? [],
    categories: fn.categories ?? [],
    varargs: fn.varargs ?? null,
    has_side_effects: fn.has_side_effects ?? false,
  });
}

// Stable sort: alias_of-having entries last, then by name.
const functions = [...grouped.values()].sort((a, b) => {
  const ax = a.alias_of ? 1 : 0;
  const bx = b.alias_of ? 1 : 0;
  if (ax !== bx) return ax - bx;
  return a.name.localeCompare(b.name);
});

const summary = functions.reduce((acc, f) => {
  acc[f.type] = (acc[f.type] ?? 0) + 1;
  return acc;
}, {});

mkdirSync(dirname(OUT), { recursive: true });
writeFileSync(OUT, JSON.stringify({
  generated_at: new Date().toISOString(),
  duckdb_bin: DUCKDB_BIN,
  miint_source: useLocalBuild ? { kind: 'local', path: MIINT_EXT } : { kind: 'released', repo: MIINT_REPO },
  summary,
  functions,
}, null, 2));

console.error(`introspect: wrote ${OUT}`);
console.error(`introspect: ${functions.length} unique function names; by type:`, summary);
