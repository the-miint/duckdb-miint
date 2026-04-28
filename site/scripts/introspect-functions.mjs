#!/usr/bin/env node
// Introspect the miint extension by diffing duckdb_functions() before and
// after the extension is loaded. Writes site/build/functions.json.
//
// Default behaviour:
//   - If a local extension build exists at
//     <repo>/build/release/extension/miint/miint.duckdb_extension, introspect
//     that. (Most likely what a developer working in the repo wants.)
//   - Otherwise install the released extension from the public miint
//     repository at https://ftp.microbio.me/pub/miint. (CI / clean checkouts
//     without a build artifact.)
//
// Inputs (env vars override defaults):
//   DUCKDB_BIN   stock duckdb client (NOT the repo's build/release/duckdb,
//                which auto-loads miint and so cannot snapshot a clean
//                baseline). Default: 'duckdb' (resolved via PATH).
//   MIINT_REPO   Custom miint extension repository to install from.
//                Default: 'https://ftp.microbio.me/pub/miint'.
//   MIINT_EXT    Explicit path to a .duckdb_extension to introspect.
//                Overrides the local-build auto-detection.
//   MIINT_FORCE_RELEASED  When set to "1" or "true", skip the local-build
//                auto-detection and always introspect the released artifact.
//                Useful for testing what the docs will look like in CI.
//
// Output:
//   site/build/functions.json
//
// Why a diff and not WHERE internal=false? The miint extension registers
// functions with internal=true (same as core extensions), so a filter on
// internal would return zero miint functions. The diff isolates exactly the
// catalog entries miint adds.

import { spawnSync } from 'node:child_process';
import { existsSync, mkdirSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = resolve(__dirname, '..', '..');
const DEFAULT_LOCAL_EXT = resolve(REPO_ROOT, 'build/release/extension/miint/miint.duckdb_extension');

const DUCKDB_BIN = process.env.DUCKDB_BIN ?? 'duckdb';
const MIINT_REPO = process.env.MIINT_REPO ?? 'https://ftp.microbio.me/pub/miint';
const FORCE_RELEASED = ['1', 'true', 'yes'].includes((process.env.MIINT_FORCE_RELEASED ?? '').toLowerCase());
// Resolution order:
//   1. MIINT_EXT explicit path (highest priority)
//   2. MIINT_FORCE_RELEASED=1 -> released repo
//   3. local build at REPO_ROOT/build/... if it exists -> local
//   4. released repo (fallback)
const MIINT_EXT = process.env.MIINT_EXT
  ?? (!FORCE_RELEASED && existsSync(DEFAULT_LOCAL_EXT) ? DEFAULT_LOCAL_EXT : null);
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

// Fields generate-reference.mjs actually consumes. Don't add fields here
// without a consumer — they bloat functions.json for no benefit.
const FIELDS = [
  'function_name', 'function_type', 'alias_of', 'description',
  'parameters', 'parameter_types', 'return_type',
  'examples', 'categories',
];

if (useLocalBuild) {
  const tag = process.env.MIINT_EXT ? 'local, MIINT_EXT explicit' : 'local, auto-detected';
  console.error(`introspect: ${DUCKDB_BIN} + LOAD '${MIINT_EXT}' (${tag})`);
} else {
  const tag = FORCE_RELEASED ? 'released, MIINT_FORCE_RELEASED' : 'released, no local build found';
  console.error(`introspect: ${DUCKDB_BIN} + INSTALL miint FROM '${MIINT_REPO}' (${tag})`);
}

// Single-session diff via a temp table: snapshot the baseline catalog,
// load miint, then return only the rows that weren't in the baseline.
// One subprocess instead of two — saves a duckdb cold start (~150-400ms)
// and the JSON serialize/parse round-trip for the unused baseline.
const newFns = runDuckDB(`
CREATE TEMP TABLE _miint_doc_baseline AS
  SELECT function_name, function_type, parameter_types FROM duckdb_functions();
${loadStmt};
SELECT ${FIELDS.join(', ')}
FROM duckdb_functions() f
WHERE NOT EXISTS (
  SELECT 1 FROM _miint_doc_baseline b
  WHERE b.function_name = f.function_name
    AND b.function_type = f.function_type
    AND b.parameter_types IS NOT DISTINCT FROM f.parameter_types
);
`);

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
