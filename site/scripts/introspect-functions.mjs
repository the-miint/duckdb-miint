#!/usr/bin/env node
// Introspect the miint extension by diffing duckdb_functions() before and
// after LOAD. Writes site/build/functions.json.
//
// Inputs (env vars override defaults):
//   MIINT_EXT  path to a built miint.duckdb_extension matching DUCKDB_BIN's version
//   DUCKDB_BIN stock duckdb client of the matching version (NOT the repo build,
//              which auto-loads miint and so cannot snapshot a clean baseline)
//
// Output:
//   site/build/functions.json
//
// Why a diff and not WHERE internal=false? The miint extension registers
// functions with internal=true (same as core extensions), so a filter on
// internal would return zero miint functions.

import { spawnSync } from 'node:child_process';
import { mkdirSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = resolve(__dirname, '..', '..');

const DUCKDB_BIN = process.env.DUCKDB_BIN
  ?? '/home/lpatel/software/miniconda3/envs/miint/bin/duckdb';
const MIINT_EXT = process.env.MIINT_EXT
  ?? resolve(REPO_ROOT, 'build/release/extension/miint/miint.duckdb_extension');
const OUT = resolve(__dirname, '..', 'build', 'functions.json');

function runDuckDB(sql) {
  const res = spawnSync(DUCKDB_BIN, ['-unsigned', '-json', '-c', sql], {
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

function snapshot(loadStmt) {
  const select = `SELECT ${FIELDS.join(', ')} FROM duckdb_functions();`;
  return runDuckDB(loadStmt ? `${loadStmt}; ${select}` : select);
}

console.error(`introspect: baseline (no extension) via ${DUCKDB_BIN}`);
const before = snapshot('');
const beforeKey = new Set(
  before.map((r) => `${r.function_name} ${r.function_type} ${(r.parameter_types ?? []).join(',')}`)
);

console.error(`introspect: after LOAD '${MIINT_EXT}'`);
const after = snapshot(`LOAD '${MIINT_EXT}'`);

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
  miint_extension: MIINT_EXT,
  summary,
  functions,
}, null, 2));

console.error(`introspect: wrote ${OUT}`);
console.error(`introspect: ${functions.length} unique function names; by type:`, summary);
