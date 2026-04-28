#!/usr/bin/env node
// Read site/build/functions.json and emit one MD page per function name into
// site/src/content/docs/reference/<type>/<category>/<name>.md, plus category
// index pages and a per-type "module landing" page.
//
// POC scope: only the slice we are documenting in C++ for Daniel's review.
// Functions outside this slice are skipped silently — they will be added once
// the convention is approved and the rest of the codebase is migrated.

import { existsSync, mkdirSync, writeFileSync, readFileSync, rmSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import yaml from 'js-yaml';

const __dirname = dirname(fileURLToPath(import.meta.url));
const IN = resolve(__dirname, '..', 'build', 'functions.json');
const OUT_ROOT = resolve(__dirname, '..', 'src', 'content', 'docs', 'reference');
const OVERVIEW_CONTENT_ROOT = resolve(__dirname, '..', 'src', 'overview-content');

// Functions whose C++ registrations have been converted to the
// RegisterDocumented* helper. Anything outside this set is silently skipped.
const POC_SLICE = new Set([
  // alignment_flag_functions.cpp
  'alignment_is_paired', 'is_paired',
  'alignment_is_proper_pair', 'is_proper_pair',
  'alignment_is_unmapped', 'is_unmapped',
  'alignment_is_mate_unmapped', 'is_munmap',
  'alignment_is_reverse', 'is_reverse',
  'alignment_is_mate_reverse', 'is_mreverse',
  'alignment_is_read1', 'is_read1',
  'alignment_is_read2', 'is_read2',
  'alignment_is_secondary', 'is_secondary',
  'alignment_is_primary',
  'alignment_is_qc_failed', 'is_qcfail',
  'alignment_is_duplicate', 'is_dup',
  'alignment_is_supplementary', 'is_supplementary',
  // alignment_functions.cpp
  'alignment_seq_identity',
  'alignment_query_length',
  'alignment_query_coverage',
  // read_alignments.cpp
  'read_alignments', 'read_sam',
]);

// Per-function-type metadata. Keys must match `function_type` values from
// duckdb_functions(). Order in TYPE_ORDER drives sidebar order.
const TYPES = {
  table:       { dir: 'table-functions',     label: 'Table functions'     },
  scalar:      { dir: 'scalar-functions',    label: 'Scalar functions'    },
  aggregate:   { dir: 'aggregate-functions', label: 'Aggregate functions' },
  table_macro: { dir: 'table-macros',        label: 'Table macros'        },
  macro:       { dir: 'macros',              label: 'Macros'              },
};
const TYPE_ORDER = ['table', 'scalar', 'aggregate', 'table_macro', 'macro'];
// Acronyms preserved as uppercase by humanizeCategory.
const ACRONYMS = new Set(['sam', 'bam', 'io', 'qc', 'cigar', 'md', 'nm', 'rna', 'dna']);
const UNCATEGORIZED = '_uncategorized';

// Sanity: every type in TYPE_ORDER has a TYPES entry.
for (const t of TYPE_ORDER) {
  if (!TYPES[t]) throw new Error(`generate: TYPE_ORDER references unknown type '${t}'`);
}

// Load a hand-authored markdown fragment that prepends or appends to an
// auto-generated overview page. Returns '' if no fragment file exists.
function loadOverviewFragment(kind, key, slot) {
  const p = resolve(OVERVIEW_CONTENT_ROOT, kind, `${key}.${slot}.md`);
  if (!existsSync(p)) return '';
  return readFileSync(p, 'utf8').trim() + '\n\n';
}

const data = JSON.parse(readFileSync(IN, 'utf8'));

// Wipe and recreate so removed/renamed functions don't linger.
for (const t of Object.values(TYPES)) {
  rmSync(resolve(OUT_ROOT, t.dir), { recursive: true, force: true });
}

const errors = [];

// Filter to POC slice and enrich each function with its picked variant +
// home category (first wins). Single source of truth for downstream passes.
const pocFunctions = data.functions
  .filter((fn) => POC_SLICE.has(fn.name))
  .map((fn) => {
    const primary = fn.variants.find((v) => v.description) ?? fn.variants[0];
    const category = (primary.categories ?? [])[0] ?? UNCATEGORIZED;
    return { ...fn, primary, category };
  });

// Bucket by type → category. Drives both the per-function page layout and
// the sidebar / index pages below — one map, no derived twins.
const fnsByTypeAndCategory = new Map();
for (const fn of pocFunctions) {
  if (!TYPES[fn.type]) {
    errors.push(`unsupported function_type '${fn.type}' for '${fn.name}'`);
    continue;
  }
  if (!fnsByTypeAndCategory.has(fn.type)) fnsByTypeAndCategory.set(fn.type, new Map());
  const m = fnsByTypeAndCategory.get(fn.type);
  if (!m.has(fn.category)) m.set(fn.category, []);
  m.get(fn.category).push(fn);
}

// Per-function pages. Categorized fns live at <type-dir>/<category>/<name>/;
// uncategorized fall back to <type-dir>/<name>/.
let written = 0;
let missingMeta = 0;
for (const fn of pocFunctions) {
  if (!TYPES[fn.type]) continue;
  if (!fn.primary.description) missingMeta++;
  const subdir = TYPES[fn.type].dir;
  const isUncategorized = fn.category === UNCATEGORIZED;
  const outDir = isUncategorized ? resolve(OUT_ROOT, subdir) : resolve(OUT_ROOT, subdir, fn.category);
  mkdirSync(outDir, { recursive: true });
  writeFileSync(resolve(outDir, `${fn.name}.md`),
                renderPage(fn, fn.primary, isUncategorized ? null : fn.category));
  written++;
}

// Per-category index pages. Skipped for the synthetic UNCATEGORIZED bucket
// (the type-level overview already lists those flat functions).
let categoryIndexes = 0;
for (const [type, cats] of fnsByTypeAndCategory) {
  for (const [category, fns] of cats) {
    if (category === UNCATEGORIZED) continue;
    const outDir = resolve(OUT_ROOT, TYPES[type].dir, category);
    mkdirSync(outDir, { recursive: true });
    writeFileSync(resolve(outDir, 'index.md'), renderCategoryIndex(category, fns));
    categoryIndexes++;
  }
}

// Per-type "module landing" pages. scikit-bio analog:
// https://scikit.bio/docs/latest/alignment.html
let typeIndexes = 0;
for (const [type, cats] of fnsByTypeAndCategory) {
  const outDir = resolve(OUT_ROOT, TYPES[type].dir);
  mkdirSync(outDir, { recursive: true });
  writeFileSync(resolve(outDir, 'index.md'), renderTypeIndex(type, cats));
  typeIndexes++;
}

// Sidebar manifest. Each type group always opens with "Overview" linking to
// the type-level index, then per-category subgroups (each with their own
// Overview), then any uncategorized functions inline. No "All" wrapper.
const sidebar = [];
for (const type of TYPE_ORDER) {
  if (!fnsByTypeAndCategory.has(type)) continue;
  const cats = fnsByTypeAndCategory.get(type);
  const subdir = TYPES[type].dir;
  const items = [{ label: 'Overview', slug: `reference/${subdir}` }];

  const realCats = [...cats.keys()].filter((c) => c !== UNCATEGORIZED).sort();
  for (const cat of realCats) {
    items.push({
      label: humanizeCategory(cat),
      collapsed: true,
      items: [
        { label: 'Overview', slug: `reference/${subdir}/${cat}` },
        ...sortedFnSlugs(cats.get(cat), `reference/${subdir}/${cat}`),
      ],
    });
  }
  if (cats.has(UNCATEGORIZED)) {
    items.push(...sortedFnSlugs(cats.get(UNCATEGORIZED), `reference/${subdir}`));
  }
  sidebar.push({ label: TYPES[type].label, items });
}
writeFileSync(resolve(__dirname, '..', 'build', 'sidebar.json'), JSON.stringify(sidebar, null, 2));

if (errors.length) {
  for (const e of errors) console.error('generate:', e);
  process.exit(1);
}

console.error(`generate: wrote ${written} pages + ${categoryIndexes} category indexes + ${typeIndexes} type indexes`);
if (missingMeta) {
  console.error(`generate: WARNING ${missingMeta} POC functions still lack description in C++`);
}

// ---------- helpers ----------

function sortedFnSlugs(fns, prefix) {
  return fns.slice()
    .sort((a, b) => a.name.localeCompare(b.name))
    .map((f) => ({ label: f.name, slug: `${prefix}/${f.name}` }));
}

// canonical-name -> [alias names]
function buildAliasesOf(fns) {
  const map = new Map();
  for (const f of fns) {
    if (!f.alias_of) continue;
    if (!map.has(f.alias_of)) map.set(f.alias_of, []);
    map.get(f.alias_of).push(f.name);
  }
  return map;
}

function shortDesc(fn, maxLen = 120) {
  const desc = fn.primary?.description?.trim().split(/\n\s*\n/)[0].replace(/\s+/g, ' ') ?? '';
  return desc.length > maxLen ? desc.slice(0, maxLen - 1) + '…' : desc;
}

function renderFunctionRows(fns, hrefPrefix = './', maxLen = 120) {
  const aliasesOf = buildAliasesOf(fns);
  const canonical = fns.filter((f) => !f.alias_of).sort((a, b) => a.name.localeCompare(b.name));
  return canonical.map((f) => {
    const aliases = (aliasesOf.get(f.name) ?? []).map((a) => `\`${a}\``).join(', ') || '—';
    return `| [\`${f.name}\`](${hrefPrefix}${f.name}/) | ${shortDesc(f, maxLen)} | ${aliases} |`;
  });
}

// YAML frontmatter — js-yaml handles all the escape edge cases (leading -,
// quotes inside, multi-line, ambiguous yes/no/null/numbers, etc.) that a
// hand-rolled escape would miss.
function frontmatter(obj) {
  return `---\n${yaml.dump(obj, { lineWidth: 0 })}---\n\n`;
}

function renderPage(fn, primary, category) {
  const summaryLine = (primary.description?.trim().split(/\n\s*\n/)[0].replace(/\s+/g, ' ') ?? '').slice(0, 160)
    || `Reference for ${fn.name}.`;

  const header = frontmatter({ title: fn.name, description: summaryLine })
    + `<!-- Auto-generated from duckdb_functions(). Edit the C++ registration to change this page. -->\n\n`;

  const aliasLine = fn.alias_of
    ? `:::note[Alias]\n\`${fn.name}\` is an alias of [\`${fn.alias_of}\`](../${fn.alias_of}/). They behave identically; this page is included so the docs match what \`SHOW FUNCTIONS\` reports.\n:::\n\n`
    : '';

  const body = primary.description
    ? `${primary.description.trim()}\n\n`
    : `_No description set yet. Add one to the C++ registration._\n\n`;

  const signatures = fn.variants.map((v, i) => {
    const params = v.parameters.map((p, j) => {
      const t = v.parameter_types[j] ?? 'ANY';
      // DuckDB returns col0/col1/... when parameter_names weren't set in
      // FunctionDescription — flag visually so the gap is obvious.
      const name = /^col\d+$/.test(p) ? `_${p}_` : p;
      return `${name} ${t}`;
    }).join(', ');
    const ret = v.return_type ?? 'TABLE';
    const heading = fn.variants.length > 1 ? `### Signature ${i + 1}\n\n` : '';
    return `${heading}\`\`\`text\n${fn.name}(${params}) → ${ret}\n\`\`\`\n`;
  }).join('\n');

  const seen = new Set();
  const allExamples = [];
  for (const v of fn.variants) {
    for (const e of (v.examples ?? [])) {
      const k = e.trim();
      if (!seen.has(k)) { seen.add(k); allExamples.push(k); }
    }
  }
  const examples = allExamples.length
    ? `## Examples\n\n${allExamples.map((e) => '```sql\n' + e + '\n```').join('\n\n')}\n`
    : '';

  // Function pages live at <type>/<category>/<name>/; the index is `../`.
  const seeAlso = category
    ? `\n## See also\n\n- [${humanizeCategory(category)} (overview)](../)\n`
    : '';

  return header + aliasLine + body + '## Signature\n\n' + signatures + '\n' + examples + seeAlso;
}

function renderCategoryIndex(category, fns) {
  const title = humanizeCategory(category);
  const intro = loadOverviewFragment('categories', category, 'intro');
  const notes = loadOverviewFragment('categories', category, 'notes');

  const lines = [
    frontmatter({ title, description: `All ${fns.length} ${title.toLowerCase()} registered by miint.` }).trimEnd(),
    '',
    `<!-- Auto-generated index for category '${category}'. Hand-authored prose, if any, lives in site/src/overview-content/categories/${category}.{intro,notes}.md -->`,
    '',
  ];
  if (intro) {
    lines.push(intro);
  } else {
    lines.push(`Functions tagged \`${category}\`. Click any name for full reference.`, '');
  }
  lines.push('| Function | Description | Aliases |', '|---|---|---|', ...renderFunctionRows(fns), '');
  if (notes) lines.push(notes);
  return lines.join('\n');
}

// Type-level overview — scikit-bio "module landing". Lists categories with
// summaries (each with its own overview link), and any uncategorized
// functions in a flat table.
function renderTypeIndex(type, catsMap) {
  const { label: title } = TYPES[type];
  const realCats = [...catsMap.keys()].filter((c) => c !== UNCATEGORIZED).sort();
  const flatFns = catsMap.get(UNCATEGORIZED) ?? [];
  const totalFns = [...catsMap.values()].reduce(
    (acc, fns) => acc + fns.filter((f) => !f.alias_of).length, 0);

  const intro = loadOverviewFragment('types', type, 'intro');
  const notes = loadOverviewFragment('types', type, 'notes');

  const lines = [
    frontmatter({ title, description: `All ${totalFns} ${title.toLowerCase()} registered by miint (POC slice).` }).trimEnd(),
    '',
    `<!-- Auto-generated index for function type '${type}'. Hand-authored prose, if any, lives in site/src/overview-content/types/${type}.{intro,notes}.md -->`,
    '',
  ];

  if (intro) lines.push(intro);

  if (realCats.length) {
    lines.push(
      `Functions are grouped by category. Each category has its own overview page; click any function name for full reference.`,
      '',
      '| Category | Description | Functions |',
      '|---|---|---|',
    );
    for (const cat of realCats) {
      const fns = catsMap.get(cat);
      const canonicalCount = fns.filter((f) => !f.alias_of).length;
      const aliasCount = fns.length - canonicalCount;
      const aliasNote = aliasCount ? ` (+${aliasCount} alias${aliasCount === 1 ? '' : 'es'})` : '';
      const sample = fns[0]?.primary?.description?.trim().split(/\n\s*\n/)[0].replace(/\s+/g, ' ') ?? '';
      const summary = sample.length > 110 ? sample.slice(0, 109) + '…' : sample;
      lines.push(`| [${humanizeCategory(cat)}](./${cat}/) | ${summary} | ${canonicalCount}${aliasNote} |`);
    }
    lines.push('');
  }

  if (flatFns.length) {
    if (realCats.length) lines.push('## Other functions', '');
    lines.push('| Function | Description | Aliases |', '|---|---|---|', ...renderFunctionRows(flatFns), '');
  }

  if (!realCats.length && !flatFns.length) {
    lines.push(`_No ${title.toLowerCase()} are documented yet._`);
  }

  if (notes) lines.push(notes);
  return lines.join('\n');
}

// sam-flags -> "SAM flags"; alignment-quality -> "Alignment quality".
// Words in ACRONYMS are uppercased; the first word otherwise gets title case.
function humanizeCategory(cat) {
  return cat.split(/[-_]+/).map((w, i) => {
    if (ACRONYMS.has(w.toLowerCase())) return w.toUpperCase();
    return i === 0 ? w.charAt(0).toUpperCase() + w.slice(1) : w;
  }).join(' ');
}
