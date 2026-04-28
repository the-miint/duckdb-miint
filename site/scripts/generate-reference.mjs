#!/usr/bin/env node
// Read site/build/functions.json and emit one MD page per function name into
// site/src/content/docs/reference/<category>/<name>.md, plus category index
// pages that group related functions in a single table.
//
// POC scope: only the slice we are documenting in C++ for Daniel's review.
// Functions outside this slice are skipped silently — they will be added once
// the convention is approved and the rest of the codebase is migrated.

import { mkdirSync, writeFileSync, readFileSync, rmSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const IN = resolve(__dirname, '..', 'build', 'functions.json');
const OUT_ROOT = resolve(__dirname, '..', 'src', 'content', 'docs', 'reference');

// Functions whose C++ registrations have been converted to use the
// RegisterDocumented* helper. Anything not in this set is silently skipped
// during generation.
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

const CATEGORY_DIR = {
  scalar: 'scalar-functions',
  aggregate: 'aggregate-functions',
  table: 'table-functions',
  table_macro: 'table-macros',
  macro: 'macros',
};

const data = JSON.parse(readFileSync(IN, 'utf8'));

// Wipe and recreate so removed/renamed functions don't linger.
for (const sub of Object.values(CATEGORY_DIR)) {
  const dir = resolve(OUT_ROOT, sub);
  rmSync(dir, { recursive: true, force: true });
}

let written = 0;
let skipped = 0;
let missingMeta = 0;
const errors = [];

// First pass: collect categories so we can render an index per category and
// link from each page.
const byCategory = new Map(); // category -> Array<fn>
const pocFunctions = data.functions.filter((fn) => POC_SLICE.has(fn.name));
for (const fn of pocFunctions) {
  const cats = fn.variants[0]?.categories ?? [];
  for (const c of cats) {
    if (!byCategory.has(c)) byCategory.set(c, []);
    byCategory.get(c).push(fn);
  }
}

// Second pass: per-function pages. Functions with a category live under
// <type-dir>/<category>/<name>.md so the category index page (next pass) is
// their true URL parent. Functions with no category fall back to <type-dir>/.
for (const fn of pocFunctions) {
  const subdir = CATEGORY_DIR[fn.type];
  if (!subdir) {
    errors.push(`unsupported function_type '${fn.type}' for '${fn.name}'`);
    continue;
  }
  const primary = fn.variants.find((v) => v.description) ?? fn.variants[0];
  if (!primary.description) missingMeta++;

  const cats = primary.categories ?? [];
  const category = cats[0] ?? null; // first category wins as the canonical home
  const outDir = category
    ? resolve(OUT_ROOT, subdir, category)
    : resolve(OUT_ROOT, subdir);
  mkdirSync(outDir, { recursive: true });
  writeFileSync(resolve(outDir, `${fn.name}.md`), renderPage(fn, primary, category));
  written++;
}

// Third pass: category index pages. Each index becomes the parent URL of the
// functions tagged with that category: <type-dir>/<category>/index.md. When a
// category spans multiple function types we'd need to revisit this, but for
// now every category is single-type.
let indexes = 0;
for (const [category, fns] of byCategory.entries()) {
  const types = new Set(fns.map((f) => f.type));
  if (types.size > 1) {
    errors.push(`category '${category}' spans multiple function types: ${[...types].join(', ')}`);
    continue;
  }
  const subdir = CATEGORY_DIR[fns[0].type];
  const outDir = resolve(OUT_ROOT, subdir, category);
  mkdirSync(outDir, { recursive: true });
  writeFileSync(resolve(outDir, 'index.md'), renderCategoryIndex(category, fns));
  indexes++;
}

if (errors.length) {
  for (const e of errors) console.error('generate:', e);
  process.exit(1);
}

// Emit a sidebar manifest the Astro config consumes. Manual groups beat
// autogenerate here because autogenerate produces a duplicate entry for any
// folder that also contains an index.md (folder name + index title both show).
const TYPE_LABEL = {
  scalar: 'Scalar functions',
  aggregate: 'Aggregate functions',
  table: 'Table functions',
  table_macro: 'Table macros',
  macro: 'Macros',
};
const sidebar = [];
const UNCATEGORIZED = '_uncategorized';
const fnsByTypeAndCategory = new Map(); // type -> category -> [fn]
for (const fn of pocFunctions) {
  const primary = fn.variants.find((v) => v.description) ?? fn.variants[0];
  // Functions with no category still belong in the sidebar — bucket them
  // under a synthetic UNCATEGORIZED group so released-mode introspection
  // (where the registered functions don't yet expose categories) still
  // produces a navigable sidebar instead of an empty Reference section.
  const cat = (primary.categories ?? [])[0] ?? UNCATEGORIZED;
  if (!fnsByTypeAndCategory.has(fn.type)) fnsByTypeAndCategory.set(fn.type, new Map());
  const m = fnsByTypeAndCategory.get(fn.type);
  if (!m.has(cat)) m.set(cat, []);
  m.get(cat).push(fn);
}
// Stable sidebar order: table functions, then scalar, then aggregates, etc.
// Matches what users tend to want (table-shape APIs are the "front door").
const TYPE_ORDER = ['table', 'scalar', 'aggregate', 'table_macro', 'macro'];
const orderedTypes = TYPE_ORDER.filter((t) => fnsByTypeAndCategory.has(t));
for (const type of orderedTypes) {
  const cats = fnsByTypeAndCategory.get(type);
  const subdir = CATEGORY_DIR[type];
  const groups = [];
  for (const [cat, fns] of cats.entries()) {
    const isUncategorized = cat === UNCATEGORIZED;
    // Per-function page paths must match what the page-write pass produced:
    // categorized -> reference/<subdir>/<cat>/<name>; uncategorized -> reference/<subdir>/<name>
    const fnItems = fns
      .slice()
      .sort((a, b) => a.name.localeCompare(b.name))
      .map((f) => ({
        label: f.name,
        slug: isUncategorized
          ? `reference/${subdir}/${f.name}`
          : `reference/${subdir}/${cat}/${f.name}`,
      }));
    if (isUncategorized) {
      // Uncategorized functions go directly under the type-level group with no
      // sub-collapse and no Overview link (no category index page exists).
      groups.push({
        label: 'All',
        collapsed: false,
        items: fnItems,
      });
    } else {
      // Starlight ManualSidebarGroupSchema is a strictObject — `link` on a group
      // is rejected. To make the category index reachable from the sidebar we
      // include an "Overview" item that points at it.
      groups.push({
        label: humanizeCategory(cat),
        collapsed: true,
        items: [{ label: 'Overview', slug: `reference/${subdir}/${cat}` }, ...fnItems],
      });
    }
  }
  // Sort but keep "All" (uncategorized) at the bottom.
  groups.sort((a, b) => {
    if (a.label === 'All') return 1;
    if (b.label === 'All') return -1;
    return a.label.localeCompare(b.label);
  });
  sidebar.push({ label: TYPE_LABEL[type] ?? type, items: groups });
}
writeFileSync(resolve(__dirname, '..', 'build', 'sidebar.json'), JSON.stringify(sidebar, null, 2));

console.error(`generate: wrote ${written} pages + ${indexes} category indexes; skipped ${skipped} (outside POC slice)`);
if (missingMeta) {
  console.error(`generate: WARNING ${missingMeta} POC functions still lack description in C++`);
}

function renderPage(fn, primary, category) {
  const summary = primary.description
    ? primary.description.trim().split(/\n\s*\n/)[0].replace(/\s+/g, ' ').slice(0, 160)
    : `Reference for ${fn.name}.`;

  const header = [
    '---',
    `title: ${fn.name}`,
    `description: ${escapeYaml(summary)}`,
    '---',
    '',
    `<!-- Auto-generated from duckdb_functions(). Edit the C++ registration to change this page. -->`,
    '',
  ].join('\n');

  // Alias points at canonical (sibling page in the same category directory).
  const aliasLine = fn.alias_of
    ? `:::note[Alias]\n\`${fn.name}\` is an alias of [\`${fn.alias_of}\`](../${fn.alias_of}/). They behave identically; this page is included so the docs match what \`SHOW FUNCTIONS\` reports.\n:::\n\n`
    : '';

  const intro = primary.description
    ? `${primary.description.trim()}\n\n`
    : `_No description set yet. Add one to the C++ registration._\n\n`;

  const signatures = fn.variants.map((v, i) => {
    const params = v.parameters.map((p, j) => {
      const t = v.parameter_types[j] ?? 'ANY';
      const name = /^col\d+$/.test(p) ? `_${p}_` : p;
      return `${name} ${t}`;
    }).join(', ');
    const ret = v.return_type ?? 'TABLE';
    const heading = fn.variants.length > 1 ? `### Signature ${i + 1}\n\n` : '';
    return `${heading}\`\`\`text\n${fn.name}(${params}) → ${ret}\n\`\`\`\n`;
  }).join('\n');

  // Examples: collect from all variants, dedup by trimmed text.
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

  // Back-link to the parent category index. Function pages live at
  // <type>/<category>/<name>/, so the index is one segment up.
  const seeAlso = category
    ? `\n## See also\n\n- [${humanizeCategory(category)} (overview)](../)\n`
    : '';

  return header + aliasLine + intro + '## Signature\n\n' + signatures + '\n' + examples + seeAlso;
}

function renderCategoryIndex(category, fns) {
  const title = humanizeCategory(category);

  // Group by canonical (alias_of: null) and attach aliases.
  const canonical = fns.filter((f) => !f.alias_of);
  const aliasesOf = new Map(); // canonical name -> [alias names]
  for (const f of fns) {
    if (f.alias_of) {
      if (!aliasesOf.has(f.alias_of)) aliasesOf.set(f.alias_of, []);
      aliasesOf.get(f.alias_of).push(f.name);
    }
  }

  const rows = canonical
    .sort((a, b) => a.name.localeCompare(b.name))
    .map((f) => {
      const desc = f.variants[0]?.description?.trim().split(/\n\s*\n/)[0].replace(/\s+/g, ' ') ?? '';
      const summary = desc.length > 120 ? desc.slice(0, 117) + '…' : desc;
      const aliases = (aliasesOf.get(f.name) ?? []).map((a) => `\`${a}\``).join(', ') || '—';
      // Index is at /<type>/<category>/, function pages are children at
      // /<type>/<category>/<name>/, so the link is `./<name>/`.
      return `| [\`${f.name}\`](./${f.name}/) | ${summary} | ${aliases} |`;
    });

  const header = [
    '---',
    `title: ${title}`,
    `description: All ${fns.length} ${title.toLowerCase()} registered by miint.`,
    '---',
    '',
    `<!-- Auto-generated index for category '${category}'. -->`,
    '',
    `Functions tagged \`${category}\`. Click any name for full reference.`,
    '',
    '| Function | Description | Aliases |',
    '|---|---|---|',
    ...rows,
    '',
  ];
  return header.join('\n');
}

function humanizeCategory(cat) {
  // sam-flags -> SAM flags; alignment-quality -> Alignment quality
  return cat.split(/[-_]/).map((w, i) => {
    if (i === 0) return w.charAt(0).toUpperCase() + w.slice(1);
    if (/^[A-Z]+$/.test(w.toUpperCase()) && w.length <= 4) return w.toUpperCase();
    return w;
  }).join(' ').replace(/\bSam\b/g, 'SAM');
}

function escapeYaml(s) {
  if (/[:#\n"']/.test(s)) {
    return `"${s.replace(/\\/g, '\\\\').replace(/"/g, '\\"')}"`;
  }
  return s;
}
