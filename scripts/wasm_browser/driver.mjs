// Run steps.json against the miint extension and compare each step to its
// expected value.
//
//   node driver.mjs --wasm <miint.duckdb_extension.wasm> --chrome <chrome binary> --duckdb-version <vX.Y.Z>
//   node driver.mjs --native <duckdb shell with miint built in>
//
// --wasm boots duckdb-wasm in headless Chrome the way the playground does
// (index.html), serving the page and the extension with no COOP/COEP headers,
// like GitHub Pages, so duckdb-wasm picks the same bundle the playground gets.
// --duckdb-version is the core the extension targets; the run fails if the
// pinned duckdb-wasm reports another.
// --native runs the same steps in a native shell. A new step's `expect` comes
// from there. Run it again after changing code a step exercises: a value that
// changes is a finding to raise, not a value to update (CLAUDE.md: expected
// values change only with a maintainer's permission), and the SQL test under
// test/sql/ that owns the behaviour is where to look first.
//
// A step is {name, sql, expect?, expect_error?}, and names are unique. Steps run
// in order on one connection. A step with `expect` must return it as the first
// cell of its first row; cast the result to VARCHAR in the SQL so both runtimes
// print it the same way. A step with `expect_error` must fail with an error
// containing that text. A step with neither (setup) must only succeed. After
// every step, in both modes, the database must still answer SELECT 'alive'; a
// step after which it does not fails, and the steps after it are not run. The
// run passes only when every step ran and passed.
import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import http from 'node:http';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const steps = JSON.parse(fs.readFileSync(path.join(here, 'steps.json'), 'utf8'));
const args = Object.fromEntries(process.argv.slice(2).reduce((acc, a, i, all) =>
  (a.startsWith('--') ? [...acc, [a.slice(2), all[i + 1]]] : acc), []));

const names = new Set();
for (const s of steps) {
  if (typeof s.name !== 'string' || typeof s.sql !== 'string') {
    console.error(`steps.json: every step needs a string name and sql: ${JSON.stringify(s).slice(0, 120)}`);
    process.exit(2);
  }
  if (!/^[A-Za-z0-9_.-]+$/.test(s.name)) {
    console.error(`steps.json: step name '${s.name}' may use only letters, digits, '_', '.' and '-'`);
    process.exit(2);
  }
  if (s.expect !== undefined && s.expect_error !== undefined) {
    console.error(`steps.json: step '${s.name}' has both expect and expect_error`);
    process.exit(2);
  }
  if (s.expect_error !== undefined && (typeof s.expect_error !== 'string' || s.expect_error === '')) {
    console.error(`steps.json: step '${s.name}' needs a non-empty string expect_error`);
    process.exit(2);
  }
  if (names.has(s.name)) {
    console.error(`steps.json: duplicate step name '${s.name}'`);
    process.exit(2);
  }
  names.add(s.name);
}

// Each step is reported once, as passed or failed; whatever is left was not run.
// `aborted` is a failure outside any step (the shell or the page), which fails the
// run even if no step did.
let passed = 0;
let failed = 0;
let aborted = null;
// A line DuckDB writes for an error: "<Kind> Error: ..." or the shell's fixed
// message for a statement with unbound parameters.
const DUCKDB_ERROR_LINE = /^([A-Za-z][A-Za-z ]*Error: |Prepared statement parameters cannot be used directly)/;

// r is {ok, value, varchar, shape, alive} or {ok: false, error, alive}. varchar
// says the step returned one VARCHAR column; alive says whether the database
// answered SELECT 'alive' after the step.
function report(step, r) {
  let why = null;
  if (step.expect_error !== undefined) {
    const errors = r.ok ? [] : r.error.split('\n').filter(l => DUCKDB_ERROR_LINE.test(l));
    if (r.ok) why = `expected an error containing "${step.expect_error}", but the step succeeded`;
    else if (!errors.some(l => l.includes(step.expect_error))) {
      why = `expected a DuckDB error containing "${step.expect_error}", got: ${r.error}`;
    }
  } else if (!r.ok) {
    why = r.error;
  } else if (step.expect !== undefined && !r.varchar) {
    why = `an expect step must return one VARCHAR column (got ${r.shape}); cast it to VARCHAR`;
  } else if (step.expect !== undefined && r.value !== step.expect) {
    why = `\n  expected ${step.expect}\n  got      ${r.value}`;
  }
  if (!r.alive) {
    const dead = 'the database no longer answers';
    if (why === null) why = r.ok ? `${dead} after this step` : `failed as expected, but ${dead}: ${r.error}`;
    else why += `\n  and ${dead}`;
  }
  if (why === null) {
    passed++;
    console.log(`ok   ${step.name}`);
  } else {
    failed++;
    console.log(`FAIL ${step.name}: ${why}`.replace(': \n', ':\n'));
  }
}

const NATIVE_TIMEOUT_MS = 900_000;
const NATIVE_MAX_BUFFER = 256 * 1024 * 1024;

// A step's result from its JSON output: one document per statement that returns
// rows, each starting on a line that begins with '['. The first such statement's
// first row is the step's result, as duckdb-wasm reports it. Throws if the
// output is not JSON. JSON null prints as NULL, as in the page. A one-column
// row whose value is a string or null is what a VARCHAR column prints; numbers
// print bare, and JSON.parse would round a BIGINT above 2^53.
function nativeResult(lines) {
  const docs = [];
  for (const line of lines.filter(l => l !== '')) {
    if (line.startsWith('[') || docs.length === 0) docs.push([line]);
    else docs.at(-1).push(line);
  }
  const parsed = docs.map(d => JSON.parse(d.join('\n')));
  const row = parsed[0]?.[0];
  if (row === undefined) return { value: null, varchar: true, shape: 'no rows' };
  const keys = Object.keys(row);
  const v = row[keys[0]];
  return { value: v === null ? 'NULL' : String(v), varchar: keys.length === 1 && (v === null || typeof v === 'string'),
           shape: `${keys.length} column(s), the first a JSON ${v === null ? 'null' : typeof v}` };
}

function runNative(duckdb) {
  // One shell session, with no ~/.duckdbrc (-no-init). Results print as JSON
  // (.mode json), so a step whose output has a line that does not start with
  // '[' or '{' failed, however the line is worded. All output goes to stderr
  // (.output stderr), so errors stay in order with the markers. After each step
  // the shell runs SELECT 'alive'. No -bail: a step expected to fail must not
  // end the run.
  //
  // -no-agent: DuckDB 2.0's shell otherwise prints errors as JSON when an AI
  // coding tool's environment variable (AI_AGENT, CLAUDECODE, ...) is set, and
  // otherwise, with stdin not a terminal, ends with a hint line outside any step.
  const script = ['.output stderr', '.mode json',
    ...steps.flatMap(s => [`.print @@step ${s.name}`, s.sql.trim().replace(/;?$/, ';'),
                           '.print @@alive', "SELECT 'alive' AS alive;"]),
    '.print @@end', ''].join('\n');
  const r = spawnSync(duckdb, ['-no-agent', '-no-init', ':memory:'], { input: script, encoding: 'utf8',
                                              timeout: NATIVE_TIMEOUT_MS, maxBuffer: NATIVE_MAX_BUFFER });
  if (r.error && r.error.code !== 'ETIMEDOUT' && r.error.code !== 'ENOBUFS') {
    aborted = `could not run ${duckdb}: ${r.error.message}`;
    return;
  }
  const sections = new Map();
  const outside = [];
  let cur = null;
  let ended = false;
  // In .mode json, DuckDB 2.0's shell prints a result's opening '[' before it
  // runs the statement, so an error follows it on the same line ("[Invalid
  // Input Error: ..."). That '[' is dropped, leaving the shell's error line.
  const unbracket = l => (l.startsWith('[') && DUCKDB_ERROR_LINE.test(l.slice(1)) ? l.slice(1) : l);
  for (const line of (r.stderr ?? '').split('\n').map(unbracket)) {
    const m = line.match(/^@@step (.*)$/);
    if (m) {
      cur = { out: [], alive: [], aliveSeen: false };
      sections.set(m[1], cur);
    } else if (line === '@@alive' && cur) {
      cur.aliveSeen = true;
    } else if (line === '@@end') {
      ended = true;
      cur = null;
    } else if (cur) {
      (cur.aliveSeen ? cur.alive : cur.out).push(line);
    } else if (line.trim() !== '') {
      outside.push(line);
    }
  }
  const stopped = r.error?.code === 'ETIMEDOUT' ? `the shell timed out after ${NATIVE_TIMEOUT_MS / 1000} s`
    : r.error?.code === 'ENOBUFS' ? `the shell printed more than ${NATIVE_MAX_BUFFER} bytes`
    : `the shell stopped (${r.signal ? `signal ${r.signal}` : `status ${r.status}`})`;
  if (sections.size === 0) {
    aborted = `${stopped} before the first step: ${(r.stderr ?? '').trim() || '(no output)'}`;
    return;
  }
  const isJson = l => l === '' || l.startsWith('[') || l.startsWith('{');
  const answered = lines => {
    try { return JSON.stringify(JSON.parse(lines.join('\n'))) === '[{"alive":"alive"}]'; } catch { return false; }
  };
  let sawError = false;
  for (const s of steps) {
    const sec = sections.get(s.name);
    if (!sec) break;
    const at = sec.out.findIndex(l => !isJson(l));
    const alive = sec.aliveSeen && answered(sec.alive);
    if (at >= 0 || (sec.aliveSeen && !alive)) sawError = true;
    if (at >= 0) report(s, { ok: false, error: sec.out.slice(at).join('\n').trim(), alive });
    else if (!sec.aliveSeen) report(s, { ok: false, error: `${stopped} during this step`, alive: false });
    else {
      let result;
      try {
        result = nativeResult(sec.out);
      } catch {
        result = null;
      }
      if (result) report(s, { ok: true, ...result, alive });
      else report(s, { ok: false, error: `output is not JSON:\n${sec.out.join('\n').trim()}`, alive });
    }
    if (!alive) break;
  }
  if (outside.length > 0) {
    aborted = `the shell printed outside any step: ${outside.slice(0, 5).map(l => JSON.stringify(l)).join(', ')}`;
  } else if (ended && (r.signal || (r.status !== 0 && !(r.status === 1 && sawError)))) {
    // The shell exits 1 after any error, expected ones included. Any other
    // nonzero status, a signal, or 1 with no error shown is unexplained.
    aborted = `${stopped}, and no step's output explains it`;
  }
}

async function runWasm(wasmFile, chrome, duckdbVersion) {
  const { default: puppeteer } = await import('puppeteer-core');
  const files = { '/': [path.join(here, 'index.html'), 'text/html'],
                  '/miint.duckdb_extension.wasm': [wasmFile, 'application/wasm'] };
  const server = http.createServer((req, res) => {
    const f = files[new URL(req.url, 'http://localhost').pathname];
    if (!f) { res.writeHead(404); return res.end(); }
    res.writeHead(200, { 'content-type': f[1], 'content-length': fs.statSync(f[0]).size, 'cache-control': 'no-store' });
    if (req.method === 'HEAD') return res.end();
    fs.createReadStream(f[0]).pipe(res);
  }).listen(0, '127.0.0.1');
  await new Promise(r => server.once('listening', r));
  const base = `http://127.0.0.1:${server.address().port}`;

  const browser = await puppeteer.launch({ executablePath: chrome, headless: true, protocolTimeout: 900_000 });
  try {
    const page = await browser.newPage();
    const log = [];
    page.on('workercreated', w => w.on('console', m => log.push(m.text())));
    await page.goto(`${base}/`, { waitUntil: 'domcontentloaded' });
    await page.waitForFunction(() => window.__ready === true, { timeout: 120_000 });
    const boot = await page.evaluate(v => window.__boot(v), duckdbVersion);
    if (!boot.ok) {
      aborted = `boot: ${boot.error}`;
      return;
    }
    console.log(`duckdb-wasm ${boot.version}, platform ${boot.platform}, bundle ${boot.mainModule}`);
    const load = await page.evaluate(u => window.__load(u), `${base}/miint.duckdb_extension.wasm`);
    if (!load.ok) {
      aborted = `LOAD: ${load.error}`;
      return;
    }
    console.log(`miint ${load.version}`);
    for (const s of steps) {
      const before = log.length;
      const r = await page.evaluate(q => window.__step(q), s.sql);
      report(s, r);
      if ((!r.ok && s.expect_error === undefined) || !r.alive) {
        const worker = log.slice(before).filter(l => l.trim());
        if (worker.length) console.log(worker.slice(0, 12).map(l => `  [worker] ${l}`).join('\n'));
      }
      if (!r.alive) break;
    }
  } finally {
    await browser.close();
    server.close();
  }
}

if (args.native) runNative(args.native);
else if (args.wasm && args.chrome && args['duckdb-version']) await runWasm(args.wasm, args.chrome, args['duckdb-version']);
else {
  console.error('usage: driver.mjs --wasm <file> --chrome <binary> --duckdb-version <vX.Y.Z> | --native <duckdb>');
  process.exit(2);
}
const notRun = steps.length - passed - failed;
if (aborted) console.log(`FAIL ${aborted}`);
if (!aborted && failed === 0 && notRun === 0) {
  console.log(`all ${steps.length} steps passed`);
  process.exit(0);
}
console.log(`${passed} passed, ${failed} failed, ${notRun} not run, of ${steps.length} steps`);
process.exit(1);
