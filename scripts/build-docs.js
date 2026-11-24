// Bundle the library (with @ratmath/core) for the docs site, assemble one HTML
// file per page from docs/src (layout.html around pages/<id>.html), and
// optionally write a preview copy (index.html without its <html> wrapper).
import { $ } from "bun";
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { renderMath, addRunbox, stampAssets, stampLinks, KATEX_LINK } from "./docs-kit.js";
import { runCode } from "../docs/assets/runbox.js";

const out = await Bun.build({ entrypoints: ["./scripts/docs-entry.js"], outdir: "./docs/assets", naming: "linalg.js", target: "browser", format: "esm", minify: { whitespace: true, syntax: true, identifiers: true, keepNames: true } });
if (!out.success) {
  for (const m of out.logs) console.error(m);
  process.exit(1);
}
console.log("docs/assets/linalg.js written");

// Asset URLs carry a hash of the assets so a browser never mixes a cached old
// script with a new page. Math in \( … \) / \[ … \] is typeset with KaTeX here,
// and every page ends with the runnable code box (docs/assets/runbox.js, shared
// with the core, analytic and algebraic docs), seeded from <template data-runbox>.
const v = stampAssets("docs/assets");
const seeds = [];

// Pages, in navigation order (the order of the layout's <select>).
const layout = readFileSync("docs/src/layout.html", "utf8");
const PAGES = [...layout.matchAll(/<option value="([a-z-]+)">([^<]+)<\/option>/g)].map(([, id, name]) => ({ id, name }));
const fileOf = (id) => (id === "intro" ? "index.html" : `${id}.html`);
const files = {};
PAGES.forEach(({ id, name }, i) => {
  const prev = PAGES[i - 1];
  const next = PAGES[i + 1];
  const pager = `<nav class="pager" aria-label="Previous and next page">${prev ? `<a class="prev" href="${fileOf(prev.id)}"><small>Previous</small>${prev.name}</a>` : "<span></span>"}${next ? `<a class="next" href="${fileOf(next.id)}"><small>Next</small>${next.name}</a>` : "<span></span>"}</nav>`;
  const boxed = addRunbox(readFileSync(`docs/src/pages/${id}.html`, "utf8"), { lib: "<code>@ratmath/linalg</code> and Core's <code>Rational</code>, <code>RationalInterval</code>" });
  seeds.push({ page: id, seed: boxed.seed });
  const math = renderMath(boxed.html);
  const content = math.html.trimEnd() + "\n" + pager;
  const html = layout
    .replaceAll("{{TITLE}}", id === "intro" ? "RatMath Linear Algebra" : `${name} · RatMath linalg`)
    .replace("{{FILE}}", id === "intro" ? "" : fileOf(id))
    .replace("{{KATEX}}", math.count ? KATEX_LINK : "")
    .replace("{{CONTENT}}", () => content)
    .replace(`<option value="${id}">`, `<option value="${id}" selected>`)
    .replace(`<a href="${fileOf(id)}">`, `<a href="${fileOf(id)}" class="on" aria-current="page">`);
  writeFileSync(`docs/${fileOf(id)}`, stampLinks(html, v));
  files[fileOf(id)] = stampLinks(html, v);
});
console.log(`${PAGES.length} pages written to docs/ (assets ?v=${v})`);

// Every page's code box must run.
const scope = await import("./docs-entry.js");
let failed = 0;
for (const { page, seed } of seeds) {
  const res = await runCode(seed, scope);
  if (res.error) {
    failed++;
    console.error(`✗ ${page}: ${res.error.name}: ${res.error.message}`);
  } else if (process.argv.includes("--verbose")) {
    console.log(`— ${page} (${res.ms.toFixed(0)} ms)`);
    for (const pr of res.prints) console.log(`  ${pr.text}`);
    if (res.hasValue) console.log(`  → ${res.type ? `${res.type} ` : ""}${res.text}`);
  }
}
if (failed) {
  console.error(`${failed} code box(es) failed`);
  process.exit(1);
}
console.log(`${seeds.length} code boxes run without errors`);
if (process.argv.includes("--preview")) {
  const dir = process.argv[process.argv.indexOf("--preview") + 1];
  mkdirSync(`${dir}/assets`, { recursive: true });
  const body = files["index.html"]
    .replace(/<!doctype html>/i, "")
    .replace(/<html[^>]*>|<\/html>|<head>|<\/head>|<body>|<\/body>/gi, "")
    .replace(/<meta charset[^>]*>/i, "")
    .replace(/<meta name="viewport"[^>]*>/i, "");
  writeFileSync(`${dir}/index.html`, body);
  for (const [f, html] of Object.entries(files)) if (f !== "index.html") writeFileSync(`${dir}/${f}`, html);
  for (const f of ["linalg.js", "explorer.js", "pages.js", "workbench.js", "workbench-eval.js", "linalg-worker.js", "runbox.js", "site.css"]) writeFileSync(`${dir}/assets/${f}`, readFileSync(`docs/assets/${f}`));
  console.log(`preview written to ${dir}`);
}
void $;
