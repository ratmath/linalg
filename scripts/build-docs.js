// Bundle the library (with @ratmath/core) for the docs site, assemble one HTML
// file per page from docs/src (layout.html around pages/<id>.html), and
// optionally write a preview copy (index.html without its <html> wrapper).
import { $ } from "bun";
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";

const out = await Bun.build({ entrypoints: ["./scripts/docs-entry.js"], outdir: "./docs/assets", naming: "linalg.js", target: "browser", format: "esm", minify: true });
if (!out.success) {
  for (const m of out.logs) console.error(m);
  process.exit(1);
}
console.log("docs/assets/linalg.js written");

// Pages, in navigation order (the order of the layout's <select>).
const layout = readFileSync("docs/src/layout.html", "utf8");
const PAGES = [...layout.matchAll(/<option value="([a-z-]+)">([^<]+)<\/option>/g)].map(([, id, name]) => ({ id, name }));
const fileOf = (id) => (id === "intro" ? "index.html" : `${id}.html`);
const files = {};
PAGES.forEach(({ id, name }, i) => {
  const prev = PAGES[i - 1];
  const next = PAGES[i + 1];
  const pager = `<nav class="pager" aria-label="Previous and next page">${prev ? `<a class="prev" href="${fileOf(prev.id)}"><small>Previous</small>${prev.name}</a>` : "<span></span>"}${next ? `<a class="next" href="${fileOf(next.id)}"><small>Next</small>${next.name}</a>` : "<span></span>"}</nav>`;
  const content = readFileSync(`docs/src/pages/${id}.html`, "utf8").trimEnd() + "\n" + pager;
  const html = layout
    .replace("{{TITLE}}", id === "intro" ? "RatMath Linear Algebra" : `${name} · RatMath linalg`)
    .replace("{{CONTENT}}", content)
    .replace(`<option value="${id}">`, `<option value="${id}" selected>`)
    .replace(`<a href="${fileOf(id)}">`, `<a href="${fileOf(id)}" class="on" aria-current="page">`);
  writeFileSync(`docs/${fileOf(id)}`, html);
  files[fileOf(id)] = html;
});
console.log(`${PAGES.length} pages written to docs/`);
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
  for (const f of ["linalg.js", "explorer.js", "pages.js", "workbench.js", "workbench-eval.js", "linalg-worker.js", "site.css"]) writeFileSync(`${dir}/assets/${f}`, readFileSync(`docs/assets/${f}`));
  console.log(`preview written to ${dir}`);
}
void $;
