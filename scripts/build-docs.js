// Bundle the library (with @ratmath/core) for the docs site, and write a
// copy of the site for previewing as a single page (no <html> wrapper).
import { $ } from "bun";
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";

const out = await Bun.build({ entrypoints: ["./index.js"], outdir: "./docs/assets", naming: "linalg.js", target: "browser", format: "esm", minify: true });
if (!out.success) {
  for (const m of out.logs) console.error(m);
  process.exit(1);
}
console.log("docs/assets/linalg.js written");
if (process.argv.includes("--preview")) {
  const dir = process.argv[process.argv.indexOf("--preview") + 1];
  mkdirSync(`${dir}/assets`, { recursive: true });
  const html = readFileSync("docs/index.html", "utf8");
  const body = html
    .replace(/<!doctype html>/i, "")
    .replace(/<html[^>]*>|<\/html>|<head>|<\/head>|<body>|<\/body>/gi, "")
    .replace(/<meta charset[^>]*>/i, "")
    .replace(/<meta name="viewport"[^>]*>/i, "");
  writeFileSync(`${dir}/index.html`, body);
  for (const f of ["linalg.js", "explorer.js", "pages.js"]) writeFileSync(`${dir}/assets/${f}`, readFileSync(`docs/assets/${f}`));
  console.log(`preview written to ${dir}`);
}
void $;
