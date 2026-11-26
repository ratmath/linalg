// Build docs/<name>.pdf from the built docs site: open every page in headless
// Chromium (so demos mount and each "Try it" box runs), take the rendered
// <main>, and print all pages as one document with a cover, a contents page
// and a PDF outline. Run after docs:build.
import { chromium } from "playwright";
import { readFileSync, writeFileSync } from "node:fs";
import { createServer } from "node:http";
import { join, extname } from "node:path";

const pkg = JSON.parse(readFileSync("package.json", "utf8"));
const short = pkg.name.replace("@ratmath/", "");
const layout = readFileSync("docs/src/layout.html", "utf8");
const PAGES = [...layout.matchAll(/<option value="([a-z-]+)">([^<]+)<\/option>/g)].map(([, id, name]) => ({ id, name, file: id === "intro" ? "index.html" : `${id}.html` }));
const SKIP = new Set((process.env.PDF_SKIP ?? "workbench").split(",")); // pages that are tools, not text

const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".svg": "image/svg+xml", ".png": "image/png" };
const server = createServer((req, res) => {
  const p = join("docs", decodeURIComponent(new URL(req.url, "http://x").pathname));
  try { res.writeHead(200, { "content-type": TYPES[extname(p)] ?? "application/octet-stream" }); res.end(readFileSync(p)); }
  catch { res.writeHead(404); res.end(); }
}).listen(0);
const base = `http://localhost:${server.address().port}/`;

const browser = await chromium.launch();
const page = await browser.newPage({ colorScheme: "light", viewport: { width: 900, height: 1200 } });
const parts = [];
for (const p of PAGES) {
  if (SKIP.has(p.id)) continue;
  await page.goto(base + p.file, { waitUntil: "networkidle" });
  await page.evaluate(async () => {
    for (const b of document.querySelectorAll(".runbox [data-run]")) b.click();
    await new Promise((r) => setTimeout(r, 300));
    // Freeze form state into attributes so it survives outerHTML.
    for (const i of document.querySelectorAll("input")) i.setAttribute("value", i.value);
    for (const t of document.querySelectorAll("textarea")) t.textContent = t.value;
    for (const s of document.querySelectorAll("select option")) s.toggleAttribute("selected", s.selected);
  });
  await page.waitForTimeout(200);
  const html = await page.evaluate(() => {
    const m = document.querySelector("main").cloneNode(true);
    m.querySelectorAll(".pager, .rb-try, .runbox-bar button, script").forEach((e) => e.remove());
    for (const t of m.querySelectorAll("textarea")) { const pre = document.createElement("pre"); pre.className = "rb-code"; pre.textContent = t.textContent; t.replaceWith(pre); }
    return m.innerHTML;
  });
  const styles = await page.evaluate(() => [...document.querySelectorAll("style")].map((s) => s.textContent).join("\n"));
  parts.push({ ...p, html, styles });
  console.log(`rendered ${p.file}`);
}

const head = readFileSync("docs/index.html", "utf8").match(/<head>([\s\S]*?)<\/head>/)[1].replace(/<title>.*<\/title>/, `<title>RatMath ${short}</title>`);
const date = new Date().toISOString().slice(0, 10);
const css = `
@page { size: A4; margin: 18mm 16mm 20mm; }
body { background: #fff; font-size: 10.5pt; }
.pdf-page { break-before: page; }
.pdf-page > section.page { margin: 0; }
pre, .runbox, figure, table, [data-demo] { break-inside: avoid; }
h1, h2, h3, h4 { break-after: avoid; }
pre, code { white-space: pre-wrap; overflow-wrap: anywhere; }
.runbox textarea { height: auto; overflow: visible; resize: none; }
.cover { height: 240mm; display: flex; flex-direction: column; justify-content: center; }
.cover h1 { font-size: 34pt; margin: 0 0 .3em; }
.toc-list a { color: inherit; text-decoration: none; }
.toc-list li { margin: .25em 0; }
`;
const doc = `<!doctype html><html lang="en"><head>${head}<base href="${base}"><style>${parts[0].styles}</style><style>${css}</style></head><body>
<section class="cover"><div class="eyebrow">RatMath Shelf</div><h1>${pkg.name}</h1><p class="lede">${pkg.description}</p><p class="muted">Version ${pkg.version} · ${date} · https://ratmath.com/${short}/</p></section>
<section class="pdf-page"><h1>Contents</h1><ol class="toc-list">${parts.map((p) => `<li><a href="#p-${p.id}">${p.name}</a></li>`).join("")}</ol></section>
${parts.map((p) => `<div class="pdf-page" id="p-${p.id}">${p.html}</div>`).join("\n")}
</body></html>`;
await page.setContent(doc, { waitUntil: "networkidle" });
await page.pdf({ path: `docs/${short}.pdf`, format: "A4", printBackground: true, outline: true, tagged: true, displayHeaderFooter: true,
  headerTemplate: "<span></span>",
  footerTemplate: `<div style="font:8px sans-serif;width:100%;padding:0 16mm;display:flex;justify-content:space-between;color:#666"><span>${pkg.name} ${pkg.version}</span><span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>` });
await browser.close();
server.close();
console.log(`docs/${short}.pdf written`);
