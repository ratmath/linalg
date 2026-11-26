// Build docs/llms.txt (a compact map for AI assistants, llmstxt.org format)
// and docs/llms-full.txt (every docs page as plain Markdown, with its runnable
// example). Generated from the source, so file:line pointers are always the
// ones in this commit. Run after docs:build.
import { readFileSync, writeFileSync, existsSync } from "node:fs";
import { execSync } from "node:child_process";

const pkg = JSON.parse(readFileSync("package.json", "utf8"));
const short = pkg.name.replace("@ratmath/", "");
const site = `https://ratmath.com/${short}/`;
let sha = "main";
try { sha = execSync("git rev-parse --short HEAD").toString().trim(); } catch {}
const gh = `https://github.com/ratmath/${short}/blob/${sha}/`;

// ---- docs pages, in nav order -------------------------------------------
const layout = readFileSync("docs/src/layout.html", "utf8");
const PAGES = [...layout.matchAll(/<option value="([a-z-]+)">([^<]+)<\/option>/g)].map(([, id, name]) => ({ id, name, file: id === "intro" ? "index.html" : `${id}.html` }));
const ent = (s) => s.replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&quot;/g, '"').replace(/&#39;/g, "'").replace(/&nbsp;/g, " ").replace(/&amp;/g, "&");
const src = (id) => readFileSync(`docs/src/pages/${id}.html`, "utf8");
const text = (h) => ent(h.replace(/<[^>]+>/g, "")).replace(/\s+/g, " ").trim();
for (const p of PAGES) {
  const h = src(p.id);
  p.lede = text(h.match(/<p[^>]*>([\s\S]*?)<\/p>/)?.[1] ?? "").split(/(?<=\.)\s/)[0];
  p.seed = ent(h.match(/<template data-runbox>([\s\S]*?)<\/template>/)?.[1] ?? "").trim();
}
function toMarkdown(h) {
  return ent(h
    .replace(/<template[\s\S]*?<\/template>/g, "")
    .replace(/<div[^>]*data-demo[^>]*><\/div>/g, "")
    .replace(/<pre[^>]*><code[^>]*>([\s\S]*?)<\/code><\/pre>/g, (_, c) => "\n```js\n" + c.trim() + "\n```\n")
    .replace(/<h([1-4])[^>]*>([\s\S]*?)<\/h\1>/g, (_, n, t) => `\n${"#".repeat(+n + 1)} ${text(t)}\n`)
    .replace(/<li[^>]*>([\s\S]*?)<\/li>/g, (_, t) => `- ${text(t)}\n`)
    .replace(/<tr[^>]*>([\s\S]*?)<\/tr>/g, (_, r) => "| " + [...r.matchAll(/<t[dh][^>]*>([\s\S]*?)<\/t[dh]>/g)].map((m) => text(m[1]).replace(/\|/g, "\\|")).join(" | ") + " |\n")
    .replace(/<code>([\s\S]*?)<\/code>/g, "`$1`")
    .replace(/<a [^>]*href="([^"]+)"[^>]*>\s*<strong>([\s\S]*?)<\/strong>([\s\S]*?)<\/a>/g, (_, u, t, d) => `- [${text(t)}](${/^https?:/.test(u) ? u : site + u}): ${text(d)}\n`)
    .replace(/<a [^>]*href="([^"]+)"[^>]*>([\s\S]*?)<\/a>/g, (_, u, t) => `[${text(t)}](${/^https?:/.test(u) ? u : site + u})`)
    .replace(/<\/p>|<br\s*\/?>/g, "\n\n")
    .replace(/<(?!\/?(?:pre|code))[^>]+>/g, ""))
    .split(/(```[\s\S]*?```)/).map((c, i) => (i % 2 ? c : c.replace(/^[ \t]+/gm, ""))).join("")
    .replace(/[ \t]+\n/g, "\n").replace(/\n{3,}/g, "\n\n").trim();
}

// ---- public API: every name exported from index.js, with its definition ---
const lineOf = (file, re) => {
  const lines = readFileSync(file, "utf8").split("\n");
  const i = lines.findIndex((l) => re.test(l));
  return i < 0 ? null : i + 1;
};
const docOf = (file, line) => { // the /** … */ or // comment right above a definition
  const lines = readFileSync(file, "utf8").split("\n");
  let i = line - 2, out = [];
  if (/\*\/\s*$/.test(lines[i] ?? "")) { while (i >= 0 && !/\/\*\*/.test(lines[i])) out.unshift(lines[i--]); out.unshift(lines[i] ?? ""); }
  else while (i >= 0 && /^\s*\/\//.test(lines[i])) out.unshift(lines[i--]);
  const t = out.join(" ").replace(/\/\*\*|\*\/|^\s*\*|\s\*\s|\/\//g, " ").replace(/@\w+[^.]*/g, "").replace(/\s+/g, " ").trim();
  const first = t.split(/(?<=\.)\s/)[0];
  return first.length > 160 ? first.slice(0, 157).replace(/\s\S*$/, "") + "…" : first;
};
const def = (name) => new RegExp(`^export\\s+(?:async\\s+)?(?:function\\*?|class|const|let)\\s+${name.replace("$", "\\$")}\\b`);
const pageFor = (name) => PAGES.find((p) => !["intro", "reference", "workbench"].includes(p.id) && new RegExp(`<code>[^<]*\\b${name}\\b`).test(src(p.id))) ?? PAGES.find((p) => new RegExp(`<code>[^<]*\\b${name}\\b`).test(src(p.id)));
const idx = readFileSync("index.js", "utf8");
const rows = [], spaces = [];
for (const m of idx.matchAll(/export\s+(\*\s+as\s+(\w+)|\{([^}]*)\})\s+from\s+"([^"]+)"/g)) {
  const [, , ns, names, from] = m;
  const file = from.startsWith(".") ? from.replace(/^\.\//, "") : null;
  if (ns) {
    const members = [...readFileSync(file, "utf8").matchAll(/^export\s+(?:async\s+)?(?:function\*?|class|const|let)\s+([\w$]+)/gm)].map((x) => x[1]);
    const re = [...readFileSync(file, "utf8").matchAll(/^export\s+\*\s+from\s+"([^"]+)"/gm)].map((x) => x[1]);
    spaces.push({ ns, file, members: [...members.map((n) => `${n}:${lineOf(file, def(n))}`), ...re.map((r) => `everything from ${r}`)] });
    continue;
  }
  for (const name of names.split(",").map((s) => s.trim()).filter(Boolean)) {
    if (!file) { rows.push({ name, where: from }); continue; }
    const line = lineOf(file, def(name));
    rows.push({ name, file, line, doc: line ? docOf(file, line) : "", page: pageFor(name) });
  }
}

// ---- signatures from index.d.ts -----------------------------------------
// One compact line per exported function or class: parameter list and return
// type, with Refinement<Q, E> written "Q ⇒ E" (quick value, then exact value).
const dts = existsSync("index.d.ts") ? readFileSync("index.d.ts", "utf8") : "";
const decls = {};
{
  const lines = dts.split("\n");
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(/^export (function|const|class|type|interface) (\w+)/);
    if (!m) continue;
    let text = lines[i], depth = 0, j = i;
    const count = (t) => { for (const c of t) depth += "({[<".includes(c) ? 1 : ")}]>".includes(c) ? -1 : 0; };
    count(text.replace(/=>/g, ""));
    while (depth > 0 && j + 1 < lines.length) { text += " " + lines[++j].trim(); count(lines[j].replace(/=>/g, "")); }
    decls[m[2]] ??= { kind: m[1], text: text.replace(/\s+/g, " "), line: i + 1 };
    if (m[1] === "class") {
      const body = text.slice(text.indexOf("{") + 1, text.lastIndexOf("}"));
      decls[m[2]].members = [...new Set([...body.matchAll(/(?:^\s*|;\s*)(?:static |readonly |get )*(\w+)\s*[(:?]/g)].map((x) => x[1]).filter((n) => n !== "constructor"))];
    }
  }
}
const refine = (t) => {
  let out = t, prev;
  do { prev = out; out = out.replace(/Refinement<((?:[^<>]|<[^<>]*>|<[^<>]*<[^<>]*>[^<>]*>)*?), ((?:[^<>]|<[^<>]*>|<[^<>]*<[^<>]*>[^<>]*>)*)>/, "($1 ⇒ $2)"); } while (out !== prev);
  return out;
};
const sigOf = (name) => {
  const d = decls[name];
  if (!d) return "";
  if (d.kind === "function") return refine(d.text.replace(/^export function \w+(<[^(]*>)?/, "").replace(/;$/, "")).replace(/\): /, ") → ");
  if (d.kind === "const") return d.text.replace(/^export const \w+: /, ": ").replace(/;$/, "");
  if (d.kind === "class") return ` class: ${d.members.join(", ")}`;
  return "";
};
const types = Object.entries(decls).filter(([n, d]) => d.kind === "type" && d.text.length < 160).map(([n, d]) => d.text.replace(/^export /, "").replace(/;$/, ""));

// ---- llms.txt --------------------------------------------------------------
const L = [];
L.push(`# ${pkg.name}`, "", `> ${pkg.description}. Part of the RatMath Shelf; numbers come from @ratmath/core (peer dependency). ESM: \`import * as ${short} from "${pkg.name}"\`.`, "");
L.push(`Built from ${pkg.version}, commit ${sha}. \`file:line\` is relative to the package root (also node_modules/${pkg.name}/); on GitHub: ${gh}<file>#L<line>. Docs pages are under ${site}. Full types: index.d.ts. \`(Q ⇒ E)\` is a Refinement: \`.quick.value\` is the float preview Q at once, \`.exact()\` the exact E.`, "");
if (types.length) L.push("Types: " + types.map((t) => `\`${t}\``).join("; "), "");
L.push("## Docs", "", `- ${site}llms-full.txt: all pages as Markdown, each with a runnable example (~${Math.round(13000 / 1000)}k tokens)`, `- ${site}${short}.pdf`, "- Pages: " + PAGES.filter((p) => p.id !== "workbench").map((p) => `${p.file.replace(".html", "")} (${p.name})`).join(", "), "");
L.push("## API", "", "name(params) → result — meaning · file:line · page", "");
for (const r of rows) {
  if (!r.file) { L.push(`- \`${r.name}\` — from ${r.where}`); continue; }
  L.push(`- \`${r.name}${sigOf(r.name)}\`${r.doc ? ` — ${r.doc}` : ""} · ${r.file}:${r.line}${r.page ? ` · ${r.page.file.replace(".html", "")}` : ""}`);
}
L.push("", "## Namespaces", "");
for (const s of spaces) L.push(`- \`${s.ns}\` (${s.file}): ${s.members.map((m) => m.replace(/:\d+$/, "")).join(", ")}`);
const ex = PAGES.find((p) => p.id === "intro")?.seed;
if (ex) L.push("", "## Example", "", "```js", `import { ${[...rows.filter((r) => r.file).map((r) => r.name), ...spaces.map((x) => x.ns)].filter((n) => new RegExp(`(?<![.\\w])${n}\\s*[.(]`).test(ex)).join(", ")} } from "${pkg.name}";`, ex.replace(/\bPRINT\(/g, "console.log("), "```");
L.push("", "## More", "", `- tests/: worked examples for every function`, "- https://ratmath.com/core/llms.txt (numbers), https://ratmath.com/shelf/llms.txt (all packages)");
writeFileSync("docs/llms.txt", L.join("\n") + "\n");

// ---- llms-full.txt ---------------------------------------------------------
const F = [`# ${pkg.name} ${pkg.version} — full documentation`, "", `> ${pkg.description}. Source: ${site}. Each page ends with a runnable example (PRINT prints).`, ""];
for (const p of PAGES) {
  if (p.id === "workbench") continue;
  F.push(`\n---\n\nPage: ${site}${p.file}\n`, toMarkdown(src(p.id)));
  if (p.seed) F.push("", "### Try it", "", "```js", p.seed, "```");
}
writeFileSync("docs/llms-full.txt", F.join("\n") + "\n");
const tok = (f) => Math.round(readFileSync(f, "utf8").length / 4);
console.log(`docs/llms.txt (~${tok("docs/llms.txt")} tokens), docs/llms-full.txt (~${tok("docs/llms-full.txt")} tokens) written`);
