# Contributing to @ratmath/linalg

```sh
git submodule update --init   # @ratmath/algebraic, until it is on npm
bun install
bun test
bun run docs:build
```

## Docs site

`docs/` is a static site: an introduction, a reference, and topic pages
(elimination, the four subspaces, determinants, inverses, eigen, Jordan
forms, orthogonality, factorizations, integer matrices, exp, nonlinear solving,
other number systems) with editable
demos showing preview and exact side by side, plus step-by-step players. The
exact workbench page is a free-form scratchpad: named matrices that take any
Core rational literal (`3/4`, `0.#3`, `2..1/2`, `3.~7~15`), an expression line
(`A^-1`, `C = rref([A, b])`, `B'B`) and buttons for every exact operation. A Numbers setting adds real algebraic
numbers, complex ones, ranges (intervals) or arithmetic mod p.
Each page is its own HTML file (`index.html`, `eigen.html`, `workbench.html`, …),
generated from `docs/src/layout.html` (head and navigation) around
`docs/src/pages/<id>.html`; the styles are `docs/assets/site.css`. Edit the
sources, then run `bun run docs:build`, which bundles the library and writes
the pages, and serve `docs/`. Old links like `index.html#eigen` redirect.
Every page ends with an editable JavaScript box that runs with the library in
scope (seed code in the page's `<template data-runbox>`); the build runs each
box and fails if one throws (`--verbose` prints what each shows). Asset links
carry a content-hash `?v=` so browsers never mix cached scripts with new pages,
and math written as `\( … \)` / `\[ … \]` is typeset with KaTeX at build time
against the shared stylesheet at `ratmath.com/assets/katex/`.
`docs/assets/runbox.js` and `scripts/docs-kit.js` are shared verbatim with the
core, analytic and algebraic docs; keep the copies in sync.

