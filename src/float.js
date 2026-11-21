// Quick previews in Float64. These are what a screen shows while a slider is
// moving: microseconds for small matrices, no allocation of BigInts. Each
// result carries `trust`, a cheap judgement of how far to believe it:
//   "good"   well conditioned; digits shown are very likely right
//   "rough"  ill conditioned; the shape is right but trailing digits are not
//   "unsure" a pivot or rank decision fell inside rounding noise; the exact
//            pass may change the kind of answer (singular vs. not, rank)

const EPS = 2 ** -52;

function copy(a, rows, cols) {
  return Array.from({ length: rows }, (_, i) => Array.from(a.subarray(i * cols, (i + 1) * cols)));
}

function maxAbs(a) {
  let m = 0;
  for (const v of a) if (Math.abs(v) > m) m = Math.abs(v);
  return m;
}

/** Row reduction with partial pivoting and a rank tolerance. */
export function rrefFloat(flat, rows, cols, columns = cols) {
  const a = copy(flat, rows, cols);
  const scale = maxAbs(flat) || 1;
  const tol = Math.max(rows, cols) * EPS * scale * 16;
  const pivots = [];
  let unsure = false;
  let r = 0;
  for (let col = 0; col < columns && r < rows; col += 1) {
    let best = r;
    for (let i = r + 1; i < rows; i += 1) if (Math.abs(a[i][col]) > Math.abs(a[best][col])) best = i;
    const size = Math.abs(a[best][col]);
    if (size <= tol) {
      if (size > 0) unsure = true;
      for (let i = r; i < rows; i += 1) a[i][col] = 0;
      continue;
    }
    if (size < tol * 1e6) unsure = true;
    [a[r], a[best]] = [a[best], a[r]];
    const p = a[r][col];
    for (let j = 0; j < cols; j += 1) a[r][j] /= p;
    for (let i = 0; i < rows; i += 1) {
      if (i === r) continue;
      const f = a[i][col];
      if (f === 0) continue;
      for (let j = 0; j < cols; j += 1) a[i][j] -= f * a[r][j];
      a[i][col] = 0;
    }
    pivots.push(col);
    r += 1;
  }
  // tidy noise-level entries so the preview does not show 1e-17
  for (const row of a) for (let j = 0; j < cols; j += 1) if (Math.abs(row[j]) <= tol) row[j] = 0;
  return { rows: a, pivots, unsure };
}

/** LU with partial pivoting, in place on a copy. */
function lu(flat, n) {
  const a = copy(flat, n, n);
  const perm = Array.from({ length: n }, (_, i) => i);
  let sign = 1;
  let singular = false;
  const scale = maxAbs(flat) || 1;
  const tol = n * EPS * scale * 16;
  let minPivot = Infinity;
  for (let k = 0; k < n; k += 1) {
    let best = k;
    for (let i = k + 1; i < n; i += 1) if (Math.abs(a[i][k]) > Math.abs(a[best][k])) best = i;
    if (best !== k) {
      [a[k], a[best]] = [a[best], a[k]];
      [perm[k], perm[best]] = [perm[best], perm[k]];
      sign = -sign;
    }
    const p = a[k][k];
    minPivot = Math.min(minPivot, Math.abs(p));
    if (Math.abs(p) <= tol) {
      singular = true;
      continue;
    }
    for (let i = k + 1; i < n; i += 1) {
      const f = (a[i][k] /= p);
      if (f === 0) continue;
      for (let j = k + 1; j < n; j += 1) a[i][j] -= f * a[k][j];
    }
  }
  return { a, perm, sign, singular, minPivot, tol };
}

export function detFloat(flat, n) {
  if (n === 0) return { value: 1, trust: "good" };
  const { a, sign, singular, minPivot, tol } = lu(flat, n);
  let d = sign;
  for (let k = 0; k < n; k += 1) d *= a[k][k];
  if (singular) return { value: 0, trust: minPivot === 0 ? "rough" : "unsure" };
  return { value: d, trust: minPivot < tol * 1e6 ? "unsure" : "good" };
}

function luSolve({ a, perm }, b, n) {
  // b: rows × k (array of arrays)
  const k = b[0]?.length ?? 0;
  const x = perm.map((p) => [...b[p]]);
  for (let i = 0; i < n; i += 1) for (let j = 0; j < i; j += 1) {
    const f = a[i][j];
    if (f) for (let c = 0; c < k; c += 1) x[i][c] -= f * x[j][c];
  }
  for (let i = n - 1; i >= 0; i -= 1) {
    for (let j = i + 1; j < n; j += 1) {
      const f = a[i][j];
      if (f) for (let c = 0; c < k; c += 1) x[i][c] -= f * x[j][c];
    }
    for (let c = 0; c < k; c += 1) x[i][c] /= a[i][i];
  }
  return x;
}

function norm1(rows) {
  let m = 0;
  const cols = rows[0]?.length ?? 0;
  for (let j = 0; j < cols; j += 1) {
    let s = 0;
    for (const row of rows) s += Math.abs(row[j]);
    m = Math.max(m, s);
  }
  return m;
}

function trustFromCondition(cond) {
  if (!Number.isFinite(cond) || cond * EPS > 1e-3) return "unsure";
  if (cond * EPS > 1e-10) return "rough";
  return "good";
}

export function inverseFloat(flat, n) {
  const f = lu(flat, n);
  if (f.singular) return { value: null, trust: f.minPivot === 0 ? "rough" : "unsure", condition: Infinity };
  const id = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  const inv = luSolve(f, id, n);
  const condition = norm1(copy(flat, n, n)) * norm1(inv);
  return { value: inv, trust: trustFromCondition(condition), condition };
}

/**
 * Solve A x = B. Square systems that look nonsingular go through LU; anything
 * else through rref of [A | B]. Same result shape as the exact solver.
 */
export function solveFloat(aFlat, rows, cols, bRows) {
  const k = bRows[0]?.length ?? 0;
  if (rows === cols && rows > 0) {
    const f = lu(aFlat, rows);
    if (!f.singular) {
      const x = luSolve(f, bRows, rows);
      // condition estimate from the factors: cheap ‖A‖₁ · ‖A⁻¹‖₁ only for small n
      let trust = f.minPivot < f.tol * 1e6 ? "unsure" : "good";
      if (rows <= 60) {
        const id = Array.from({ length: rows }, (_, i) => Array.from({ length: rows }, (_, j) => (i === j ? 1 : 0)));
        const cond = norm1(copy(aFlat, rows, cols)) * norm1(luSolve(f, id, rows));
        const t = trustFromCondition(cond);
        if (t !== "good") trust = t;
      }
      return { value: { status: "unique", x, basis: [] }, trust };
    }
  }
  const width = cols + k;
  const wide = new Float64Array(rows * width);
  for (let i = 0; i < rows; i += 1) {
    for (let j = 0; j < cols; j += 1) wide[i * width + j] = aFlat[i * cols + j];
    for (let j = 0; j < k; j += 1) wide[i * width + cols + j] = bRows[i][j];
  }
  const { rows: r, pivots, unsure } = rrefFloat(wide, rows, width, cols);
  const trust = unsure ? "unsure" : "good";
  for (let i = pivots.length; i < rows; i += 1) {
    if (r[i].slice(cols).some((v) => v !== 0)) return { value: { status: "none", x: null, basis: basisOf(r, pivots, cols) }, trust };
  }
  const x = Array.from({ length: cols }, () => new Array(k).fill(0));
  pivots.forEach((col, i) => {
    for (let j = 0; j < k; j += 1) x[col][j] = r[i][cols + j];
  });
  const basis = basisOf(r, pivots, cols);
  return { value: { status: basis.length ? "infinite" : "unique", x, basis }, trust };
}

function basisOf(r, pivots, n) {
  const isPivot = new Set(pivots);
  const basis = [];
  for (let f = 0; f < n; f += 1) {
    if (isPivot.has(f)) continue;
    const v = new Array(n).fill(0);
    v[f] = 1;
    pivots.forEach((col, i) => (v[col] = -r[i][f]));
    basis.push(v);
  }
  return basis;
}

/**
 * Best rational approximation p/q of x with q ≤ maxDen (continued fractions).
 * Used to guess an exact answer from a float preview; the guess is then
 * checked exactly, never trusted.
 */
export function guessRational(x, maxDen = 1e6) {
  if (!Number.isFinite(x)) return null;
  const sign = x < 0 ? -1n : 1n;
  let v = Math.abs(x);
  let [h0, h1, k0, k1] = [0n, 1n, 1n, 0n];
  for (let i = 0; i < 40; i += 1) {
    const a = Math.floor(v);
    if (!Number.isSafeInteger(a)) return null;
    const A = BigInt(a);
    const h2 = A * h1 + h0;
    const k2 = A * k1 + k0;
    if (k2 > BigInt(maxDen)) break;
    [h0, h1, k0, k1] = [h1, h2, k1, k2];
    const frac = v - a;
    if (frac < 1e-12 * Math.max(1, Math.abs(x))) break;
    v = 1 / frac;
  }
  return k1 === 0n ? null : [sign * h1, k1];
}
