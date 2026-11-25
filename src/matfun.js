// Matrix functions beyond exp: sin, cos, sinh, cosh, the principal square
// root, the principal logarithm and real powers. Every exact stage returns an
// interval for each entry that certainly contains the true value.
//
// The arithmetic is fixed point: a number is a BigInt counted in units of
// 2^-P, an interval is [lo, hi] of those, and every operation rounds outward.
//
//   sin, cos   Taylor series of A/2^s with a tail bound, then s double-angle
//              steps cos 2X = 2cos²X − I, sin 2X = 2 sinX cosX.
//   sqrt       a float guess (Denman–Beavers), then a Krawczyk test on
//              X² − A = 0 in n² unknowns proves one root in a small box, and
//              a Lyapunov test proves every eigenvalue of that root has
//              positive real part, so the root is the principal one.
//   log        square roots until the matrix is near I, the series of
//              log(I + X), then doubling back.
//   power      exp(t · log A).

const abs = (a) => (a < 0n ? -a : a);
const floorDiv = (a, b) => {
  const q = a / b;
  return (a % b !== 0n && (a < 0n) !== (b < 0n)) ? q - 1n : q;
};
const ceilDiv = (a, b) => -floorDiv(-a, b);

const ZERO = [0n, 0n];
const add = (x, y) => [x[0] + y[0], x[1] + y[1]];
const sub = (x, y) => [x[0] - y[1], x[1] - y[0]];
const neg = (x) => [-x[1], -x[0]];
const mag = (x) => (abs(x[0]) > abs(x[1]) ? abs(x[0]) : abs(x[1]));

function kernel(P) {
  const p = BigInt(P);
  const one = 1n << p;
  const down = (a) => a >> p; // floor
  const up = (a) => -((-a) >> p); // ceil
  const mul = (x, y) => {
    const a = x[0] * y[0];
    const b = x[0] * y[1];
    const c = x[1] * y[0];
    const d = x[1] * y[1];
    let lo = a;
    let hi = a;
    if (b < lo) lo = b; else if (b > hi) hi = b;
    if (c < lo) lo = c; else if (c > hi) hi = c;
    if (d < lo) lo = d; else if (d > hi) hi = d;
    return [down(lo), up(hi)];
  };
  // sum of products, rounded once
  const dot = (xs, ys) => {
    let lo = 0n;
    let hi = 0n;
    for (let k = 0; k < xs.length; k += 1) {
      const x = xs[k];
      const y = ys[k];
      if (x[0] === x[1]) {
        if (x[0] === 0n) continue;
        if (x[0] > 0n) { lo += x[0] * y[0]; hi += x[0] * y[1]; } else { lo += x[0] * y[1]; hi += x[0] * y[0]; }
        continue;
      }
      const a = x[0] * y[0];
      const b = x[0] * y[1];
      const c = x[1] * y[0];
      const d = x[1] * y[1];
      let l = a;
      let h = a;
      if (b < l) l = b; else if (b > h) h = b;
      if (c < l) l = c; else if (c > h) h = c;
      if (d < l) l = d; else if (d > h) h = d;
      lo += l;
      hi += h;
    }
    return [down(lo), up(hi)];
  };
  const matmul = (X, Y) => {
    const n = Y.length;
    const cols = Y[0].map((_, j) => Array.from({ length: n }, (_, k) => Y[k][j]));
    return X.map((row) => cols.map((col) => dot(row, col)));
  };
  const fromPair = ([n, d]) => [floorDiv(n << p, d), ceilDiv(n << p, d)];
  const identity = (n) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? [one, one] : ZERO)));
  /** Upper bound of the ∞-norm, in units of 2^-P. */
  const norm = (X) => {
    let best = 0n;
    for (const row of X) {
      let s = 0n;
      for (const x of row) s += mag(x);
      if (s > best) best = s;
    }
    return best;
  };
  const fromFloat = (x) => {
    // a double is a dyadic rational; this is exact when P is large enough, and a valid point either way
    if (!Number.isFinite(x)) throw new RangeError("matrix function: the float stage overflowed");
    const scaled = BigInt(Math.round(x * 2 ** 40));
    return P >= 40 ? scaled << BigInt(P - 40) : scaled >> BigInt(40 - P);
  };
  const toNumber = (a) => {
    const shift = Math.max(0, a.toString(2).length - 1000);
    return Number(a >> BigInt(shift)) * 2 ** (shift - P);
  };
  return { P, p, one, down, up, mul, dot, matmul, fromPair, identity, norm, fromFloat, toNumber };
}

const madd = (X, Y) => X.map((row, i) => row.map((x, j) => add(x, Y[i][j])));
const msub = (X, Y) => X.map((row, i) => row.map((x, j) => sub(x, Y[i][j])));
const mneg = (X) => X.map((row) => row.map(neg));
const widen = (X, t) => X.map((row) => row.map(([lo, hi]) => [lo - t, hi + t]));
const shiftRight = (X, s) => (s === 0 ? X : X.map((row) => row.map(([lo, hi]) => [lo >> BigInt(s), -((-hi) >> BigInt(s))])));
const shiftLeft = (X, s) => (s === 0 ? X : X.map((row) => row.map(([lo, hi]) => [lo << BigInt(s), hi << BigInt(s)])));
const divInt = (X, k) => X.map((row) => row.map(([lo, hi]) => [floorDiv(lo, k), ceilDiv(hi, k)]));
const widest = (X) => X.flat().reduce((w, [lo, hi]) => (hi - lo > w ? hi - lo : w), 0n);

/** How many halvings bring the norm (units of 2^-P) to at most 2^-e. */
function halvings(norm, P, e = 1) {
  let s = 0;
  while ((norm >> BigInt(s)) > (1n << BigInt(P - e))) s += 1;
  return s;
}

function* taylorPowers(K, B) {
  // yields the terms B^k/k! for k = 1, 2, … as interval matrices
  let term = K.identity(B.length);
  for (let k = 1; ; k += 1) {
    term = divInt(K.matmul(term, B), BigInt(k));
    yield term;
  }
}

/** Σ_{k>N} b^k/k! ≤ 2·b^(N+1)/(N+1)! when b ≤ 1/2. */
function tailBound(K, b, terms) {
  // 2·b^(terms+1)/(terms+1)!, rounded up, in units of 2^-P
  let t = 2n * K.one;
  for (let k = 1; k <= terms + 1; k += 1) t = ceilDiv(t * b, K.one * BigInt(k)) + 1n;
  return t;
}

function termsFor(K, b) {
  for (let n = 1; n < 2000; n += 1) if (tailBound(K, b, n) <= 2n) return n;
  return 2000;
}

/** exp of an interval matrix, as an interval matrix. Generator (yields progress). */
function* expIv(K, A) {
  const s = halvings(K.norm(A), K.P);
  const B = shiftRight(A, s);
  const b = K.norm(B);
  const N = termsFor(K, b);
  let sum = K.identity(A.length);
  let k = 0;
  for (const term of taylorPowers(K, B)) {
    k += 1;
    sum = madd(sum, term);
    if (k % 6 === 0) yield { step: "series", done: k, of: N + s };
    if (k >= N) break;
  }
  let E = widen(sum, tailBound(K, b, N));
  for (let i = 0; i < s; i += 1) {
    E = K.matmul(E, E);
    yield { step: "square", done: N + i + 1, of: N + s };
  }
  return E;
}

/** [cos A, sin A] of an interval matrix. Generator. */
function* cosSinIv(K, A) {
  const n = A.length;
  const s = halvings(K.norm(A), K.P);
  const B = shiftRight(A, s);
  const b = K.norm(B);
  const N = termsFor(K, b);
  let C = K.identity(n);
  let S = Array.from({ length: n }, () => Array.from({ length: n }, () => ZERO));
  let k = 0;
  for (const term of taylorPowers(K, B)) {
    k += 1;
    // k mod 4: 1 → +sin, 2 → −cos, 3 → −sin, 0 → +cos
    if (k % 4 === 1) S = madd(S, term);
    else if (k % 4 === 2) C = msub(C, term);
    else if (k % 4 === 3) S = msub(S, term);
    else C = madd(C, term);
    if (k % 6 === 0) yield { step: "series", done: k, of: N + s };
    if (k >= N) break;
  }
  const t = tailBound(K, b, N);
  C = widen(C, t);
  S = widen(S, t);
  const I = K.identity(n);
  for (let i = 0; i < s; i += 1) {
    const C2 = msub(shiftLeft(K.matmul(C, C), 1), I);
    S = shiftLeft(K.matmul(S, C), 1);
    C = C2;
    yield { step: "double", done: N + i + 1, of: N + s };
  }
  return [C, S];
}

// ---- floats ------------------------------------------------------------------------------

const fmul = (X, Y) => X.map((row) => Y[0].map((_, j) => row.reduce((acc, x, k) => acc + x * Y[k][j], 0)));
const fadd = (X, Y, c = 1) => X.map((row, i) => row.map((x, j) => x + c * Y[i][j]));
const feye = (n) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
const fnorm = (X) => Math.max(0, ...X.map((r) => r.reduce((s, x) => s + Math.abs(x), 0)));

/** Inverse by Gauss–Jordan with partial pivoting, or null. */
function finv(A) {
  const n = A.length;
  const M = A.map((row, i) => [...row, ...feye(n)[i]]);
  for (let c = 0; c < n; c += 1) {
    let best = c;
    for (let r = c + 1; r < n; r += 1) if (Math.abs(M[r][c]) > Math.abs(M[best][c])) best = r;
    if (!(Math.abs(M[best][c]) > 1e-300)) return null;
    [M[c], M[best]] = [M[best], M[c]];
    const piv = M[c][c];
    for (let j = c; j < 2 * n; j += 1) M[c][j] /= piv;
    for (let r = 0; r < n; r += 1) {
      if (r === c || M[r][c] === 0) continue;
      const f = M[r][c];
      for (let j = c; j < 2 * n; j += 1) M[r][j] -= f * M[c][j];
    }
  }
  return M.map((row) => row.slice(n));
}

/** Principal square root in floats (Denman–Beavers), or null when the iteration does not settle. */
export function sqrtmFloat(A) {
  const n = A.length;
  let Y = A.map((r) => [...r]);
  let Z = feye(n);
  for (let k = 0; k < 80; k += 1) {
    const Yi = finv(Y);
    const Zi = finv(Z);
    if (!Yi || !Zi) return null;
    const Yn = fadd(Y, Zi).map((r) => r.map((x) => x / 2));
    Z = fadd(Z, Yi).map((r) => r.map((x) => x / 2));
    const change = fnorm(fadd(Yn, Y, -1));
    Y = Yn;
    if (!Y.flat().every(Number.isFinite)) return null;
    if (change <= 1e-15 * Math.max(1, fnorm(Y))) break;
  }
  const res = fnorm(fadd(fmul(Y, Y), A, -1));
  return res <= 1e-8 * Math.max(1, fnorm(A)) ? Y : null;
}

/** Principal logarithm in floats: square roots until near I, then the series. */
export function logmFloat(A) {
  const n = A.length;
  let R = A;
  let k = 0;
  while (fnorm(fadd(R, feye(n), -1)) > 0.25) {
    R = sqrtmFloat(R);
    if (!R || k > 60) return null;
    k += 1;
  }
  const X = fadd(R, feye(n), -1);
  let term = feye(n);
  let sum = X.map((r) => r.map(() => 0));
  for (let j = 1; j <= 40; j += 1) {
    term = fmul(term, X);
    sum = fadd(sum, term, (j % 2 ? 1 : -1) / j);
  }
  return sum.map((r) => r.map((x) => x * 2 ** k));
}

export function expFloat(A) {
  const n = A.length;
  let s = 0;
  const norm = fnorm(A);
  while (norm / 2 ** s > 0.5) s += 1;
  const B = A.map((r) => r.map((x) => x / 2 ** s));
  let term = feye(n);
  let sum = feye(n);
  for (let k = 1; k <= 18; k += 1) {
    term = fmul(term, B).map((r) => r.map((x) => x / k));
    sum = fadd(sum, term);
  }
  for (let k = 0; k < s; k += 1) sum = fmul(sum, sum);
  return sum;
}

export function cosSinFloat(A) {
  const n = A.length;
  let s = 0;
  const norm = fnorm(A);
  while (norm / 2 ** s > 0.5) s += 1;
  const B = A.map((r) => r.map((x) => x / 2 ** s));
  let term = feye(n);
  let C = feye(n);
  let S = A.map((r) => r.map(() => 0));
  for (let k = 1; k <= 20; k += 1) {
    term = fmul(term, B).map((r) => r.map((x) => x / k));
    if (k % 4 === 1) S = fadd(S, term);
    else if (k % 4 === 2) C = fadd(C, term, -1);
    else if (k % 4 === 3) S = fadd(S, term, -1);
    else C = fadd(C, term);
  }
  for (let i = 0; i < s; i += 1) {
    const C2 = fadd(fmul(C, C).map((r) => r.map((x) => 2 * x)), feye(n), -1);
    S = fmul(S, C).map((r) => r.map((x) => 2 * x));
    C = C2;
  }
  return [C, S];
}

// ---- certified square root -------------------------------------------------------------------

/** Solve the n²×n² float system for the inverse of X ↦ XS + SX (row-major), or null. */
function sylvesterInverse(S) {
  const n = S.length;
  const N = n * n;
  const J = Array.from({ length: N }, () => new Array(N).fill(0));
  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j < n; j += 1) {
      const r = i * n + j;
      // (XS + SX)_ij = Σ_l X_il S_lj + Σ_k S_ik X_kj
      for (let l = 0; l < n; l += 1) J[r][i * n + l] += S[l][j];
      for (let k = 0; k < n; k += 1) J[r][k * n + j] += S[i][k];
    }
  }
  return finv(J);
}

/** Interval LDLᵀ without pivoting: true when every symmetric matrix in the interval matrix is positive definite. */
function provenPositiveDefinite(K, M) {
  const n = M.length;
  const W = M.map((row) => row.map((x) => [x[0], x[1]]));
  for (let c = 0; c < n; c += 1) {
    const piv = W[c][c];
    if (piv[0] <= 0n) return false;
    for (let r = c + 1; r < n; r += 1) {
      // factor = W[r][c] / piv
      const num = W[r][c];
      const cands = [floorDiv(num[0] * K.one, piv[0]), floorDiv(num[0] * K.one, piv[1]), floorDiv(num[1] * K.one, piv[0]), floorDiv(num[1] * K.one, piv[1])];
      let lo = cands[0];
      let hi = cands[0];
      for (const v of cands) { if (v < lo) lo = v; if (v > hi) hi = v; }
      const f = [lo, hi + 1n];
      for (let j = c + 1; j < n; j += 1) W[r][j] = sub(W[r][j], K.mul(f, W[c][j]));
    }
  }
  return true;
}

/**
 * Principal square root of an interval matrix A: an interval matrix that
 * contains the principal root of every matrix in A. Throws RangeError when
 * the proof does not go through.
 */
function sqrtIv(K, A) {
  const n = A.length;
  const N = n * n;
  const mid = A.map((row) => row.map(([lo, hi]) => K.toNumber((lo + hi) >> 1n)));
  const S = sqrtmFloat(mid);
  const Cf = S && sylvesterInverse(S);
  if (!S || !Cf) throw new RangeError("no principal square root was found: the matrix is singular or too close to one with a negative real eigenvalue");
  const C = Cf.map((row) => row.map((x) => { const v = K.fromFloat(x); return [v, v]; }));
  const point = (X) => X.map((row) => row.map((v) => [v, v]));
  let m = S.map((row) => row.map(K.fromFloat));
  const residual = (mp) => msub(K.matmul(point(mp), point(mp)), A); // m² − A as intervals
  const apply = (Mat, vec) => Mat.map((row) => K.dot(row, vec));
  // Newton with the fixed approximate inverse: m ← m − C·(m² − mid A)
  let step = null;
  let last = null;
  for (let it = 0; it < 60; it += 1) {
    step = apply(C, residual(m).flat());
    const size = step.reduce((s, x) => (mag(x) > s ? mag(x) : s), 0n);
    if (last !== null && size * 2n >= last) break;
    last = size;
    m = m.map((row, i) => row.map((v, j) => v - ((step[i * n + j][0] + step[i * n + j][1]) >> 1n)));
  }
  step = apply(C, residual(m).flat());
  let r = step.reduce((s, x) => (mag(x) > s ? mag(x) : s), 0n) * 2n + 4n;
  for (let attempt = 0; attempt < 12; attempt += 1, r *= 4n) {
    const X = m.map((row) => row.map((v) => [v - r, v + r]));
    // I − C·J(X), with J(X)[ij, il] += X_lj and J(X)[ij, kj] += X_ik
    const G = Array.from({ length: N }, (_, a) => Array.from({ length: N }, (_, b) => (a === b ? [K.one, K.one] : ZERO)));
    for (let a = 0; a < N; a += 1) {
      const Ca = C[a];
      const acc = Array.from({ length: N }, () => [0n, 0n]);
      for (let i = 0; i < n; i += 1) {
        for (let j = 0; j < n; j += 1) {
          const c = Ca[i * n + j][0];
          if (c === 0n) continue;
          for (let l = 0; l < n; l += 1) {
            const x = X[l][j];
            const t = acc[i * n + l];
            if (c > 0n) { t[0] += c * x[0]; t[1] += c * x[1]; } else { t[0] += c * x[1]; t[1] += c * x[0]; }
          }
          for (let k = 0; k < n; k += 1) {
            const x = X[i][k];
            const t = acc[k * n + j];
            if (c > 0n) { t[0] += c * x[0]; t[1] += c * x[1]; } else { t[0] += c * x[1]; t[1] += c * x[0]; }
          }
        }
      }
      for (let b = 0; b < N; b += 1) G[a][b] = sub(G[a][b], [K.down(acc[b][0]), K.up(acc[b][1])]);
    }
    const ball = Array.from({ length: N }, () => [-r, r]);
    const spread = apply(G, ball);
    let inside = true;
    const Kbox = [];
    for (let a = 0; a < N; a += 1) {
      const centre = m[Math.floor(a / n)][a % n];
      const lo = centre - step[a][1] + spread[a][0];
      const hi = centre - step[a][0] + spread[a][1];
      if (!(lo > centre - r && hi < centre + r)) inside = false;
      Kbox.push([lo, hi]);
    }
    if (!inside) continue;
    const root = Array.from({ length: n }, (_, i) => Kbox.slice(i * n, i * n + n));
    // principal: find P ≻ 0 with X P + P Xᵀ ≻ 0 for every X in the box
    const Pf = lyapunov(S);
    if (!Pf) break;
    const Pm = point(Pf.map((row) => row.map(K.fromFloat)));
    const XP = K.matmul(root, Pm);
    const Msym = XP.map((row, i) => row.map((x, j) => add(x, XP[j][i])));
    if (!provenPositiveDefinite(K, Pm) || !provenPositiveDefinite(K, Msym)) break;
    return root;
  }
  throw new RangeError("the principal square root could not be proven at this precision");
}

/** Float solution P = Pᵀ of S P + P Sᵀ = I, or null. */
function lyapunov(S) {
  const n = S.length;
  const N = n * n;
  const L = Array.from({ length: N }, () => new Array(N).fill(0));
  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j < n; j += 1) {
      const r = i * n + j;
      for (let k = 0; k < n; k += 1) {
        L[r][k * n + j] += S[i][k]; // (S P)_ij
        L[r][i * n + k] += S[j][k]; // (P Sᵀ)_ij
      }
    }
  }
  const Li = finv(L);
  if (!Li) return null;
  const rhs = feye(n).flat();
  const p = Li.map((row) => row.reduce((s, x, k) => s + x * rhs[k], 0));
  const Pm = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (p[i * n + j] + p[j * n + i]) / 2));
  return Pm.flat().every(Number.isFinite) ? Pm : null;
}

/** Principal logarithm of an interval matrix. Generator. */
function* logIv(K, A) {
  const n = A.length;
  const I = K.identity(n);
  let R = A;
  let k = 0;
  const quarter = K.one >> 2n;
  while (K.norm(msub(R, I)) > quarter) {
    if (k >= 64) throw new RangeError("the logarithm did not come within reach of its series");
    R = sqrtIv(K, R);
    k += 1;
    yield { step: "root", done: k };
  }
  const X = msub(R, I);
  const x = K.norm(X); // ≤ 1/4
  // tail of Σ (−1)^(j+1) X^j / j after J terms ≤ x^(J+1) / ((J+1)(1 − x)) ≤ 2·x^(J+1)
  let J = 1;
  let pow = ceilDiv(x * x, K.one) + 1n;
  while (2n * pow > 1n && J < 4000) {
    J += 1;
    pow = ceilDiv(pow * x, K.one) + 1n;
  }
  let term = I;
  let sum = Array.from({ length: n }, () => Array.from({ length: n }, () => ZERO));
  for (let j = 1; j <= J; j += 1) {
    term = K.matmul(term, X);
    const part = divInt(term, BigInt(j));
    sum = j % 2 ? madd(sum, part) : msub(sum, part);
    if (j % 6 === 0) yield { step: "series", done: j, of: J };
  }
  return shiftLeft(widen(sum, 2n * pow + 2n), k);
}

// ---- public generators -------------------------------------------------------------------

const toPairs = (K, X) => X.map((row) => row.map(([lo, hi]) => [[lo, K.one], [hi, K.one]]));

/**
 * Run `compute(kernel, A as intervals)` at rising precision until every
 * entry is narrower than 2^-bits. Returns pairs intervals [[n, d], [n, d]].
 */
function* sharpen(pairs, bits, extra, compute) {
  let P = bits + 24 + extra;
  let lastError = null;
  for (let attempt = 0; attempt < 6; attempt += 1) {
    const K = kernel(P);
    const A = pairs.map((row) => row.map(K.fromPair));
    let E;
    try {
      E = yield* compute(K, A);
    } catch (e) {
      if (!(e instanceof RangeError) || attempt >= 3) throw e;
      lastError = e;
      P += 64;
      continue;
    }
    const w = widest(E);
    if (w <= (1n << BigInt(P - bits))) return toPairs(K, E);
    P += w.toString(2).length - (P - bits) + 16;
  }
  throw lastError ?? new Error("matrix function: could not reach the requested width");
}

const growth = (pairs) => {
  // bits lost to the size of the answer: about ‖A‖/ln 2
  const norm = Math.max(0, ...pairs.map((row) => row.reduce((s, [n, d]) => s + Math.abs(Number(n) / Number(d)), 0)));
  if (!Number.isFinite(norm)) throw new RangeError("matrix function: entries too large");
  return Math.ceil(norm * 1.45) + 2 * Math.ceil(Math.log2(norm + 2)) * Math.ceil(Math.log2(pairs.length + 1));
};

/** Generators returning interval matrices as pairs. */
export function* cosmPairs(pairs, bits = 53) {
  return yield* sharpen(pairs, bits, growth(pairs), function* (K, A) { return (yield* cosSinIv(K, A))[0]; });
}
export function* sinmPairs(pairs, bits = 53) {
  return yield* sharpen(pairs, bits, growth(pairs), function* (K, A) { return (yield* cosSinIv(K, A))[1]; });
}
export function* coshmPairs(pairs, bits = 53) {
  return yield* sharpen(pairs, bits, growth(pairs), function* (K, A) {
    return shiftRight(madd(yield* expIv(K, A), yield* expIv(K, mneg(A))), 1);
  });
}
export function* sinhmPairs(pairs, bits = 53) {
  return yield* sharpen(pairs, bits, growth(pairs), function* (K, A) {
    return shiftRight(msub(yield* expIv(K, A), yield* expIv(K, mneg(A))), 1);
  });
}
export function* sqrtmPairs(pairs, bits = 53) {
  return yield* sharpen(pairs, bits, 16, function* (K, A) { return sqrtIv(K, A); });
}
export function* logmPairs(pairs, bits = 53) {
  return yield* sharpen(pairs, bits, 40, function* (K, A) { return yield* logIv(K, A); });
}
/** A^t for a rational t = [n, d], as exp(t·log A). */
export function* powmPairs(pairs, t, bits = 53) {
  const tf = Math.abs(Number(t[0]) / Number(t[1]));
  return yield* sharpen(pairs, bits, 40 + Math.ceil(tf * 8), function* (K, A) {
    const L = yield* logIv(K, A);
    const tv = K.fromPair(t);
    return yield* expIv(K, L.map((row) => row.map((x) => K.mul(tv, x))));
  });
}
