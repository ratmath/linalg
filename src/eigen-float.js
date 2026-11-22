// Float64 eigenvalues (Hessenberg + Francis double-shift QR), eigenvectors
// by inverse iteration, and the SVD by one-sided Jacobi. These are previews:
// quick, usually accurate, never certified. The exact side lives in spectral.js.

/** 1-indexed copy of a row-major n×n Float64Array, as the classic routines expect. */
function oneIndexed(flat, n) {
  const a = Array.from({ length: n + 1 }, () => new Float64Array(n + 1));
  for (let i = 0; i < n; i += 1) for (let j = 0; j < n; j += 1) a[i + 1][j + 1] = flat[i * n + j];
  return a;
}

/** Reduce to upper Hessenberg form by stabilised elementary similarities. */
function elmhes(a, n) {
  for (let m = 2; m < n; m += 1) {
    let x = 0;
    let i = m;
    for (let j = m; j <= n; j += 1) if (Math.abs(a[j][m - 1]) > Math.abs(x)) { x = a[j][m - 1]; i = j; }
    if (i !== m) {
      for (let j = m - 1; j <= n; j += 1) [a[i][j], a[m][j]] = [a[m][j], a[i][j]];
      for (let j = 1; j <= n; j += 1) [a[j][i], a[j][m]] = [a[j][m], a[j][i]];
    }
    if (x) {
      for (i = m + 1; i <= n; i += 1) {
        let y = a[i][m - 1];
        if (y !== 0) {
          y /= x;
          a[i][m - 1] = y;
          for (let j = m; j <= n; j += 1) a[i][j] -= y * a[m][j];
          for (let j = 1; j <= n; j += 1) a[j][m] += y * a[j][i];
        }
      }
    }
  }
  for (let i = 3; i <= n; i += 1) for (let j = 1; j < i - 1; j += 1) a[i][j] = 0;
}

const SIGN = (a, b) => (b >= 0 ? Math.abs(a) : -Math.abs(a));

/** Eigenvalues of an upper Hessenberg matrix (destroyed). Returns { re, im } or null if it fails to converge. */
function hqr(a, n) {
  const wr = new Float64Array(n + 1);
  const wi = new Float64Array(n + 1);
  let anorm = 0;
  for (let i = 1; i <= n; i += 1) for (let j = Math.max(i - 1, 1); j <= n; j += 1) anorm += Math.abs(a[i][j]);
  let nn = n;
  let t = 0;
  let p = 0, q = 0, r = 0, s = 0, w = 0, x = 0, y = 0, z = 0;
  while (nn >= 1) {
    let its = 0;
    let l;
    do {
      for (l = nn; l >= 2; l -= 1) {
        s = Math.abs(a[l - 1][l - 1]) + Math.abs(a[l][l]);
        if (s === 0) s = anorm;
        if (Math.abs(a[l][l - 1]) + s === s) {
          a[l][l - 1] = 0;
          break;
        }
      }
      x = a[nn][nn];
      if (l === nn) {
        wr[nn] = x + t;
        wi[nn] = 0;
        nn -= 1;
      } else {
        y = a[nn - 1][nn - 1];
        w = a[nn][nn - 1] * a[nn - 1][nn];
        if (l === nn - 1) {
          p = 0.5 * (y - x);
          q = p * p + w;
          z = Math.sqrt(Math.abs(q));
          x += t;
          if (q >= 0) {
            z = p + SIGN(z, p);
            wr[nn - 1] = wr[nn] = x + z;
            if (z) wr[nn] = x - w / z;
            wi[nn - 1] = wi[nn] = 0;
          } else {
            wr[nn - 1] = wr[nn] = x + p;
            wi[nn] = z;
            wi[nn - 1] = -z;
          }
          nn -= 2;
        } else {
          if (its === 60) return null;
          if (its === 10 || its === 20 || its === 40) {
            t += x;
            for (let i = 1; i <= nn; i += 1) a[i][i] -= x;
            s = Math.abs(a[nn][nn - 1]) + Math.abs(a[nn - 1][nn - 2]);
            y = x = 0.75 * s;
            w = -0.4375 * s * s;
          }
          its += 1;
          let m;
          for (m = nn - 2; m >= l; m -= 1) {
            z = a[m][m];
            r = x - z;
            s = y - z;
            p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
            q = a[m + 1][m + 1] - z - r - s;
            r = a[m + 2][m + 1];
            s = Math.abs(p) + Math.abs(q) + Math.abs(r);
            p /= s;
            q /= s;
            r /= s;
            if (m === l) break;
            const u = Math.abs(a[m][m - 1]) * (Math.abs(q) + Math.abs(r));
            const v = Math.abs(p) * (Math.abs(a[m - 1][m - 1]) + Math.abs(z) + Math.abs(a[m + 1][m + 1]));
            if (u + v === v) break;
          }
          for (let i = m + 2; i <= nn; i += 1) {
            a[i][i - 2] = 0;
            if (i !== m + 2) a[i][i - 3] = 0;
          }
          for (let k = m; k <= nn - 1; k += 1) {
            if (k !== m) {
              p = a[k][k - 1];
              q = a[k + 1][k - 1];
              r = 0;
              if (k !== nn - 1) r = a[k + 2][k - 1];
              if ((x = Math.abs(p) + Math.abs(q) + Math.abs(r)) !== 0) {
                p /= x;
                q /= x;
                r /= x;
              }
            }
            if ((s = SIGN(Math.sqrt(p * p + q * q + r * r), p)) !== 0) {
              if (k === m) {
                if (l !== m) a[k][k - 1] = -a[k][k - 1];
              } else a[k][k - 1] = -s * x;
              p += s;
              x = p / s;
              y = q / s;
              z = r / s;
              q /= p;
              r /= p;
              for (let j = k; j <= nn; j += 1) {
                p = a[k][j] + q * a[k + 1][j];
                if (k !== nn - 1) {
                  p += r * a[k + 2][j];
                  a[k + 2][j] -= p * z;
                }
                a[k + 1][j] -= p * y;
                a[k][j] -= p * x;
              }
              const mmin = nn < k + 3 ? nn : k + 3;
              for (let i = l; i <= mmin; i += 1) {
                p = x * a[i][k] + y * a[i][k + 1];
                if (k !== nn - 1) {
                  p += z * a[i][k + 2];
                  a[i][k + 2] -= p * r;
                }
                a[i][k + 1] -= p * q;
                a[i][k] -= p;
              }
            }
          }
        }
      }
    } while (l < nn - 1);
  }
  return { re: Array.from(wr.subarray(1)), im: Array.from(wi.subarray(1)) };
}

/**
 * Eigenvalues as { re, im } pairs, sorted: real ones ascending first, then
 * complex pairs by real part. null if QR fails to converge (rare).
 */
export function eigenvaluesFloat(flat, n) {
  if (n === 0) return [];
  const a = oneIndexed(flat, n);
  // balance-free; scale to unit size to keep the tests in hqr meaningful
  elmhes(a, n);
  const res = hqr(a, n);
  if (!res) return null;
  const out = res.re.map((re, i) => ({ re, im: res.im[i] }));
  return out.sort((u, v) => (u.im !== 0) - (v.im !== 0) || u.re - v.re || u.im - v.im);
}

/** Solve M x = b by Gaussian elimination with partial pivoting (small helper). */
function gaussSolve(M, b) {
  const n = M.length;
  const a = M.map((row, i) => [...row, b[i]]);
  for (let k = 0; k < n; k += 1) {
    let best = k;
    for (let i = k + 1; i < n; i += 1) if (Math.abs(a[i][k]) > Math.abs(a[best][k])) best = i;
    [a[k], a[best]] = [a[best], a[k]];
    if (a[k][k] === 0) a[k][k] = 1e-300;
    for (let i = k + 1; i < n; i += 1) {
      const f = a[i][k] / a[k][k];
      if (f) for (let j = k; j <= n; j += 1) a[i][j] -= f * a[k][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i -= 1) {
    let s = a[i][n];
    for (let j = i + 1; j < n; j += 1) s -= a[i][j] * x[j];
    x[i] = s / a[i][i];
  }
  return x;
}

/** Normalise so the entry of largest size is 1. */
export function normaliseMax(v) {
  let k = 0;
  for (let i = 1; i < v.length; i += 1) if (Math.abs(v[i]) > Math.abs(v[k])) k = i;
  const s = v[k] || 1;
  return v.map((x) => x / s);
}

/** Eigenvector for a real eigenvalue by inverse iteration. */
export function eigenvectorFloat(flat, n, lambda) {
  const shift = lambda + (Math.abs(lambda) + 1) * 1e-10;
  const M = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => flat[i * n + j] - (i === j ? shift : 0)));
  let v = Array.from({ length: n }, (_, i) => 1 + i / (n + 1));
  for (let it = 0; it < 3; it += 1) v = normaliseMax(gaussSolve(M, v));
  return v;
}

/**
 * Thin SVD by one-sided Jacobi: A (m×n) = U·diag(s)·Vᵀ, s descending,
 * U m×k, V n×k with k = min(m, n). Returns { U, s, V } as plain arrays.
 */
export function svdFloat(flat, m, n) {
  if (m < n) {
    const t = new Float64Array(m * n);
    for (let i = 0; i < m; i += 1) for (let j = 0; j < n; j += 1) t[j * m + i] = flat[i * n + j];
    const { U, s, V } = svdFloat(t, n, m);
    return { U: V, s, V: U };
  }
  const U = Array.from({ length: m }, (_, i) => Array.from(flat.subarray(i * n, (i + 1) * n)));
  const V = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  for (let sweep = 0; sweep < 60; sweep += 1) {
    let off = 0;
    for (let j = 0; j < n - 1; j += 1) {
      for (let k = j + 1; k < n; k += 1) {
        let alpha = 0, beta = 0, gamma = 0;
        for (let i = 0; i < m; i += 1) {
          alpha += U[i][j] * U[i][j];
          beta += U[i][k] * U[i][k];
          gamma += U[i][j] * U[i][k];
        }
        if (gamma === 0 || Math.abs(gamma) <= 1e-15 * Math.sqrt(alpha * beta)) continue;
        off = Math.max(off, Math.abs(gamma) / Math.sqrt(alpha * beta));
        const zeta = (beta - alpha) / (2 * gamma);
        const t = Math.sign(zeta || 1) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
        const c = 1 / Math.sqrt(1 + t * t);
        const s = c * t;
        for (let i = 0; i < m; i += 1) {
          const x = U[i][j];
          U[i][j] = c * x - s * U[i][k];
          U[i][k] = s * x + c * U[i][k];
        }
        for (let i = 0; i < n; i += 1) {
          const x = V[i][j];
          V[i][j] = c * x - s * V[i][k];
          V[i][k] = s * x + c * V[i][k];
        }
      }
    }
    if (off < 1e-15) break;
  }
  const s = Array.from({ length: n }, (_, j) => Math.hypot(...U.map((row) => row[j])));
  const order = s.map((_, j) => j).sort((a, b) => s[b] - s[a]);
  return {
    s: order.map((j) => s[j]),
    U: U.map((row) => order.map((j) => (s[j] ? row[j] / s[j] : 0))),
    V: V.map((row) => order.map((j) => row[j])),
  };
}

/** Householder QR, thin: A (m×n, m ≥ n) = Q R, Q m×n orthonormal columns. */
export function qrFloat(flat, m, n) {
  const R = Array.from({ length: m }, (_, i) => Array.from(flat.subarray(i * n, (i + 1) * n)));
  const Q = Array.from({ length: m }, (_, i) => Array.from({ length: m }, (_, j) => (i === j ? 1 : 0)));
  for (let k = 0; k < Math.min(m - 1, n); k += 1) {
    let norm = 0;
    for (let i = k; i < m; i += 1) norm += R[i][k] ** 2;
    norm = Math.sqrt(norm);
    if (norm === 0) continue;
    const alpha = R[k][k] > 0 ? -norm : norm;
    const v = new Array(m).fill(0);
    for (let i = k; i < m; i += 1) v[i] = R[i][k];
    v[k] -= alpha;
    const vv = v.reduce((a, x) => a + x * x, 0);
    if (vv === 0) continue;
    for (let j = 0; j < n; j += 1) {
      let d = 0;
      for (let i = k; i < m; i += 1) d += v[i] * R[i][j];
      d = (2 * d) / vv;
      for (let i = k; i < m; i += 1) R[i][j] -= d * v[i];
    }
    for (let i = 0; i < m; i += 1) {
      let d = 0;
      for (let j = k; j < m; j += 1) d += Q[i][j] * v[j];
      d = (2 * d) / vv;
      for (let j = k; j < m; j += 1) Q[i][j] -= d * v[j];
    }
  }
  // make R's diagonal non-negative
  for (let k = 0; k < Math.min(m, n); k += 1) {
    if (R[k][k] < 0) {
      for (let j = 0; j < n; j += 1) R[k][j] = -R[k][j];
      for (let i = 0; i < m; i += 1) Q[i][k] = -Q[i][k];
    }
  }
  const k = Math.min(m, n);
  return { Q: Q.map((row) => row.slice(0, k)), R: R.slice(0, k) };
}
