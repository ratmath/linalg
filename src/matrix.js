// An immutable matrix of exact rationals. Entries are stored as two flat
// arrays of BigInt (numerators, denominators); `get` hands out Core
// Rationals. Matrices also carry a Float64Array view, made on first use,
// that the quick (preview) algorithms read.

import { toPair, toRational, pairToNumber, pairToString, addQ, subQ, mulQ, reduce } from "./q.js";

export class Matrix {
  #num;
  #den;
  #float;

  /** @internal use Matrix.from */
  constructor(rows, cols, num, den) {
    this.rows = rows;
    this.cols = cols;
    this.#num = num;
    this.#den = den;
    Object.freeze(this);
  }

  /**
   * Build from nested rows: Matrix.from([[1, "1/2"], [R`3/4`, 2n]]).
   * A flat array is read as a column vector.
   */
  static from(input) {
    if (input instanceof Matrix) return input;
    if (!Array.isArray(input)) throw new TypeError("Matrix.from expects an array of rows");
    const rows = input.map((row) => (Array.isArray(row) ? row : [row]));
    const r = rows.length;
    const c = r ? rows[0].length : 0;
    const num = new Array(r * c);
    const den = new Array(r * c);
    rows.forEach((row, i) => {
      if (row.length !== c) throw new RangeError(`Row ${i} has ${row.length} entries, expected ${c}`);
      row.forEach((x, j) => {
        const [n, d] = toPair(x);
        num[i * c + j] = n;
        den[i * c + j] = d;
      });
    });
    return new Matrix(r, c, num, den);
  }

  /** @internal build from [n, d] pairs, already reduced */
  static fromPairs(pairs) {
    const r = pairs.length;
    const c = r ? pairs[0].length : 0;
    const num = new Array(r * c);
    const den = new Array(r * c);
    for (let i = 0; i < r; i += 1) {
      for (let j = 0; j < c; j += 1) {
        num[i * c + j] = pairs[i][j][0];
        den[i * c + j] = pairs[i][j][1];
      }
    }
    return new Matrix(r, c, num, den);
  }

  static identity(n) {
    return Matrix.fromPairs(Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => [i === j ? 1n : 0n, 1n])));
  }

  static zeros(rows, cols = rows) {
    return Matrix.fromPairs(Array.from({ length: rows }, () => Array.from({ length: cols }, () => [0n, 1n])));
  }

  get isSquare() {
    return this.rows === this.cols;
  }

  /** The entry as a Core Rational. */
  get(i, j) {
    return toRational(this.pair(i, j));
  }

  /** @internal the entry as an [n, d] pair */
  pair(i, j) {
    const k = i * this.cols + j;
    return [this.#num[k], this.#den[k]];
  }

  /** Rows of [n, d] pairs (a fresh copy). */
  pairs() {
    return Array.from({ length: this.rows }, (_, i) => Array.from({ length: this.cols }, (_, j) => this.pair(i, j)));
  }

  /** Rows of Core Rationals. */
  toArray() {
    return this.pairs().map((row) => row.map(toRational));
  }

  /** Rows of strings like "3/4". */
  toStrings() {
    return this.pairs().map((row) => row.map(pairToString));
  }

  /** Rows of nearest doubles. */
  toNumbers() {
    const f = this.float();
    return Array.from({ length: this.rows }, (_, i) => Array.from(f.subarray(i * this.cols, (i + 1) * this.cols)));
  }

  /** @internal row-major Float64Array, cached */
  float() {
    if (!this.#float) {
      const f = new Float64Array(this.rows * this.cols);
      for (let k = 0; k < f.length; k += 1) f[k] = pairToNumber(this.#num[k], this.#den[k]);
      this.#float = f;
    }
    return this.#float;
  }

  map(fn) {
    return Matrix.fromPairs(this.pairs().map((row, i) => row.map((p, j) => fn(p, i, j))));
  }

  transpose() {
    return Matrix.fromPairs(Array.from({ length: this.cols }, (_, j) => Array.from({ length: this.rows }, (_, i) => this.pair(i, j))));
  }

  add(other) {
    other = Matrix.from(other);
    sameShape(this, other, "add");
    return this.map((p, i, j) => addQ(p, other.pair(i, j)));
  }

  subtract(other) {
    other = Matrix.from(other);
    sameShape(this, other, "subtract");
    return this.map((p, i, j) => subQ(p, other.pair(i, j)));
  }

  scale(x) {
    const s = toPair(x);
    return this.map((p) => mulQ(p, s));
  }

  multiply(other) {
    other = Matrix.from(other);
    if (this.cols !== other.rows) throw new RangeError(`Cannot multiply ${this.rows}×${this.cols} by ${other.rows}×${other.cols}`);
    const out = [];
    for (let i = 0; i < this.rows; i += 1) {
      const row = [];
      for (let j = 0; j < other.cols; j += 1) {
        // sum over a common denominator, reduce once
        let n = 0n;
        let d = 1n;
        for (let k = 0; k < this.cols; k += 1) {
          const [a, b] = this.pair(i, k);
          if (a === 0n) continue;
          const [c, e] = other.pair(k, j);
          if (c === 0n) continue;
          const den = b * e;
          if (den === d) n += a * c;
          else {
            n = n * den + a * c * d;
            d *= den;
          }
        }
        row.push(reduce(n, d));
      }
      out.push(row);
    }
    return Matrix.fromPairs(out);
  }

  equals(other) {
    other = Matrix.from(other);
    if (this.rows !== other.rows || this.cols !== other.cols) return false;
    for (let k = 0; k < this.#num.length; k += 1) {
      const i = Math.floor(k / this.cols);
      const [n, d] = other.pair(i, k % this.cols);
      if (this.#num[k] !== n || this.#den[k] !== d) return false;
    }
    return true;
  }

  column(j) {
    return Matrix.fromPairs(Array.from({ length: this.rows }, (_, i) => [this.pair(i, j)]));
  }

  /** Columns side by side: [this | other]. */
  augment(other) {
    other = Matrix.from(other);
    if (other.rows !== this.rows) throw new RangeError("augment: row counts differ");
    const a = this.pairs();
    const b = other.pairs();
    return Matrix.fromPairs(a.map((row, i) => [...row, ...b[i]]));
  }

  toString() {
    return this.toStrings().map((row) => `[${row.join(", ")}]`).join("\n");
  }
}

function sameShape(a, b, op) {
  if (a.rows !== b.rows || a.cols !== b.cols) throw new RangeError(`Cannot ${op} ${a.rows}×${a.cols} and ${b.rows}×${b.cols}`);
}
