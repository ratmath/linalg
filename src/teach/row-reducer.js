// An interactive row reducer: the learner chooses the row operations, the
// reducer applies them exactly, keeps history for undo, says whether the
// matrix is in echelon or reduced form yet, and can suggest the next step.

import { Matrix } from "../matrix.js";
import { toPair, subQ, mulQ, divQ, addQ } from "../q.js";

const fmt = ([n, d]) => (d === 1n ? `${n}` : `${n}/${d}`);
const isZ = ([n]) => n === 0n;
const sub = (k) => String(k).split("").map((c) => "₀₁₂₃₄₅₆₇₈₉"[c]).join("");
const R = (i) => `R${sub(i + 1)}`;

export class RowReducer {
  constructor(input, { bar = null } = {}) {
    this.M = Matrix.from(input).pairs();
    this.bar = bar;
    this.history = [];
  }

  get matrix() {
    return this.M.map((r) => r.map(fmt));
  }

  get columns() {
    return this.bar ?? this.M[0].length;
  }

  #push(label, next) {
    this.history.push({ label, M: this.M });
    this.M = next;
    return label;
  }

  /** Rᵢ ↔ Rⱼ (0-based rows). */
  swap(i, j) {
    const next = this.M.map((r) => [...r]);
    [next[i], next[j]] = [next[j], next[i]];
    return this.#push(`${R(i)} ↔ ${R(j)}`, next);
  }

  /** Rᵢ ← c·Rᵢ, c ≠ 0. */
  scale(i, c) {
    const k = toPair(c);
    if (isZ(k)) throw new RangeError("Scaling a row by 0 is not allowed: it would lose information");
    const next = this.M.map((r, t) => (t === i ? r.map((v) => mulQ(v, k)) : [...r]));
    return this.#push(`${R(i)} ← ${fmt(k)}·${R(i)}`, next);
  }

  /** Rᵢ ← Rᵢ + c·Rⱼ, i ≠ j. */
  add(i, j, c) {
    if (i === j) throw new RangeError("Add a multiple of a different row");
    const k = toPair(c);
    const next = this.M.map((r, t) => (t === i ? r.map((v, col) => addQ(v, mulQ(k, this.M[j][col]))) : [...r]));
    return this.#push(`${R(i)} ← ${R(i)} ${k[0] < 0n ? "−" : "+"} ${fmt([k[0] < 0n ? -k[0] : k[0], k[1]])}·${R(j)}`, next);
  }

  undo() {
    const last = this.history.pop();
    if (last) this.M = last.M;
    return last?.label ?? null;
  }

  #leads() {
    return this.M.map((r) => {
      const j = r.slice(0, this.columns).findIndex((v) => !isZ(v));
      return j;
    });
  }

  /** Row echelon form: zero rows at the bottom, leading entries move right. */
  isEchelon() {
    const leads = this.#leads();
    let last = -1;
    let seenZero = false;
    for (const j of leads) {
      if (j < 0) {
        seenZero = true;
        continue;
      }
      if (seenZero || j <= last) return false;
      last = j;
    }
    return true;
  }

  /** Reduced: echelon, leading entries 1, and alone in their columns. */
  isReduced() {
    if (!this.isEchelon()) return false;
    return this.#leads().every((j, i) => j < 0 || (this.M[i][j][0] === 1n && this.M[i][j][1] === 1n && this.M.every((r, t) => t === i || isZ(r[j]))));
  }

  /** The next operation Gauss–Jordan would do from here, as { label, apply() }. */
  hint() {
    const m = this.M.length;
    let r = 0;
    for (let col = 0; col < this.columns && r < m; col += 1) {
      let p = r;
      while (p < m && isZ(this.M[p][col])) p += 1;
      if (p === m) continue;
      const done = this.M[r][col][0] === 1n && this.M[r][col][1] === 1n && this.M.every((row, i) => i === r || isZ(row[col]));
      if (done && p === r) {
        r += 1;
        continue;
      }
      if (p !== r && isZ(this.M[r][col])) return { label: `${R(r)} ↔ ${R(p)}`, why: `Get a nonzero pivot into row ${r + 1}, column ${col + 1}.`, apply: () => this.swap(r, p) };
      const piv = this.M[r][col];
      if (!(piv[0] === 1n && piv[1] === 1n)) return { label: `${R(r)} ← ${fmt(divQ([1n, 1n], piv))}·${R(r)}`, why: "Make the pivot 1.", apply: () => this.scale(r, divQ([1n, 1n], piv)) };
      const i = this.M.findIndex((row, t) => t !== r && !isZ(row[col]));
      const f = this.M[i][col];
      return { label: `${R(i)} ← ${R(i)} − ${fmt(f)}·${R(r)}`, why: `Clear column ${col + 1} in row ${i + 1}.`, apply: () => this.add(i, r, subQ([0n, 1n], f)) };
    }
    return null;
  }
}
