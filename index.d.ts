// Type declarations for @ratmath/linalg.
import type { Rational, RationalInterval, Integer } from "@ratmath/core";

/** Anything one matrix entry can be read from. */
export type Entry = number | bigint | string | Rational | Integer | { numerator: bigint; denominator: bigint };
/** Rows of entries, a flat column, or a Matrix. */
export type MatrixInput = Matrix | Entry[][] | Entry[];
export type Trust = "good" | "rough" | "unsure";

export class Matrix {
  readonly rows: number;
  readonly cols: number;
  readonly isSquare: boolean;
  static from(input: MatrixInput): Matrix;
  static identity(n: number): Matrix;
  static zeros(rows: number, cols?: number): Matrix;
  get(i: number, j: number): Rational;
  toArray(): Rational[][];
  toStrings(): string[][];
  toNumbers(): number[][];
  transpose(): Matrix;
  add(other: MatrixInput): Matrix;
  subtract(other: MatrixInput): Matrix;
  scale(x: Entry): Matrix;
  multiply(other: MatrixInput): Matrix;
  equals(other: MatrixInput): boolean;
  column(j: number): Matrix;
  augment(other: MatrixInput): Matrix;
  toString(): string;
}

export interface FloatStage<V> { stage: "float"; exact: false; value: V; trust: Trust }
export interface ProgressStage { stage: "progress"; exact: false; step: string; done: number; of: number; [extra: string]: unknown }
export interface ExactStage<V> { stage: "exact"; exact: true; value: V; via: string }

export interface SettleOptions<Q, E> {
  signal?: AbortSignal;
  sliceMs?: number;
  onStage?: (stage: FloatStage<Q> | ExactStage<E>) => void;
  onProgress?: (stage: ProgressStage) => void;
}

/** A float preview now, the exact answer later. */
export class Refinement<Q, E> {
  constructor(quick: () => { value: Q; trust: Trust }, exact: () => Generator<Record<string, unknown>, { value: E; via: string }>);
  readonly quick: FloatStage<Q>;
  readonly settled: boolean;
  stages(): Generator<FloatStage<Q> | ProgressStage | ExactStage<E>>;
  exact(): E;
  exactStage(): ExactStage<E>;
  settle(options?: SettleOptions<Q, E>): Promise<E>;
}

export interface LiveUpdate<A extends unknown[], Q, E> {
  (...args: A): Refinement<Q, E>;
  cancel(): void;
}
export function live<A extends unknown[], Q, E>(
  compute: (...args: A) => Refinement<Q, E>,
  onStage: (stage: FloatStage<Q> | ExactStage<E>, job: Refinement<Q, E>) => void,
  options?: { idleMs?: number; sliceMs?: number; onProgress?: (p: ProgressStage) => void; onError?: (e: unknown) => void },
): LiveUpdate<A, Q, E>;

/** A real algebraic number: a root of a square-free integer polynomial in an isolating interval. */
export class AlgebraicReal {
  readonly poly: bigint[];
  lo: [bigint, bigint];
  hi: [bigint, bigint];
  static rational(pair: [bigint, bigint]): AlgebraicReal;
  refine(bits?: number): this;
  readonly rational: [bigint, bigint] | null;
  readonly isRational: boolean;
  interval(bits?: number): [Rational, Rational];
  toRational(): Rational | null;
  toNumber(): number;
  toString(digits?: number): string;
}

/** A complex algebraic number: the only root of `poly` in a certified box. */
export class ComplexAlgebraic {
  readonly poly: bigint[];
  conjugate(): ComplexAlgebraic;
  refine(bits?: number): this;
  box(bits?: number): { re: [Rational, Rational]; im: [Rational, Rational] };
  toNumbers(): { re: number; im: number };
  toString(digits?: number): string;
}

/** An eigenvector with entries in Q(λ): polynomials in λ. */
export class AlgebraicVector {
  readonly lambda: AlgebraicReal | ComplexAlgebraic;
  toStrings(variable?: string): string[];
  approximate(bits?: number): Matrix;
  toNumbers(): number[] | { re: number; im: number }[];
}

export interface SolveValue<X> { status: "unique" | "infinite" | "none"; x: X | null; basis: X[] }

export function det(A: MatrixInput): Refinement<number, Rational>;
export function rank(A: MatrixInput): Refinement<number, number>;
export function rref(A: MatrixInput, options?: { columns?: number }): Refinement<{ matrix: number[][]; pivots: number[] }, { matrix: Matrix; pivots: number[] }>;
export function solve(A: MatrixInput, b: MatrixInput, options?: { guessDenominator?: number }): Refinement<SolveValue<number[][]>, SolveValue<Matrix>>;
export function inverse(A: MatrixInput, options?: { guessDenominator?: number }): Refinement<number[][] | null, Matrix | null>;
export function nullspace(A: MatrixInput): Refinement<number[][], Matrix[]>;
export function multiply(A: MatrixInput, B: MatrixInput): Matrix;

export interface Complex { re: number; im: number }
export interface EigenValues {
  real: { value: AlgebraicReal; multiplicity: number }[];
  complex: { value: ComplexAlgebraic | null; multiplicity: number; approx?: Complex }[];
  complexCertified: boolean;
  complexCount: number;
  charpoly: Rational[];
}
export interface EigenVectorEntry {
  value: AlgebraicReal | ComplexAlgebraic;
  multiplicity: number;
  exact: boolean;
  vectors: Matrix[];
  algebraic?: AlgebraicVector | null;
  approximate: ((bits?: number) => Matrix | null) | null;
}
export function charpoly(A: MatrixInput): Refinement<number[] | null, Rational[]>;
export function eigenvalues(A: MatrixInput): Refinement<Complex[] | null, EigenValues>;
export function eigenvectors(A: MatrixInput): Refinement<{ value: number; vector: number[] }[] | null, EigenVectorEntry[]>;
export function singularValues(A: MatrixInput): Refinement<number[], { value: AlgebraicReal; multiplicity: number }[]>;
export interface FloatSVD { U: number[][]; s: number[]; V: number[][] }
export function svd(A: MatrixInput): Refinement<FloatSVD, FloatSVD & { values: { value: AlgebraicReal; multiplicity: number }[] }>;

export function lu(A: MatrixInput): Refinement<{ perm: number[]; L: number[][]; U: number[][] }, { perm: number[]; L: Matrix; U: Matrix }>;
export function ldl(A: MatrixInput): Refinement<{ L: number[][]; D: number[]; positiveDefinite: boolean }, { L: Matrix; D: Rational[]; positiveDefinite: boolean } | null>;
export function qr(A: MatrixInput): Refinement<{ Q: number[][]; R: number[][] } | null, { Q: Matrix; R: Matrix; normsSquared: Rational[] }>;
export function gramSchmidt(A: MatrixInput): Refinement<number[][] | null, Matrix[]>;
export function power(A: MatrixInput, k: number | bigint): Refinement<number[][] | null, Matrix | null>;
export function trace(A: MatrixInput): Rational;
export function pseudoinverse(A: MatrixInput): Refinement<number[][], Matrix>;
export function leastSquares(A: MatrixInput, b: MatrixInput): Refinement<{ x: number[][]; residual: number[][] }, { x: Matrix; residual: Matrix }>;
export function spaces(A: MatrixInput): { rank: number; nullity: number; columnSpace: Matrix[]; rowSpace: Matrix[] };
export function hermite(A: MatrixInput): { H: Matrix; U: Matrix; pivots: number[] };
export function smith(A: MatrixInput): { S: Matrix; U: Matrix; V: Matrix; diagonal: bigint[] };

export interface Similarity {
  invariantFactors: Rational[][];
  minimalPolynomial: Rational[];
  rationalCanonical: Matrix;
  jordan: { value: AlgebraicReal | ComplexAlgebraic; multiplicity: number; blocks: number[] }[];
  J: Matrix | null;
  P: Matrix | null;
}
export function similarity(A: MatrixInput): Refinement<{ eigenvalues: Complex[] | null }, Similarity>;
export const jordan: typeof similarity;
export function expm(A: MatrixInput, options?: { bits?: number }): Refinement<number[][], RationalInterval[][]>;

type IntervalLike = RationalInterval | [Entry, Entry] | Entry;
export function certifyRoot(options: {
  f: (x: Rational[]) => (Rational | RationalInterval)[];
  jacobian: (box: RationalInterval[]) => IntervalLike[][];
  guess: Entry[];
  bits?: number;
  steps?: number;
}): Refinement<number[], { status: "unique" | "none" | "unknown"; box: RationalInterval[] | null }>;
export function solveIntervals(A: IntervalLike[][], b: IntervalLike[], options?: { bits?: number }): Refinement<number[] | null, { status: "unique" | "unknown"; box: RationalInterval[] | null }>;

/* Lower-level modules; see the README for what each exports. */
export const exact: Record<string, any>;
export const float: Record<string, any>;
export const field: Record<string, any>;
export const poly: Record<string, any>;
export const interval: Record<string, any>;
export const teach: Record<string, any>;
