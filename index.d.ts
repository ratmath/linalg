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
  /** For a rational eigenvalue: a basis of the eigenspace. Otherwise rational approximations of `basis`. */
  vectors: Matrix[];
  /** For an irrational or complex eigenvalue: an exact basis of the eigenspace over Q(λ). */
  basis?: AlgebraicVector[];
  /** The first vector of `basis`. */
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
  /** J and P with exact algebraic entries for any eigenvalues; null if a complex eigenvalue was not certified. Computed when called. */
  algebraicForm(): { J: AlgebraicMatrix; P: AlgebraicMatrix } | null;
}

/** An exact algebraic number, real or complex. */
export type Algebraic = AlgebraicReal | ComplexAlgebraic;
/** Rows of exact algebraic numbers. */
export type AlgebraicMatrix = Algebraic[][];
export type NormKind = 1 | 2 | "inf" | "fro";
/** A matrix norm, exactly: 1 and "inf" are rational, "fro" a square root, 2 the largest singular value. */
export function norm(A: MatrixInput, p?: NormKind): Refinement<number, Algebraic>;
/** The condition number ‖A‖·‖A⁻¹‖; Infinity for a singular matrix. */
export function cond(A: MatrixInput, p?: NormKind): Refinement<number | null, Algebraic | number>;
export interface Definiteness {
  kind: "positive-definite" | "positive-semidefinite" | "negative-definite" | "negative-semidefinite" | "indefinite" | "zero";
  /** Counts of positive, negative and zero eigenvalues. */
  positive: number; negative: number; zero: number;
}
export function definiteness(A: MatrixInput): Refinement<null, Definiteness>;
/** A = L·Lᵀ with exact square roots, or null when A is not positive definite. D holds the pivots of the LDLᵀ form. */
export function cholesky(A: MatrixInput): Refinement<{ L: number[][] } | null, { L: AlgebraicMatrix; D: Rational[] } | null>;
/** A = Q·diag(values)·Qᵀ with orthonormal columns in Q, exactly (symmetric A). */
export function symmetricEigen(A: MatrixInput): Refinement<null, { values: Algebraic[]; Q: AlgebraicMatrix }>;
/** A basis of the integer solutions of A·x = 0 (columns, LLL-reduced). */
export function integerKernel(A: MatrixInput): Matrix;
/** A basis of the lattice spanned by the columns of A (columns, LLL-reduced). */
export function integerImage(A: MatrixInput): Matrix;
/** The rows of an integer matrix replaced by an LLL-reduced basis of the same lattice. */
export function lllReduce(A: MatrixInput): Matrix;
/** A basis of the eigenspace of an irrational or complex eigenvalue over Q(λ); `pairs` is Matrix.pairs(). */
export function eigenspace(pairs: [bigint, bigint][][], value: Algebraic): AlgebraicVector[];
export function jordanAlgebraic(pairs: [bigint, bigint][][], eigen: unknown): { J: AlgebraicMatrix; P: AlgebraicMatrix } | null;
/** The operations record that the functions in `field` take. */
export interface FieldOps<T> { zero: T; one: T; add(a: T, b: T): T; sub(a: T, b: T): T; mul(a: T, b: T): T; div(a: T, b: T): T; neg(a: T): T; sign(a: T): 0 | 1 | null }
/** Exact algebraic numbers, real or complex, as a field. */
export const complexField: FieldOps<Algebraic>;
/** Rows of algebraic numbers from rationals and algebraic numbers. */
export function algebraicMatrix(rows: (Entry | Algebraic)[][]): AlgebraicMatrix;
export function numberFieldOps(K: unknown): FieldOps<any>;
export function similarity(A: MatrixInput): Refinement<{ eigenvalues: Complex[] | null }, Similarity>;
export const jordan: typeof similarity;
export function expm(A: MatrixInput, options?: { bits?: number }): Refinement<number[][], RationalInterval[][]>;
/** sin(A), cos(A), sinh(A), cosh(A): every entry as a certified interval narrower than 2^-bits. */
export function sinm(A: MatrixInput, options?: { bits?: number }): Refinement<number[][], RationalInterval[][]>;
export function cosm(A: MatrixInput, options?: { bits?: number }): Refinement<number[][], RationalInterval[][]>;
export function sinhm(A: MatrixInput, options?: { bits?: number }): Refinement<number[][], RationalInterval[][]>;
export function coshm(A: MatrixInput, options?: { bits?: number }): Refinement<number[][], RationalInterval[][]>;
/**
 * The principal square root and logarithm. Both need a matrix with no
 * eigenvalue on the negative real axis or at 0 and throw a RangeError
 * otherwise. The preview is `null` when the float iteration does not settle.
 */
export function sqrtm(A: MatrixInput, options?: { bits?: number }): Refinement<number[][] | null, RationalInterval[][]>;
export function logm(A: MatrixInput, options?: { bits?: number }): Refinement<number[][] | null, RationalInterval[][]>;
/** A^t for a rational t, as exp(t·log A) with the principal logarithm. */
export function powm(A: MatrixInput, t: Entry, options?: { bits?: number }): Refinement<number[][] | null, RationalInterval[][]>;

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
/** Exact arithmetic on AlgebraicReal: add, sub, mul, div, inv, neg, pow, sqrt, sign, compare, realRoots, fromRational. */
export const algebraic: Record<string, any>;
export const poly: Record<string, any>;
export const interval: Record<string, any>;
export const teach: Record<string, any>;

/** One step of an iterative solver: the double iterate and three exact distances (largest-entry norm). */
export interface MethodRow {
  n: number;
  x: number[];
  /** from the double iterate to the exact solution */
  error: number;
  /** from the exact-rational iterate to the exact solution (null once exact iterates are dropped) */
  methodError: number | null;
  /** from the double iterate to the exact-rational iterate */
  rounding: number | null;
  residual: number;
}
export interface StationaryRun {
  method: string;
  solution: string[];
  rows: MethodRow[];
  iterationMatrix: string[][];
  spectralRadius: number;
  converges: boolean;
  strictlyDiagonallyDominant: boolean;
  ratio: number | null;
}
export interface MethodOptions {
  x0?: unknown[];
  iterations?: number;
  exactIterations?: number;
  tolerance?: number;
  maxBits?: number;
}
export interface EigenRow { n: number; x: number[]; estimate: number; exactEstimate: string | null; error: number | null; methodError: number | null }
export interface EigenRun {
  method: string;
  eigenvalues: { re: number; im: number }[];
  target: number | null;
  note: string | null;
  convergenceFactor: number | null;
  rows: EigenRow[];
}
/** Classical numerical methods run in doubles beside exact arithmetic (see the "Iterative and eigenvalue methods" page). */
export const methods: {
  exactDouble(x: number): [bigint, bigint];
  jacobi(A: unknown, b: unknown[], options?: MethodOptions): StationaryRun;
  gaussSeidel(A: unknown, b: unknown[], options?: MethodOptions): StationaryRun;
  sor(A: unknown, b: unknown[], options?: MethodOptions & { omega?: number | string }): StationaryRun;
  conjugateGradient(A: unknown, b: unknown[], options?: { iterations?: number }): { method: string; solution: string[]; rows: MethodRow[]; exactSteps: number };
  powerMethod(A: unknown, options?: MethodOptions): EigenRun;
  inverseIteration(A: unknown, options?: MethodOptions & { shift?: number | string }): EigenRun;
  qrAlgorithm(A: unknown, options?: { iterations?: number }): {
    method: string;
    eigenvalues: { re: number; im: number }[];
    rows: { n: number; diagonal: number[]; subdiagonal: number }[];
    final: number[][];
  };
  pivoting(A: unknown, b: unknown[]): {
    exact: string[];
    none: { x: number[]; error: number; growth: number };
    partial: { x: number[]; error: number; growth: number };
  };
};
