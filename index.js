// @ratmath/linalg: exact rational linear algebra with fast previews.
export { Matrix } from "./src/matrix.js";
export { Refinement, live } from "./src/refine.js";
export { det, rank, rref, solve, inverse, nullspace, multiply } from "./src/ops.js";
export {
  charpoly, eigenvalues, eigenvectors, singularValues, svd,
  lu, ldl, qr, gramSchmidt, power, trace,
  pseudoinverse, leastSquares, spaces, hermite, smith,
  certifyRoot, solveIntervals, similarity, jordan, expm,
} from "./src/ops2.js";
export { AlgebraicReal } from "./src/poly.js";
export { ComplexAlgebraic } from "./src/complex.js";
export { AlgebraicVector } from "./src/structure.js";
export * as poly from "./src/poly.js";
export * as interval from "./src/interval.js";
export * as teach from "./src/teach/index.js";
export * as exact from "./src/exact.js";
export * as float from "./src/float.js";
export * as field from "./src/field.js";
export * as algebraic from "./src/algebraic.js";
