// Docs bundle entry: the library plus the Core pieces the pages use directly.
export * from "../index.js";
export { Rational, Integer, RationalInterval, parseNumber } from "@ratmath/core";
export { serveLinalg, workerClient } from "../src/worker.js";
