// Run the exact stage in a Web Worker so even long computations never touch
// the page's main thread. The preview still runs on the calling side (it is
// cheap and must be instant); only the exact work is shipped out.
//
// In the worker script:
//   import { serveLinalg } from "@ratmath/linalg/worker";
//   serveLinalg(self);
//
// On the page:
//   import { workerClient } from "@ratmath/linalg/worker";
//   const la = workerClient(new Worker("linalg-worker.js", { type: "module" }));
//   const job = la.solve(A, b);   // job.quick is ready now
//   job.settle({ onStage, signal }).then(...)  // exact arrives from the worker
//
// Matrices travel as strings, so anything Matrix.from accepts may be passed.

import * as ops from "../index.js";
import { toPlain, fromPlain } from "./serialize.js";

const OPERATIONS = [
  "det", "rank", "rref", "solve", "inverse", "nullspace", "charpoly", "eigenvalues", "eigenvectors",
  "singularValues", "svd", "lu", "ldl", "qr", "gramSchmidt", "power", "pseudoinverse", "leastSquares",
  "similarity", "jordan", "expm", "sinm", "cosm", "sinhm", "coshm", "sqrtm", "logm", "powm",
];

/** Answer requests { id, op, args } with { id, stage } messages, and { id, cancel } to stop. */
export function serveLinalg(scope) {
  const running = new Map();
  scope.addEventListener("message", async (event) => {
    const { id, op, args, cancel } = event.data ?? {};
    if (cancel) {
      running.get(id)?.abort();
      return;
    }
    if (!OPERATIONS.includes(op)) {
      scope.postMessage({ id, error: `Unknown operation ${op}` });
      return;
    }
    const controller = new AbortController();
    running.set(id, controller);
    try {
      const job = ops[op](...fromPlain(args));
      await job.settle({
        signal: controller.signal,
        onProgress: (p) => scope.postMessage({ id, progress: toPlain(p) }),
        onStage: (s) => s.exact && scope.postMessage({ id, stage: toPlain(s) }),
      });
    } catch (e) {
      if (e?.name !== "AbortError") scope.postMessage({ id, error: String(e?.message ?? e) });
    } finally {
      running.delete(id);
    }
  });
}

/**
 * Same operations as the library, but each returns a Refinement-like object
 * whose `quick` is computed locally and whose exact stage comes from the worker.
 */
export function workerClient(worker) {
  let next = 1;
  const waiting = new Map();
  worker.addEventListener("message", (event) => {
    const { id, stage, progress, error } = event.data ?? {};
    const w = waiting.get(id);
    if (!w) return;
    if (progress) w.onProgress?.(progress);
    else if (error) {
      waiting.delete(id);
      w.reject(new Error(error));
    } else if (stage) {
      waiting.delete(id);
      w.resolve(fromPlain(stage));
    }
  });
  const client = {};
  for (const op of OPERATIONS) {
    client[op] = (...args) => {
      const local = ops[op](...args);
      let exactStage = null;
      return {
        get quick() {
          return local.quick;
        },
        get settled() {
          return exactStage !== null;
        },
        settle({ signal, onStage, onProgress } = {}) {
          onStage?.(local.quick);
          if (exactStage) {
            onStage?.(exactStage);
            return Promise.resolve(exactStage.value);
          }
          const id = next++;
          return new Promise((resolve, reject) => {
            const abort = () => {
              worker.postMessage({ id, cancel: true });
              waiting.delete(id);
              const e = new Error("Refinement aborted");
              e.name = "AbortError";
              reject(e);
            };
            if (signal?.aborted) return abort();
            signal?.addEventListener("abort", abort, { once: true });
            waiting.set(id, {
              onProgress,
              reject,
              resolve: (s) => {
                exactStage = s;
                onStage?.(s);
                resolve(s.value);
              },
            });
            worker.postMessage({ id, op, args: toPlain(args.map((a) => (Array.isArray(a) ? ops.Matrix.from(a) : a))) });
          });
        },
      };
    };
  }
  return client;
}

export { toPlain, fromPlain };
