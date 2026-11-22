// Fast first, exact last.
//
// A Refinement is a computation that can answer at several levels of
// sharpness. `quick` is a Float64 preview, computed synchronously the moment
// it is read. The exact answer is computed by a generator that yields after
// each small unit of work, so it can be run all at once (`exact()`), in time
// slices that keep a page responsive (`settle()`), or abandoned mid-way when a
// slider moves and the question changes (`live()`).
//
// This shape is meant to be shared by the RatMath libraries (linear algebra,
// graphics, analytic): anything with `quick`, `stages()`, `exact()` and
// `settle()` can be driven by the same `live()` loop.
//
// Stage records:
//   { stage: "float",    exact: false, value, trust }   trust: good | rough | unsure
//   { stage: "progress", exact: false, done, of, step }
//   { stage: "exact",    exact: true,  value, via }      via: names the method

const now = typeof performance !== "undefined" ? () => performance.now() : () => Date.now();

function breathe() {
  if (typeof MessageChannel !== "undefined") {
    return new Promise((resolve) => {
      const { port1, port2 } = new MessageChannel();
      port1.onmessage = () => {
        port1.close();
        port2.close();
        resolve();
      };
      port2.postMessage(0);
    });
  }
  return new Promise((resolve) => setTimeout(resolve, 0));
}

function abortError() {
  const e = new Error("Refinement aborted");
  e.name = "AbortError";
  return e;
}

export class Refinement {
  #quickFn;
  #exactGen;
  #quick;
  #exact;
  #hasExact = false;

  /**
   * @param {() => {value, trust}} quick       cheap float preview
   * @param {() => Generator<progress, {value, via}>} exact  exact computation
   */
  constructor(quick, exact) {
    this.#quickFn = quick;
    this.#exactGen = exact;
  }

  /** The float preview, computed on first read. */
  get quick() {
    if (!this.#quick) {
      const { value, trust } = this.#quickFn();
      this.#quick = Object.freeze({ stage: "float", exact: false, value, trust });
    }
    return this.#quick;
  }

  /** True once the exact answer is known (it is cached). */
  get settled() {
    return this.#hasExact;
  }

  /** Every stage in order: the float preview, progress ticks, the exact answer. */
  *stages() {
    yield this.quick;
    if (!this.#hasExact) {
      const gen = this.#exactGen();
      let step;
      while (!(step = gen.next()).done) yield { stage: "progress", exact: false, ...step.value };
      this.#exact = Object.freeze({ stage: "exact", exact: true, value: step.value.value, via: step.value.via });
      this.#hasExact = true;
    }
    yield this.#exact;
  }

  /** The exact answer's stage record, computed synchronously. */
  exactStage() {
    if (!this.#hasExact) for (const _ of this.stages());
    return this.#exact;
  }

  /** The exact answer, computed synchronously. */
  exact() {
    return this.exactStage().value;
  }

  /**
   * Run to the exact answer in slices of about `sliceMs`, giving the event
   * loop a turn between slices. Calls onStage for the float preview and the
   * exact answer, onProgress for ticks. Rejects with an AbortError if
   * `signal` aborts first.
   */
  async settle({ signal, sliceMs = 8, onStage, onProgress } = {}) {
    if (signal?.aborted) throw abortError();
    let start = now();
    for (const s of this.stages()) {
      if (signal?.aborted) throw abortError();
      if (s.stage === "progress") onProgress?.(s);
      else onStage?.(s);
      if (s.exact) return s.value;
      if (now() - start > sliceMs) {
        await breathe();
        if (signal?.aborted) throw abortError();
        start = now();
      }
    }
  }
}

/**
 * Wire a computation to a control that changes often (a slider, a drag).
 *
 *   const update = live((t) => solve(A(t), b), (stage) => draw(stage));
 *   slider.oninput = () => update(slider.value);
 *
 * Each update shows the float preview at once, then sharpens to the exact
 * answer in the background. A newer update cancels older unfinished work, so
 * only the latest question is ever finished. Exact work starts after the
 * input has been still for `idleMs`.
 *
 * Returns update(...args) (which returns the Refinement), with
 * update.cancel() to stop pending work.
 */
export function live(compute, onStage, { idleMs = 30, sliceMs = 8, onProgress, onError } = {}) {
  let controller = null;
  let timer = null;
  const cancel = () => {
    controller?.abort();
    controller = null;
    if (timer !== null) clearTimeout(timer);
    timer = null;
  };
  const update = (...args) => {
    cancel();
    const job = compute(...args);
    onStage(job.quick, job);
    const mine = (controller = new AbortController());
    const go = () => {
      timer = null;
      job
        .settle({ signal: mine.signal, sliceMs, onProgress, onStage: (s) => s.exact && onStage(s, job) })
        .catch((e) => {
          if (e?.name !== "AbortError") (onError ?? ((err) => { throw err; }))(e);
        });
    };
    if (idleMs > 0) timer = setTimeout(go, idleMs);
    else go();
    return job;
  };
  update.cancel = cancel;
  return update;
}
