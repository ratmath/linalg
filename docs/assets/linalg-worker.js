// The docs' worker: answers exact computations off the page's main thread.
import { serveLinalg } from "./linalg.js";

serveLinalg(self);
