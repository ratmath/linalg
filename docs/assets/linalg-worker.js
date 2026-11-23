// The docs' worker: answers exact computations off the page's main thread.
import { serveLinalg } from "./linalg.js?v=91b64861e5";

serveLinalg(self);
