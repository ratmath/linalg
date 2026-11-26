// The docs' worker: answers exact computations off the page's main thread.
import { serveLinalg } from "./linalg.js?v=fed8d93433";

serveLinalg(self);
