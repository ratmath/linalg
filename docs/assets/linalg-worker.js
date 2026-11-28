// The docs' worker: answers exact computations off the page's main thread.
import { serveLinalg } from "./linalg.js?v=5591b17725";

serveLinalg(self);
