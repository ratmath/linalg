// The docs' worker: answers exact computations off the page's main thread.
import { serveLinalg } from "./linalg.js?v=41dbefb424";

serveLinalg(self);
