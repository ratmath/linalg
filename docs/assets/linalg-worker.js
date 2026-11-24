// The docs' worker: answers exact computations off the page's main thread.
import { serveLinalg } from "./linalg.js?v=df9e75031d";

serveLinalg(self);
