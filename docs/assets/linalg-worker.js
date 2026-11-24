// The docs' worker: answers exact computations off the page's main thread.
import { serveLinalg } from "./linalg.js?v=268bc0bc7e";

serveLinalg(self);
