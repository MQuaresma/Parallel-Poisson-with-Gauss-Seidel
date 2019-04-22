# Parallel-Poisson-with-Gauss-Seidel
Parallel implementation, using a [Red-Black](https://www.cs.cornell.edu/~bindel/class/cs5220-s10/slides/lec14.pdf) scheme, of an 
iterative approach to the [Poisson equation](https://en.wikipedia.org/wiki/Poisson%27s_equation) using the 
[Gauss-Seidel method](https://en.wikipedia.org/wiki/Gauss–Seidel_method).


## Goals
- [x] Barebones (sequential) [implementation](src/SEQ) with no Red-Black scheme
- [x] Parallel-ready [implementation](src/RB-SEQ) using Red-Black scheme
- [x] Actual parallel [implementation](src/RB-PAR) in a [distributed memory](https://en.wikipedia.org/wiki/Distributed_memory) paradigm
