# lattice-gas
Lattice gas code to simulate diffusive processes in 2D
---
Features
* 1x1 Shape Tracer
* 2x2 Shape Tracer
* Output:
  * mean squared distance travelled
  * step rate
  * waiting time distribution
  * conditional probability for combinations of up to 3 steps
---
TODO:
* trapping
* step rejection logging (tracer type & relative position)
---
Build:
to build using the supplied makefile, do:

make build
---
Usage:
the program has to be supplied command line arguments, to use:

lattice_gas xdim ydim t_sim t_warm n_1x1 n_2x2 wtd_max wtd_res r_1x1 r_2x2 n_l

* xdim    size of lattice in x direction
* ydim    size of lattice in y direction
* t_sim   time max of the simulation
* t_warm  time warmup before simulation (no logging during warmup)
* n_1x1   number of tracers 1x1
* n_2x2   number of tracers 2x2
* wtd_max max waiting time logged (e.g. if 20, the wtd histogram goes up to 20)
* wtd_res resolution of the waiting time histogram (bin_width = 1/wtd_res)
* r_1x1   step attempt rate 1x1
* r_2x2   step attempt rate 2x2
* n_l     number of lattices (parallel runs of the sim), which are averaged for results
