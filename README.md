# Cancer BarCode (CBC) Simulation 

CBC is an agent-based model implemented in [julia](https://julialang.org/) that simulates the experimental evolution of uniquely barcoded cancer cells. 

## Packages and Versions

CBC requires the following julia packages: 
`Distributions`,`DataFrames`,`RCall`,`CSV`

CBC currently runs on [julia](https://julialang.org/) version >= 1.4.0

CBC can be loaded on a local or remote machine by running 
```julia
include CBC_Sim_Init.jl
```

## Examples

The main experiment simulation is run using the function `Run_Exp_save_output`

This function can be found in `CBC_Sim_Experiments.jl`

Ensure that a directory exists in the CBC project folder named `Outputs`

The following are a list of parameters used with the function:
* `N` 
* `b`
* `d`
* `p`
* `mu`
* `sig`
* `del`
* `R_real`
* `n_pulse`
* `Nmax`
* `N_seed`
* `t_CO`
* `t_DT`
* `Nsim`
* `Passage`
* `insta_kill`
* `lim_probs`
* `psi`
* `al` 

