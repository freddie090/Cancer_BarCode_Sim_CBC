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
* `N` - The number of uniquely barcoded cells when the simulation begins.
* `b` - The birth rates of all cells.
* `d` - The death rates of all cells.
* `p` - The probability any given cell is resistant when the simulation begins. 
* `mu` - The probability of a sensitive cell acquiring the resistant phenotype per cell division. 
* `sig` - The probability of a resistant cell acquiring the senstiive phenotype per cell division. 
* `del` - The relative fitness cost of the resistant phenotype. 
* `R_real` - How the fitness cost is incurred. Takes the following values: 
  * `b` - 
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

