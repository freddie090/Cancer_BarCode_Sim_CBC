# Cancer BarCode (CBC) Simulation 

CBC is an agent-based model implemented in [julia](https://julialang.org/) that simulates the experimental evolution of uniquely barcoded cancer cells. 

## Packages and Versions

CBC requires the following julia packages: 
`Distributions`,`DataFrames`,`RCall`,`CSV`

CBC currently runs on [julia](https://julialang.org/) version >= 1.4.0

n.b. that because CBC requires `RCall`, it also requires a version of [R](https://www.r-project.org/). 
Details on how to get `RCall` working with `julia` can be found [here](https://juliainterop.github.io/RCall.jl/stable/installation/).

CBC can be loaded on a local or remote machine by downloading the 6 CBC `.jl` files and running 
```julia
include CBC_Sim_Init.jl
```
in julia from within the CBC directory. 

## Simulation Design 

The simulation recapitulates key features of an evolutionary experiment as follows: 
1. Cells are tagged with a unique molecular 'barcode' that is completely heritable.
2. All cells share a mutual expansion step. 
3. The expanded cells are split into 8 replicate sub-populations. 
   * 4 are control (`CO`) replicates
   * 4 are drug-treatment (`DT`) replicates
4. The sub-populations are grown until they reach a given population size. 
   * The control replicates are simply grown with their assigned birth and death rates. 
   * The drug-treatment replicates are periodically exposed to drug-treatment which kills cells given a resistant phenotype (`R`).
5. This process is repeated for a given number of passages, where cells from the previous passage are used to seed the subsequent.
6. The number of cells and the distribution of barcode lineages are stored for each replicate sub-population for statistical analysis.

Various parameters control the experimental setup and the evolution of a resistant phenotype.

More details and an example of how to run the experiment function are given in the Examples section.


## Examples

The main experiment simulation is run using the function `Run_Exp_save_output`

This function can be found in `CBC_Sim_Experiments.jl`

Ensure that a directory exists in the CBC project folder named `Outputs`

The following are a list of parameters used with the function:
* `N` - the number of uniquely barcoded cells when the simulation begins.
* `b` - the birth rates of all cells (currently can only take one value).
* `d` - the death rates of all cells (currently can only take one value).
* `p` - the probability any given cell is resistant when the simulation begins. 
* `mu` - the probability of a sensitive cell acquiring the resistant phenotype per cell division. 
* `sig` - the probability of a resistant cell acquiring the senstiive phenotype per cell division. 
* `del` - the relative fitness cost of the resistant phenotype. 
* `R_real` - how the fitness cost is incurred. Can take one of the following values: 
  * `b` - the cost is realised as a lower birth date (default). 
  * `d` - the cost is realised as a higher death rate.
  * `l` - the cost is realised as a lower net growth rate (`b - d`). 
* `n_pulse` - the number of times to 'pulse' the drug-treatment replicates with a drug-killing step. 
* `Nmax` - the maximum population size for each replicate sub-population. 
* `N_seed` - the number of cells used to seed reach replicate sub-population. 
* `t_CO` - the maximum time the control replicates can run for. 
* `t_DT` - the maximum time the drug-treatment replicates can run for. 
* `Nsim` - the number of simulations to run for the current parameter set. 
* `Passage` - the number of passages replicate sub-populations will experience.
* `insta_kill` - should cells be killed instantly after assignment to drug-treatment replicates, or allowed to grow for 1x `t_DT`.
* `lim_probs` - should the resistant phenotype be assigned to cells when the simulation begins according to the parameter equilibrium frequencies? 

The variable types for each function argument are as follows: 

```julia 
function Run_Exp_save_output(N::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64,
    R_real::String, n_pulse::Int64, Nmax::Int64,
    N_seed::Int64, t_CO::Float64, t_DT::Float64, Nsim::Int64, Passage::Int64,
    insta_kill::Bool, lim_probs::Bool)
```

For example, the following runs a given simulation with the specified parameter values in julia: 

```julia 
Run_Exp_save_output(1000000, 0.8, 0.2, 0.0, 10^-4, 10^-4, 0.9, "l", 60, 6400000, 100000, 120.0, 120.0, 1, 4, true, true)
```


