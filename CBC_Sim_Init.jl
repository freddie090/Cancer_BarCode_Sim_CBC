
# Cancer Barcode Simulation
# CBC.
# Freddie Whiting - 2021

# Initialize Simulation.

# Load in Euler's
e = Base.MathConstants.e

using Distributions
using DataFrames
using RCall
using CSV
using Base.Threads
using Dates

include("CBC_Sim_Structs.jl")
include("CBC_Sim_Functions.jl")
include("CBC_Sim_Data_Coll.jl")
include("CBC_Sim_Experiments.jl")
#include("CBC_Sim_Plotting.jl")
