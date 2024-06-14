module GCMC_2D_Fluids

# Write your package code here.
export GCMC_System, show, simulate, Histogram, GCMC_Simulation, try_insert!, try_delete!, try_move!, energy

include("structs.jl")
include("utils.jl")
include("sim_utils.jl")
include("simulation.jl")
end
