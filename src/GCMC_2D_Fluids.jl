module GCMC_2D_Fluids

# Write your package code here.
export GCMC_System, show, simulate, GCMC_Simulation, create_training_data, GCMC_TrainingData

include("structs.jl")
include("utils.jl")
include("sim_utils.jl")
include("simulation.jl")

end
