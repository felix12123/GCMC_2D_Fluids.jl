module GCMC_2D_Fluids


include("visualisation.jl")

include("structs.jl")
include("utils.jl")
include("potentials.jl")
include("sim_utils.jl")
include("simulation.jl")
include("training_data.jl")
include("training_data_1D.jl")
include("reservoir.jl")

# Write your package code here.
export GCMC_System,
	show,
	simulate,
	GCMC_Simulation,
	GCMC_TrainingData,
	create_training_data,
	create_training_data_1D,
	generate_random_potential,
	generate_random_potential_1D,
	PotentialOptions,
	get_c1,
	make_res_filter,
	filter_reservoir,
	plot_data_folder

end
