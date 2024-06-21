# include("src/GCMC_2D_Fluids.jl")
using Statistics, ImageFiltering, Plots, LinearAlgebra, PProf
import PlotlyJS: surface

include("src/potentials.jl")
include("src/structs.jl")
include("src/utils.jl")
include("src/sim_utils.jl")
include("src/simulation.jl")
include("src/training_data.jl")

# Create a system
L = 16
dx = 0.05

Vext_sin(x; n::Int, A::Number, φ::Number, L::Number) = A * sin(2π * x * n / L + φ)

Vext_lin(x; x1::Number, x2::Number, E1::Number, E2::Number) = x > x1 && x < x2 ? E1 + (x - x1) * (E2 - E1) / (x2 - x1) : 0

Vext_wall(x; xw::Number, L::Number) = x < xw || x > L - xw ? Inf : 0

function generate_Vext(L::Number; num_sin=4, num_lin=rand(1:5), wall=true)
    Avar = 1.0
    sin_parameters = []
    for n in 1:num_sin  # Generate random parameters for periodic sine functions with increasing frequency
        push!(sin_parameters, (n = n, A = randn() * Avar, φ = rand() * 2π, L = L))
    end
    Evar = 1.0
    lin_parameters = []
    for _ in 1:num_lin  # Generate random parameters for discontinuous linear segments
        push!(lin_parameters, (x1 = round(rand() * L, digits=2), x2 = round(rand() * L, digits=2), E1 = randn() * Evar, E2 = randn() * Evar))
    end
    xwmax = 1.0
    wall_params = (xw = ceil(rand() * xwmax, digits=2), L = L)  # Set a random wall width
    function (x)  # Return a method which evaluates a combination of all functions with the chosen parameters above
        result = 0.0
        for sin_params in sin_parameters
            result += Vext_sin(x; sin_params...)
        end
        for lin_params in lin_parameters
            result += Vext_lin(x; lin_params...)
        end
        if wall
            result += Vext_wall(x; wall_params...)
        end
        result
    end
end

function Vext_gen()
	Vext1 = generate_Vext(L)
	Vext2 = generate_Vext(L)

	Vext(x) = Vext1(x[1]) + Vext2(x[2])
	return Vext
end


to = GCMC_TrainingData("box_data2",
	L=L,
	dx=dx,
	Vext_generator=Vext_gen,
	steps=1e7,
	num_systems=50,
	repetitions=25,
	threads=(4,Threads.nthreads()-4),
	μ_range=(-5, 20),
)

create_training_data(to, true)

function convert_2D_to1D(folder::String, target_folder::String)
	files = filter(x -> occursin(r".*\.dat", x), readdir(folder, join=true))
	if !isdir(target_folder)
		mkdir(target_folder)
	end
	for file in files
		data = readdlm(file, ';')
		c1 = data[:,1]
		rho = data[:,2]
		
		L = Int(sqrt(length(rho)))
		
		rho = reshape(rho, (L, L))
		c1 = reshape(c1, (L, L))
		replace!(c1, -Inf => NaN)
		replace!(c1, Inf => NaN)
		# bmloc = log.(rho) .- c1 
		# plot(heatmap(c1), heatmap(rho), heatmap(bmloc), layout=(1,3), size=(1200,400)) |> display
		# sleep(2)
		# bmloc1D = mean(bmloc, dims=2) |> vec
		rho1D = mean(rho, dims=2) |> vec
		c11D = log.(rho1D) .- bmloc1D
		# c11D = mean(c1, dims=2) |> vec
		# plot(plot(c11D), plot(rho1D), plot(bmloc1D), layout=(1,3), size=(1200,400)) |> display
		# sleep(2)
		writedlm(joinpath(target_folder, basename(file)), [c11D rho1D], ';')
	end
end

convert_2D_to1D("1D_data", "1D_data_1D")

# convert_2D_to1D("test_folder", "test_folder_1D")


nothing;