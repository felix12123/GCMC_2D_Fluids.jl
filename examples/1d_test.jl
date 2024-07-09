using Statistics

include("src/potentials.jl")
include("src/structs.jl")
include("src/utils.jl")
include("src/sim_utils.jl")
include("src/simulation.jl")
include("src/training_data.jl")
include("src/training_data_1D.jl")


L = 10
d = 1
dx = 0.03


accept_condition = function(rho)
	# check if the density is too low
	if mean(rho) < 0.1
		println("rejected because of mean, mean=$(mean(rho))")
		return false
	end
	# check if the density is too high
	if mean(rho)*pi/4 > 0.5
		println("rejected because of mean, mean=$(mean(rho))")
		return false
	end
	srho = filter(x -> isfinite(x) && !iszero(x), sort(rho |> vec))
	n = length(srho)÷20
	maxmin = mean(srho[end-n:end]) / mean(srho[1:n])
	if maxmin > 20
		println("rejected because of maxmin.maximin=$maxmin")
		return false
	end
	return true
end

Vext_gen = function()
	po = PotentialOptions(L=L, num_sin=(5, 10), sin_amp=(0, 0.5), sin_periods=(0.5, 4), periodic=false, wall=:walls)
	generate_random_potential_1D(po)
end

to = GCMC_TrainingData("1d_dx0_03",
	L=L,
	dx=dx,
	Vext_generator=Vext_gen,
	steps=2e7,
	num_systems=10,
	repetitions=250,
	# accept_condition=accept_condition,
	# rho_smooth_func=smooth_rho,
	μ_range=(-5, 5),
	threads=(2, Threads.nthreads()),
)


# create_training_data_1D(to, true)

function show_data_folder(folder::String, pause = 1)
	files = filter(x -> occursin(".dat", x), readdir(folder, join=true))
	for file in files
		d = readdlm(file, ';')
		L = Int(sqrt(length(d[:, 1])))
		c1 = reshape(d[:, 1], L, L)
		rho = reshape(d[:, 2], L, L)
		plot(heatmap(c1, title="c1"), heatmap(rho, title="rho"),plot_title=basename(file), size=(800, 400)) |> display
		# set overall title to filename
		sleep(pause)
	end
end

show_data_folder_1D("1d_dx0_03_1D")

function show_data_folder_1D(folder::String, pause = 1)
	files = filter(x -> occursin(".dat", x), readdir(folder, join=true))
	for file in files
		d = readdlm(file, ';')
		c1 = d[:, 1]
		rho = d[:, 2]
		c1 = mean(reshape(c1, (Int(sqrt(length(c1))), Int(sqrt(length(c1))))), dims=2)[:, 1]
		rho = mean(reshape(rho, (Int(sqrt(length(rho))), Int(sqrt(length(rho))))), dims=2)[:, 1]
		plot(plot(c1, title="c1"), plot(rho, title="rho"), size=(800, 400)) |> display
		sleep(pause)	
	end
end

using ImageFiltering, Plots, LinearAlgebra
function inhomogenity(rho)
	if length(size(rho)) == 1
		rho = reshape(rho, (Int(sqrt(length(rho))), Int(sqrt(length(rho)))))
	end

	kernel = [0.0 1.0 0.0; 1.0 -4.0 1.0; 0.0 1.0 0.0]
	return imfilter(rho ./ mean(rho), kernel) .|> abs #|> mean
end

function my_runstd(x::Matrix, rad::Int)
	res = zeros(size(x))
	for ij in Iterators.product(1:size(x, 1), 1:size(x, 2))
		i, j = ij
		res[i, j] = std(view(x, max(1, i-rad):min(size(x, 1), i+rad), max(1, j-rad):min(size(x, 2), j+rad)))
	end
	res
end


files = filter(x -> occursin(".dat", x), readdir("1d_dx0_03", join=true))
for file in files
	d = readdlm(file, ';')
	c1 = d[:, 1]
	rho = d[:, 2]
	L = Int(sqrt(length(c1)))
	c1 = reshape(c1, L, L)
	rho = reshape(rho, L, L)
	if mean(rho) == 0
		rm(file)
	end
end


# i = 3
# d = readdlm("1d_dx0_03/data$i.dat", ';')
# c1 = d[:, 1]
# rho = d[:, 2]
# L = Int(sqrt(length(c1)))
# c1 = reshape(c1, L, L)
# rho = reshape(rho, L, L)
# rho_i = inhomogenity(rho)
# rho_s = my_runstd(rho, 2)
# plot(heatmap(rho, title="rho"), heatmap(rho_s[3:end-3, 3:end-3], title="local std"), size=(800, 400)) |> display
