using Plots, ImageFiltering, DelimitedFiles
import Base: show
function show(io::IO, sys::GCMC_System) # FIXME does not work properly yet with periodic boundary conditions
	ms = sys.σ/sys.L*220
	p = plot(legend=false, xlim=(-sys.σ/2, sys.L+sys.σ/2), ylim=(-sys.σ/2, sys.L+sys.σ/2), size=(500, 500))
	for pos in reverse.(sys.positions[1:sys.N, :] |> eachrow)
		scatter!(p, [pos[1]], [pos[2]], markercolor=:blue, markersize=ms)
		if pos[1] < sys.σ/2
			scatter!(p, [pos[1]+sys.L], [pos[2]], markercolor=:blue, markersize=ms)
		end
		if pos[1] > sys.L - sys.σ/2
			scatter!(p, [pos[1]-sys.L], [pos[2]], markercolor=:blue, markersize=ms)
		end
		if pos[2] < sys.σ/2
			scatter!(p, [pos[1]], [pos[2]+sys.L], markercolor=:blue, markersize=ms)
		end
		if pos[2] > sys.L - sys.σ/2
			scatter!(p, [pos[1]], [pos[2]-sys.L], markercolor=:blue, markersize=ms)
		end
	end
	# draw box
	p = plot!(p, [-sys.σ/2, sys.L+sys.σ/4, sys.L+sys.σ/4, -sys.σ/4, -sys.σ/4], [-sys.σ/4, -sys.σ/4, sys.L+sys.σ/4, sys.L+sys.σ/4, -sys.σ/2], linewidth=sys.σ/sys.L*230, linecolor=:black)

	display(p)
end

function my_runstd(x::Matrix, rad::Int)
	res = zeros(size(x))
	for ij in Iterators.product(1:size(x, 1), 1:size(x, 2))
		i, j = ij
		res[i, j] = std(view(x, max(1, i-rad):min(size(x, 1), i+rad), max(1, j-rad):min(size(x, 2), j+rad)))
	end
	res
end
function my_runmean(x::Matrix, rad::Int)
	res = zeros(size(x))
	for ij in Iterators.product(1:size(x, 1), 1:size(x, 2))
		i, j = ij
		res[i, j] = mean(view(x, max(1, i-rad):min(size(x, 1), i+rad), max(1, j-rad):min(size(x, 2), j+rad)))
	end
	res
end
function smooth_rho(rho, rad::Int)
	# smooth the density
	rho = imfilter(rho, Kernel.gaussian(rad))
	return rho
end
function analyze_data_folder(folder::String)
	files = filter(x -> contains(x, ".dat"), readdir(folder, join=true))
	get_noise(rh0) = mean(my_runstd(rh0 .- smooth_rho(rh0, 5), 1))
	function minmax(rho)
		srho = sort(vec(rho))
		N = length(srho) * 0.05 |> round |> Int
		return mean(srho[end-N:end]) / mean(srho[1:N])
	end
	
	noises = []
	means = []
	minmaxs = []
	for file in files
		d = readdlm(file, ';')
		rho = d[:, 2]
		L = Int(sqrt(length(rho)))
		rho = reshape(rho, L, L)

		push!(means, mean(rho))
		push!(minmaxs, minmax(rho))
		push!(noises, get_noise(rho ./ mean(rho)))
	end
	
	println("Mean = ", mean(means), " ± ", std(means))
	println("Minmax = ", mean(minmaxs), " ± ", std(minmaxs))
	println("Noise = ", mean(noises), " ± ", std(noises))

	plt1 = histogram(means, bins=20, label="mean")
	plt2 = histogram(minmaxs, bins=20, label="minmax")
	plt3 = histogram(noises, bins=20, label="noise")

	plot(plt1, plt2, plt3, layout=(3, 1), size=(500, 1000))
end