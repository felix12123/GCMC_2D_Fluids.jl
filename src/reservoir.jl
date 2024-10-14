# if we have created a reservoir of training data, we might want to filter out specific cofigurations

using Statistics, DelimitedFiles, ProgressMeter, Plots

function total_uncertainty_of_error(xs::Vector{<:Real}, σs::Vector{<:Real})
	σ_stat = std(xs) / sqrt(length(xs))
	σ_prop = sqrt(sum(σs .^ 2)) / length(xs)
	return sqrt(σ_stat^2 + σ_prop^2)
end



function make_res_filter(;rho_lims=(0.0,1.0), c1_lims=(-Inf, Inf), minmax_rho=(0.0, Inf), n::Real=0.05)
	function res_filter(rho, c1)
		if !(rho_lims[1] <= mean(rho) <= rho_lims[2])
			println("rejected beacause rho is ", mean(rho), " but should be in ", rho_lims)
			return false
		end

		if c1_lims != (-Inf, Inf) && !(c1_lims[1] <= mean(c1) <= c1_lims[2])
			println("rejected beacause c1 is ", mean(c1), " but should be in ", c1_lims)
			return false
		end

		if minmax_rho != (0.0, Inf)
			srho = sort(rho |> vec)
			N = length(srho) * n |> round |> Int
			minmax = mean(srho[end-N:end]) / mean(srho[1:N])
			if !(minmax_rho[1] <= minmax <= minmax_rho[2])
				println("rejected beacause minmax_rho is ", minmax, " but should be in ", minmax_rho)
				return false
			end
		end
		
		return true
	end
	return res_filter
end

function filter_reservoir(res_folder::String, filter_func::Function, out_folder::String, downscale::Int=1)
	files = filter(x -> occursin(".dat", x), readdir(res_folder, join=true))
	if !isdir(out_folder)
		mkdir(out_folder)
	end
	if !isdir(out_folder*"/uncertainty")
		mkdir(out_folder*"/uncertainty")
	end

	@showprogress for file in files
		data = readdlm(file, ';')
		number = collect(basename(file))[isnumeric.(collect(basename(file)))] |> String
		rho_std_0 = readdlm(res_folder * "/uncertainty/uncert" * number * ".dat", ';')
		c1_0  = data[:, 1]
		rho_0 = data[:, 2]
		L = Int(sqrt(length(c1_0)))
		rho_0 = reshape(rho_0, L, L)
		c1_0  = reshape(c1_0,  L, L)
		rho_std_0 = reshape(rho_std_0, L, L)

		# we need to get βμloc to recalculate c1 after we scaled down the density profile.
		# c1 = log(ρ) - βμloc
		# βμloc = log(ρ) - c1
		βμloc_0 = log.(rho_0) .- c1_0

		l = Int(L÷downscale)
		rho     = zeros(l, l)
		rho_std = zeros(l, l)
		βμloc   = zeros(l, l)

		rho_window = zeros(downscale * downscale)
		bmloc_window = zeros(downscale * downscale)
		for ij in Iterators.product(1:l, 1:l)
			i, j = ij

			# save the window of the original data that we want to mean
			rho_window = vec(rho_0[(i-1)*downscale+1:i*downscale, (j-1)*downscale+1:j*downscale])
			rho_std_window = vec(rho_std_0[(i-1)*downscale+1:i*downscale, (j-1)*downscale+1:j*downscale])
			bmloc_window = vec(βμloc_0[(i-1)*downscale+1:i*downscale, (j-1)*downscale+1:j*downscale])

			# filter out any Infs or NaNs of bmloc_window
			if any(isfinite.(bmloc_window))
				βμloc[i, j] = mean(bmloc_window[isfinite.(bmloc_window)])
			else
				βμloc[i, j] = -Inf
			end
			rho[i, j]     = mean(rho_window)
			rho_std[i, j] = total_uncertainty_of_error(rho_window, rho_std_window)
		end
		c1 = log.(rho) .- βμloc

		if filter_func(rho, c1)
			writedlm(out_folder*"/$(basename(file)[1:end-4]).dat", [vec(c1) vec(rho)], ';')
			writedlm(out_folder*"/uncertainty/uncert$(number).dat", vec(rho_std), ';')
		end
	end
	return out_folder
end



function asses_reservoir(folder::String)
	N = length(filter(x -> occursin(".dat", x), readdir(folder, join=true)))
	rho_files = folder .* "/data" .* string.(1:N) .* ".dat"
	std_files = folder .* "/uncertainty/uncert" .* string.(1:N) .* ".dat"
	
	avg_std = zeros(N)
	avg_rho = zeros(N)
	avg_c1  = zeros(N)

	for i in 1:N
		data = readdlm(rho_files[i], ';')
		std_data = readdlm(std_files[i], ';')
		avg_std[i] = mean(std_data[.!iszero.(std_data)])
		avg_rho[i] = mean(data[:, 2][.!iszero.(data[:, 2])])
		avg_c1[i]  = mean(data[:, 1][isfinite.(data[:, 1])])
	end

	std_hist = histogram(avg_std, title="std histogram", xlabel="std", ylabel="count", bins=15)
	rho_hist = histogram(avg_rho, title="rho histogram", xlabel="rho", ylabel="count", bins=15)
	c1_hist  = histogram(avg_c1, title="c1 histogram", xlabel="c1", ylabel="count", bins=15)
	plt = plot(std_hist, rho_hist, c1_hist, layout=(1, 3), size=(1200, 400))
	display(plt)

	println("average std: ", mean(avg_std), " ± ", std(avg_std), "\tstd_range: ", minimum(avg_std), " - ", maximum(avg_std))
	println("average rho: ", mean(avg_rho), " ± ", std(avg_rho), "\trho_range: ", minimum(avg_rho), " - ", maximum(avg_rho))
	println("average c1: ", mean(avg_c1), " ± ", std(avg_c1), "\tc1_range: ", minimum(avg_c1), " - ", maximum(avg_c1))
	return avg_std, avg_rho, avg_c1
end