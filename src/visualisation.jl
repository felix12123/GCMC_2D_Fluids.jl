function plot_data_folder(folder::String, plot_folder::String=folder*"/plots"; project=:yes)
	files = filter(x->occursin(".dat", x), readdir(folder, join=true))
	if isempty(files)
		println("No data files found in folder $folder")
		return
	end
	if !isdir(plot_folder)
		mkdir(plot_folder)
	end

	function plot_rho(rho)
		if project == :yes
			rho_1d = mean(rho, dims=2)
			return plot(rho_1d, xlabel="x", ylabel="ρ(x)", legend=false)
		else
			return heatmap(rho, aspect_ratio=1)
		end
	end
	
	for file in files
		d = readdlm(file, ';')
		L = Int(sqrt(size(d, 1)))
		rho = reshape(d[:, 2], L, L)
		c1 = reshape(d[:, 1], L, L)
		plt_r = plot_rho(rho)
		plt_c = plot_rho(c1)
		plt = plot(plt_r, plt_c, layout=(1, 2), size=(800, 400))
		savefig(plt, plot_folder*"/$(basename(file)[1:end-4]).png")
	end
	
end

