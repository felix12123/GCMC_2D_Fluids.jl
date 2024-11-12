using Plots, DelimitedFiles, ProgressMeter

function plot_file(file::String)
	d = readdlm(file, ';')
	L = Int(sqrt(size(d, 1)))
	rho = reshape(d[:, 2], L, L)
	c1 = reshape(d[:, 1], L, L)
	V = c1 .- log.(rho)

   plt_r = heatmap(rho, aspect_ratio=1, title="Density")
	plt_c = heatmap(c1, aspect_ratio=1, title="c1")
	plt_V = heatmap(V, aspect_ratio=1, title="Potential")
	plt = plot(plt_r, plt_c, plt_V, layout=(1, 3), size=(1200, 400))
	
	return plt
end
function plot_data_folder(folder::String, plot_folder::String=folder*"/plots")
	files = filter(x->occursin(".dat", x), readdir(folder, join=true))
	if isempty(files)
		println("No data files found in folder $folder")
		return
	end
	if !isdir(plot_folder)
		mkdir(plot_folder)
	end

	folder_name = split(folder, "/")[end]
	@showprogress "Plotting folder $(folder_name)" for file in files
		plt = plot_file(file)
		savefig(plt, plot_folder*"/$(basename(file)[1:end-4]).png")
	end
	
end

