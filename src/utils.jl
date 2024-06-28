using Plots
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
