using Plots
import Base: show
function show(io::IO, sys::GCMC_System) # FIXME does not work properly yet with periodic boundary conditions
	p = plot(legend=false, xlim=(-sys.σ/2, sys.L+sys.σ/2), ylim=(-sys.σ/2, sys.L+sys.σ/2), size=(500, 500))
	for pos in sys.positions
		scatter!(p, [pos[1]], [pos[2]], markercolor=:blue, markersize=sys.σ/sys.L*230)
		if pos[1] < sys.σ/2
			scatter!(p, [pos[1]+sys.L], [pos[2]], markercolor=:blue, markersize=sys.σ/sys.L*230)
		end
		if pos[1] > sys.L - sys.σ/2
			scatter!(p, [pos[1]-sys.L], [pos[2]], markercolor=:blue, markersize=sys.σ/sys.L*230)
		end
		if pos[2] < sys.σ/2
			scatter!(p, [pos[1]], [pos[2]+sys.L], markercolor=:blue, markersize=sys.σ/sys.L*230)
		end
		if pos[2] > sys.L - sys.σ/2
			scatter!(p, [pos[1]], [pos[2]-sys.L], markercolor=:blue, markersize=sys.σ/sys.L*230)
		end
	end
			# p = plot!(p, [0, sys.L, sys.L, 0, 0], [0, 0, sys.L, sys.L, 0], linewidth=2, linecolor=:black)
			p = plot!(p, [0-sys.σ/4, sys.L+sys.σ/4, sys.L+sys.σ/4, -sys.σ/4, -sys.σ/4], [-sys.σ/4, -sys.σ/4, sys.L+sys.σ/4, sys.L+sys.σ/4, -sys.σ/4], linewidth=sys.σ/sys.L*230, linecolor=:black)

	display(p)
end
