### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 2f4dd2f2-3c25-11ee-2f9d-8d445e1bd0c7
# ╠═╡ show_logs = false
begin
	cd(joinpath(@__DIR__, ".."))
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using CairoMakie
	using CairoMakie.Makie.PlotUtils
	# 
	using PlutoUI
	using LaTeXStrings
end

# ╔═╡ 0c56f2a6-d146-4d73-afc4-ff6b4b00ca9e
md"""
# Four-sheeted amplitude. Plots

Particle scattering, particularly in the context of resonance phenomena, necessitates a rigorous exploration of the complex analysis of scattering amplitudes. This notebook systematically investigates the multivalued nature of this amplitude, attributed to the complex branch points associated with the Mandelstam variables. We elucidate the emergence of Riemann sheets, stemming from inherent square-root singularities, and delineate between the physical and unphysical regions in the complex energy plane. By transitioning our analysis from the complex $s$-plane to the $\omega$-plane, we provide a comprehensive understanding of the interplay between singularities, poles, and branch points. This analytical approach, bolstered by visual representations, offers a robust framework for interpreting the intricate structures governing particle interactions.

See the [Review on Resonances](https://pdg.lbl.gov/2021/reviews/rpp2021-rev-resonances.pdf) for more details.

"""

# ╔═╡ d48c02ef-4b0b-4b69-b068-818dbda10119
(path = joinpath(@__DIR__, "..", "plots"); !isdir(path) && mkdir(path))

# ╔═╡ 4a4f9272-b3d2-455e-8e1b-1598ab6bf241
Base.sqrt(z,ϕ) = sqrt(z*cis(ϕ))*cis(-ϕ/2)

# ╔═╡ 3f823f5a-1e14-4e83-963b-ebbeeebb0fb4
md"""
Defining the ω function for the complex energy plane.
"""

# ╔═╡ 6bb69d21-4470-4cb3-88c7-556e3ff2afc7
ω(ϵ) = sqrt(ϵ,-π) + sqrt(ϵ-1,-π)

# ╔═╡ ccbaf000-5f4c-446e-b5a8-6a963fc1795c
md"""
Defining the ω functions for the different Riemann sheets.
"""

# ╔═╡ 46453bf3-770a-487a-bcf3-c4a34423a4f1
begin
	ω_II(ϵ) = -sqrt(ϵ,-π) + sqrt(ϵ-1,-π)
	ω_III(ϵ) = sqrt(ϵ,-π) - sqrt(ϵ-1,-π)
	ω_IV(ϵ) = -sqrt(ϵ,-π) - sqrt(ϵ-1,-π)
end

# ╔═╡ f808b8b9-5dd3-48cc-959b-bbcab0f0d04c
md"""
Function to generate a semicircle in the complex plane for visualization purposes.
"""

# ╔═╡ ceee5325-c223-4b94-9ad0-e1359efd52b4
function semicirc(c0,ch) 
	vcat(
		range(c0,2,101) .+ 1e-7im,
		range(2,ch,11) .+ 1e-7im,
		ch .* cis.(range(0,2π,20)[2:end-1]),
		range(ch,2,11) .- 1e-7im,
		range(2,c0,101) .- 1e-7im
		)
end

# ╔═╡ 9d3dae8f-d380-4254-a976-9d03cb8bf876
md"""
## Plot of energy plane in 3d
"""

# ╔═╡ 87c5322e-68b3-4fbf-9550-bdd7674bdcff
# fig1 = let
function four_sheets_s_plane_3d(grid_position)
	ym = 1.4
	x = range(-1,2, 20)
	y1 = range(-ym,-1e-5,10)
	y2 = range(1e-5,ym, 10)
	y = vcat(y1, y2)
	# 
	ax = Axis3(grid_position; aspect = (1, 1, 0.5),
		    perspectiveness = 0.5, elevation = 0.47, azimuth=1.0,
		    xlabel = L"\mathrm{Re}(s)",
    		ylabel = L"\mathrm{Im}(s)",
			zlabel = "", # -\mathrm{Re}(\omega(s))
			xtickformat = xv->latexstring.(string.(xv)),
			ytickformat = xv->latexstring.(string.(xv)),
			ztickformat = xv->latexstring.(string.(xv)))
	# 
	translate!(ax.scene, (0,0,2.5))
	scale!(ax.scene, 1.2, 1.2, 1.2)
	# 
	function addsurface!(x,y, f; c=:green)
		f′ = real ∘ f ∘ complex
		z = -f′.(x,y')
		surface!(ax, x, y, z;
			transparency = false, colormap=PlotUtils.cgrad([c,c]))
		wireframe!(ax, x, y, z, color=:black)
	end
	addsurface!(x,y1,ω_IV; c=:white)
	addsurface!(x,y1,ω_II; c=:red)
	addsurface!(x,y1,ω_III; c=:orange)
	addsurface!(x,y1,ω; c=RGB(0.2,0.4,1))
	addsurface!(x,y2,ω; c=RGB(0.2,0.4,1))
	# 
	xv = range(-1,0,2)
	yv = zero.(xv)
	zv = map(real ∘ ω, xv .- 1e-6im)
	lines!(Point3f.(zip(xv,yv,zv .+ 0.03)), linewidth=5, color=:black)
	# 
	xv = range(0,2,101)
	yv = zero.(xv)
	zv = map(real ∘ ω, xv .- 1e-6im)
	lines!(Point3f.(zip(xv,yv,zv .+ 0.03)), linewidth=5, color=:lime)
	scatter!(Point3f[(0,0.03,0.03)], markersize=23, color=:lime)
	# 
	sv = ["21", "12", "11", "22"]
	xy = Point3f[(0.8,-0.1,-0.2),(1.9,-0.1,1.2),(1.1,1.1,0.1),(1.9,-0.1,-1.5)]
	scatter!(ax, xy, markersize=80, color=:white)
	annotations!(sv, xy, fontsize=30, align=(:center,:center))
end

# ╔═╡ 769b7491-4b29-4d42-b9e6-ad5cd0bb72e9
begin
	f = Figure(fontsize = 25, resolution = (1200, 600))
	four_sheets_s_plane_3d(f[1,1])
	try 
		save(joinpath(@__DIR__, "..", "plots", "foursheets_s_plane.pdf"), f)
		save(joinpath(@__DIR__, "..", "plots", "foursheets_s_plane.png"), f)
	catch e
		@info "Could not save: $e.f"
	end
	Box(f[1, 1], color = (:red, 0.2), strokewidth = 0)
	f
end

# ╔═╡ 59e98547-1bb6-424c-9be7-362f8aee4ab0
md"""
## Plot of ω plane with sheets colored
"""

# ╔═╡ cd2215e2-2332-465d-a3e3-24b07b99ee8f
function foursheets_omega_plane(grid_position)
	ax = Axis(grid_position, aspect=1,
		    xlabel = L"\mathrm{Re}(\omega)",
    		ylabel = L"\mathrm{Im}(\omega)",
			xtickformat = xv->latexstring.(string.(xv)),
			ytickformat = xv->latexstring.(string.(xv)))
	# plot(xlim=(-1.5,1.5),ylim=(-1.5,1.5), aspect_ratio=1)
	poly!(ax, Point2f.(reim.(ω_II.(semicirc(0,1e5)))), color=:red)
	poly!(ax, Point2f.(reim.(ω_III.(semicirc(0,1e5)))), color=:orange)
	poly!(ax, Point2f.(reim.(ω_IV.(semicirc(0,2)))), color=RGB(0.95,0.95,0.95))
	poly!(ax, Point2f.(reim.(ω.(semicirc(0,2)))), color=:blue)
	# 
	hlines!(ax, [0], linewidth=2, color=:black)
	vlines!(ax, [0], linewidth=2, color=:black)
	# 
	lines!(ax, Point2f.(reim.(ω.(range(-1,0,2)))), linewidth=6, color=:black)
	lines!(ax, Point2f.(reim.(ω.(range(0,2,61)))), linewidth=6, color=:lime)
	scatter!(ax, Point2f[(0,1)], markersize=30, color=:lime)
	# 
	limits!(ax, (-1.5,1.5), (-1.5,1.5))
	sv = ["21", "12", "11", "22"]
	xy = Point2f[(0.4,0.4),(0.4,-0.4),(1.1,1.1),(1.1,-1.1)]
	scatter!(ax, xy, markersize=90, color=:white)
	annotations!(sv, xy, fontsize=30, align=(:center,:center))
end

# ╔═╡ 95c399cc-09c3-46b2-a3d5-9d56592931b2
let
	f = Figure(fontsize = 25, resolution = (1500, 700))
	four_sheets_s_plane_3d(f[1,1])
	foursheets_omega_plane(f[1,2])
	# 
	# Box(f[1, 1], color = (:red, 0.2), strokewidth = 0)
	# Box(f[1, 2], color = (:red, 0.2), strokewidth = 0)
	# 
	colsize!(f.layout, 1, Relative(0.62))
	# colgap!(f.layout, -30)
	#
	try
		save(joinpath(@__DIR__, "..", "plots", "foursheets_maps.pdf"), f)
		save(joinpath(@__DIR__, "..", "plots", "foursheets_maps.png"), f)
	catch e
		@info "Could not save: $e.f"
	end
	f
end

# ╔═╡ Cell order:
# ╟─0c56f2a6-d146-4d73-afc4-ff6b4b00ca9e
# ╠═2f4dd2f2-3c25-11ee-2f9d-8d445e1bd0c7
# ╠═d48c02ef-4b0b-4b69-b068-818dbda10119
# ╠═4a4f9272-b3d2-455e-8e1b-1598ab6bf241
# ╟─3f823f5a-1e14-4e83-963b-ebbeeebb0fb4
# ╠═6bb69d21-4470-4cb3-88c7-556e3ff2afc7
# ╟─ccbaf000-5f4c-446e-b5a8-6a963fc1795c
# ╠═46453bf3-770a-487a-bcf3-c4a34423a4f1
# ╟─f808b8b9-5dd3-48cc-959b-bbcab0f0d04c
# ╠═ceee5325-c223-4b94-9ad0-e1359efd52b4
# ╟─9d3dae8f-d380-4254-a976-9d03cb8bf876
# ╠═87c5322e-68b3-4fbf-9550-bdd7674bdcff
# ╠═769b7491-4b29-4d42-b9e6-ad5cd0bb72e9
# ╟─59e98547-1bb6-424c-9be7-362f8aee4ab0
# ╠═cd2215e2-2332-465d-a3e3-24b07b99ee8f
# ╠═95c399cc-09c3-46b2-a3d5-9d56592931b2
