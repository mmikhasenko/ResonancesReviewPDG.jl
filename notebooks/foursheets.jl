### A Pluto.jl notebook ###
# v0.19.22

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
	using WGLMakie
	using WGLMakie.Makie.PlotUtils
	# 
	using PlutoUI
	using LaTeXStrings
end

# ╔═╡ 4a4f9272-b3d2-455e-8e1b-1598ab6bf241
Base.sqrt(z,ϕ) = sqrt(z*cis(ϕ))*cis(-ϕ/2)

# ╔═╡ 6bb69d21-4470-4cb3-88c7-556e3ff2afc7
ω(ϵ) = sqrt(ϵ,-π) + sqrt(ϵ-1,-π)

# ╔═╡ 46453bf3-770a-487a-bcf3-c4a34423a4f1
begin
	ω_II(ϵ) = -sqrt(ϵ,-π) + sqrt(ϵ-1,-π)
	ω_III(ϵ) = sqrt(ϵ,-π) - sqrt(ϵ-1,-π)
	ω_IV(ϵ) = -sqrt(ϵ,-π) - sqrt(ϵ-1,-π)
end

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

# ╔═╡ 87c5322e-68b3-4fbf-9550-bdd7674bdcff
fig = let
	ym = 1.4
	x = range(-1,2, 20)
	y1 = range(-ym,-1e-5,10)
	y2 = range(1e-5,ym, 10)
	y = vcat(y1, y2)
	# 
	fig = Figure(fontsize = 25, resolution = (1000, 800))
	ax = Axis3(fig[1, 1]; aspect = (1, 1, 0.5),
		    perspectiveness = 0.5, elevation = 0.47, azimuth=1.0,
		    xlabel = L"\mathrm{Re}(s)",
    		ylabel = L"\mathrm{Im}(s)",
			zlabel = "", # -\mathrm{Re}(\omega(s))
			xtickformat = xv->latexstring.(string.(xv)),
			ytickformat = xv->latexstring.(string.(xv)),
			ztickformat = xv->latexstring.(string.(xv)))
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
	xv = range(0,2,101)
	yv = zero.(xv)
	zv = map(real ∘ ω, xv .- 1e-6im)
	lines!(Point3f.(zip(xv,yv,zv .+ 0.03)), linewidth=4, color=:lime)
	scatter!(Point3f[(0,0.03,0.03)], markersize=23, color=:lime)
	# 
	sv = ["21", "12", "11", "22"]
	xy = Point3f[(0.9,-0.1,-0.1),(1.9,-0.2,1.5),(1.1,1.1,0.1),(1.9,-0.1,-1.0)]
	scatter!(ax, xy, markersize=80, color=:white)
	annotations!(sv, xy, fontsize=30, align=(:center,:center))	
	#
	try 
		s = save("foursheets_s_plane.pdf", fig)
	catch e
		@info "Could not save: $e.f"
	end
	fig
end

# ╔═╡ cd2215e2-2332-465d-a3e3-24b07b99ee8f
begin
	f = Figure(fontsize = 25, resolution = (1000, 800))
	ax = Axis(f[1, 1], aspect=1,
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
	lines!(ax, Point2f.(reim.(ω.(range(0,2,61)))), linewidth=6, color=:lime)
	scatter!(ax, Point2f[(0,1)], markersize=30, color=:lime)
	# 
	limits!(ax, (-1.5,1.5), (-1.5,1.5))
	sv = ["21", "12", "11", "22"]
	xy = Point2f[(0.4,0.4),(0.4,-0.4),(1.1,1.1),(1.1,-1.1)]
	scatter!(ax, xy, markersize=80, color=:white)
	annotations!(sv, xy, fontsize=30, align=(:center,:center))
	#
	try 
		s = save("foursheets_omega_plane.pdf", fig)
	catch e
		@info "Could not save: $e.f"
	end
	fig
	f
end

# ╔═╡ Cell order:
# ╠═2f4dd2f2-3c25-11ee-2f9d-8d445e1bd0c7
# ╠═4a4f9272-b3d2-455e-8e1b-1598ab6bf241
# ╠═6bb69d21-4470-4cb3-88c7-556e3ff2afc7
# ╠═46453bf3-770a-487a-bcf3-c4a34423a4f1
# ╠═ceee5325-c223-4b94-9ad0-e1359efd52b4
# ╠═cd2215e2-2332-465d-a3e3-24b07b99ee8f
# ╠═87c5322e-68b3-4fbf-9550-bdd7674bdcff
