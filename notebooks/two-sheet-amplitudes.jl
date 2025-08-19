### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 3fa905de-2848-11f0-2de5-8b4eebc36f0f
# ╠═╡ show_logs = false
begin
	cd(joinpath(@__DIR__, ".."))
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using CairoMakie
	# using WGLMakie
	using CairoMakie.Makie.PlotUtils
	# 
	using PlutoUI
	using LaTeXStrings
	using Parameters
	using HadronicLineshapes
end

# ╔═╡ 39797575-ff41-4eeb-91b1-fc2d4bbeb5c5
md"""
# Two sheet amplitude. Plots

Particle scattering, particularly in the context of resonance phenomena, necessitates a rigorous exploration of the complex analysis of scattering amplitudes. This notebook systematically investigates the multivalued nature of this amplitude, attributed to the complex branch points associated with the Mandelstam variables. We elucidate the emergence of Riemann sheets, stemming from inherent square-root singularities, and delineate between the physical and unphysical regions in the complex energy plane. By transitioning our analysis from the complex $s$-plane to the $\omega$-plane, we provide a comprehensive understanding of the interplay between singularities, poles, and branch points. This analytical approach, bolstered by visual representations, offers a robust framework for interpreting the intricate structures governing particle interactions.

See the [Review on Resonances](https://pdg.lbl.gov/2021/reviews/rpp2021-rev-resonances.pdf) for more details.

"""

# ╔═╡ 98047efc-9326-4f8f-a92c-c26ee663d60d
begin
	@with_kw struct myBW{T}
		m::T
		Γ::T
		mth::T
	end
	sheetI(bw::myBW, s::Number) =
		1/(bw.m^2-s-1im*bw.m*bw.Γ*sqrt(s-bw.mth^2)/sqrt(bw.m^2-bw.mth^2))
	sheetII(bw::myBW, s::Number) =
		1/(bw.m^2-s+1im*bw.m*bw.Γ*sqrt(s-bw.mth^2)/sqrt(bw.m^2-bw.mth^2))
end

# ╔═╡ bc8e0cb2-608f-4282-9cb1-bfdb6cb18550
bw = myBW(; m=0.77, Γ=0.15, mth=2*0.13)

# ╔═╡ 2e595ec1-de80-453a-8f7e-0789bc6c1231
begin
	xv = range(-0.2, 1.1, 31)
	yv = range(-0.3, 0.3, 61)
	sv = xv' .+ 1im .* yv
end;

# ╔═╡ 98db5124-17c7-433a-9c0d-acea8942cff1
begin
	yv_up = yv[div(length(yv),2):end]
	yv_dn = yv[1:div(length(yv),2)+1]
end

# ╔═╡ 3f1e3cb1-6a71-4820-b8cc-d4e7b8ba218f
set_theme!(transparency=false)

# ╔═╡ 3fe6c7f3-737d-411c-ac9d-3099df298089
function surface_and_wireframe!(ax, x, y, z; col=:blue, kw...)
	surface!(ax, x, y, z; colormap=cgrad([col, col]), kw...)
	# wireframe!(ax, x, y, z; color=:black, transparency=false, linewidth=0.1, kw...)
end

# ╔═╡ 34576736-848a-4cae-8181-6a659f172f4a
let
    fig = Figure(size=(700, 500))
    ax = Axis3(fig[1, 1];
        aspect = (1, 1, 0.6),
        perspectiveness = 0.5,
        elevation = 0.5, azimuth = -0.9
    )
	

	zv1 = sheetI.(bw |> Ref, xv' .+ 1im .* yv_up) .|> imag
	zv2 = sheetII.(bw |> Ref, xv' .+ 1im .* yv_dn) .|> imag
	# 
    surface_and_wireframe!(ax, xv, yv_dn, zv2'; col=:red)
    surface_and_wireframe!(ax, xv, yv_up, zv1'; col=:red)
	# 
	zlims!(-20,20)
	# 
    fig
end

# ╔═╡ df668f80-18b0-40ac-b9b3-6bb634954298
let
    fig = Figure(size=(700, 500))
    ax = Axis3(fig[1, 1];
        aspect = (1, 1, 0.6),
        perspectiveness = 0.5,
        elevation = 0.5, azimuth = -0.9
    )
	

	zv1 = sheetI.(bw |> Ref, xv' .+ 1im .* yv_up) .|> imag
	zv2 = sheetI.(bw |> Ref, xv' .+ 1im .* yv_dn) .|> imag
	clip(z) = abs(z)<15 ? z : 15*sign(z)
	zv1 = map(clip, zv1)
	zv2 = map(clip, zv2)
	# 
    surface_and_wireframe!(ax, xv, yv_up, zv1'; col=:red, transparency=false)
    surface_and_wireframe!(ax, xv, yv_dn, zv2'; col=:green, transparency=false)
	# 
	zlims!(-20,20)
	# 
    fig
end

# ╔═╡ Cell order:
# ╟─39797575-ff41-4eeb-91b1-fc2d4bbeb5c5
# ╠═3fa905de-2848-11f0-2de5-8b4eebc36f0f
# ╠═98047efc-9326-4f8f-a92c-c26ee663d60d
# ╠═bc8e0cb2-608f-4282-9cb1-bfdb6cb18550
# ╠═2e595ec1-de80-453a-8f7e-0789bc6c1231
# ╠═98db5124-17c7-433a-9c0d-acea8942cff1
# ╠═3f1e3cb1-6a71-4820-b8cc-d4e7b8ba218f
# ╠═3fe6c7f3-737d-411c-ac9d-3099df298089
# ╠═34576736-848a-4cae-8181-6a659f172f4a
# ╠═df668f80-18b0-40ac-b9b3-6bb634954298
