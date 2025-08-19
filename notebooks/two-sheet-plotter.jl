### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# ╔═╡ 3fa905de-2848-11f0-2de5-8b4eebc36f0f
# ╠═╡ show_logs = false
begin
    cd(joinpath(@__DIR__, ".."))
    using Pkg: Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    # 
    using CairoMakie
    using CairoMakie.Makie.PlotUtils
    # 
    using PlutoUI
    using LaTeXStrings
    using Parameters
    using HadronicLineshapes
end

# ╔═╡ 39797575-ff41-4eeb-91b1-fc2d4bbeb5c5
md"""
# Two-sheet Breit-Wigner amplitude. Professional plots

Particle scattering, particularly in the context of resonance phenomena, necessitates a rigorous exploration of the complex analysis of scattering amplitudes. This notebook systematically investigates the multivalued nature of the Breit-Wigner amplitude, attributed to the complex branch points associated with the Mandelstam variables. We elucidate the emergence of Riemann sheets, stemming from inherent square-root singularities, and delineate between the physical and unphysical regions in the complex energy plane. By analyzing the complex $s$-plane, we provide a comprehensive understanding of the interplay between singularities, poles, and branch points. This analytical approach, bolstered by visual representations, offers a robust framework for interpreting the intricate structures governing particle interactions.

See the [Review on Resonances](https://pdg.lbl.gov/2021/reviews/rpp2021-rev-resonances.pdf) for more details.

"""

# ╔═╡ d48c02ef-4b0b-4b69-b068-818dbda10119
(path = joinpath(@__DIR__, "..", "plots"); !isdir(path) && mkdir(path))

# ╔═╡ 98047efc-9326-4f8f-a92c-c26ee663d60d
begin
    @with_kw struct myBW{T}
        m::T
        Γ::T
        mth::T
    end
    sheetI(bw::myBW, s::Number) =
        1 / (bw.m^2 - s - 1im * bw.m * bw.Γ * sqrt(s - bw.mth^2) / sqrt(bw.m^2 - bw.mth^2))
    sheetII(bw::myBW, s::Number) =
        1 / (bw.m^2 - s + 1im * bw.m * bw.Γ * sqrt(s - bw.mth^2) / sqrt(bw.m^2 - bw.mth^2))
end

# ╔═╡ bc8e0cb2-608f-4282-9cb1-bfdb6cb18550
md"""
### Breit-Wigner parameters

The resonance is characterized by:
- **m**: resonance mass (GeV)
- **Γ**: resonance width (GeV) 
- **mth**: threshold mass (GeV) - the mass of the decay products
"""

# ╔═╡ 4a4f9272-b3d2-455e-8e1b-1598ab6bf241
bw = myBW(; m = 0.77, Γ = 0.15, mth = 2 * 0.13)

# ╔═╡ 6bb69d21-4470-4cb3-88c7-556e3ff2afc7
md"""
**Current parameters:**
- Mass: $(bw.m) GeV
- Width: $(bw.Γ) GeV  
- Threshold: $(bw.mth) GeV
- Threshold squared: $(round(bw.mth^2, digits=3)) GeV²
"""

# ╔═╡ 2e595ec1-de80-453a-8f7e-0789bc6c1231
begin
    xv = range(-0.2, 1.1, 31)
    yv = range(-0.3, 0.3, 61)
    sv = xv' .+ 1im .* yv
end;

# ╔═╡ 98db5124-17c7-433a-9c0d-acea8942cff1
begin
    yv_up = yv[div(length(yv), 2):end]
    yv_dn = yv[1:div(length(yv), 2)+1]
end

# ╔═╡ 3f1e3cb1-6a71-4820-b8cc-d4e7b8ba218f
set_theme!(transparency = false)

# ╔═╡ 3fe6c7f3-737d-411c-ac9d-3099df298089
function surface_and_wireframe!(ax, x, y, z; col = :blue, kw...)
    surface!(ax, x, y, z; colormap = cgrad([col, col]), kw...)
    wireframe!(ax, x, y, z; color = :black, transparency = false, linewidth = 0.5, kw...)
end

# ╔═╡ 87c5322e-68b3-4fbf-9550-bdd7674bdcff
md"""
## 3D visualization of two-sheet Breit-Wigner amplitude

The function below creates a professional 3D plot showing the imaginary part of the Breit-Wigner amplitude on both Riemann sheets. Sheet I (blue) corresponds to the physical sheet, while Sheet II (red) corresponds to the unphysical sheet. The branch cut is shown as a black line starting from the threshold point (green dot).
"""

# ╔═╡ 34576736-848a-4cae-8181-6a659f172f4a
function two_sheets_s_plane_3d(grid_position)
    ax = Axis3(grid_position; aspect = (1, 1, 0.5),
        perspectiveness = 0.5, elevation = 0.47, azimuth = -1.2,
        xlabel = L"\mathrm{Re}(s)",
        ylabel = L"\mathrm{Im}(s)",
        zlabel = L"\mathrm{Im}(A)",
        zlabelrotation = π * (0.55),
        xlabeloffset = 20,
        ylabeloffset = 20,
        xtickformat = xv -> latexstring.(string.(xv)),
        ytickformat = xv -> latexstring.(string.(xv)),
        ztickformat = xv -> latexstring.(string.(xv)),
    )
    hidedecorations!(ax, ticks = true, ticklabels = true, label = false, grid = true)

    # Calculate amplitudes for both sheets
    zv1 = sheetI.(bw |> Ref, xv' .+ 1im .* yv_up) .|> imag
    zv2 = sheetI.(bw |> Ref, xv' .+ 1im .* yv_dn) .|> imag

    # Clip values for better visualization
    clip(z) = abs(z) < 15 ? z : 15 * sign(z)
    zv1 = map(clip, zv1)
    zv2 = map(clip, zv2)

    # Add surfaces with wireframes
    surface_and_wireframe!(ax, xv, yv_dn, zv2'; col = :red, transparency = false)
    surface_and_wireframe!(ax, xv, yv_up, zv1'; col = RGB(0.2, 0.4, 1), transparency = false)

    # Add branch cut line
    xv_cut = range(bw.mth^2, -0.2, 101)
    yv_cut = zero.(xv_cut)
    zv_cut = zero.(xv_cut)
    lines!(Point3f.(zip(xv_cut, yv_cut, zv_cut .+ 0.1)), linewidth = 5, color = :black)

    # Add threshold point
    scatter!(Point3f[(bw.mth^2, 0, 0.1)], markersize = 23, color = :lime)

    # Add sheet labels
    sv = ["I", "II"]
    xy = Point3f[(0.5, 0.15, 5), (0.5, -0.15, 5)]
    scatter!(ax, xy, markersize = 60, color = :white)
    text!(xy; text = sv, fontsize = 25, align = (:center, :center))

    # Set z limits
    zlims!(-20, 20)
end

# ╔═╡ 769b7491-4b29-4d42-b9e6-ad5cd0bb72e9
let
    f = Figure(fontsize = 25, size = (800, 600))
    two_sheets_s_plane_3d(f[1, 1])
    try
        save(joinpath(@__DIR__, "..", "plots", "two_sheets_s_plane.pdf"), f)
        save(joinpath(@__DIR__, "..", "plots", "two_sheets_s_plane.png"), f)
    catch e
        @info "Could not save: $(e)"
    end
    f
end

# ╔═╡ Cell order:
# ╟─39797575-ff41-4eeb-91b1-fc2d4bbeb5c5
# ╠═3fa905de-2848-11f0-2de5-8b4eebc36f0f
# ╠═d48c02ef-4b0b-4b69-b068-818dbda10119
# ╠═98047efc-9326-4f8f-a92c-c26ee663d60d
# ╟─bc8e0cb2-608f-4282-9cb1-bfdb6cb18550
# ╠═4a4f9272-b3d2-455e-8e1b-1598ab6bf241
# ╠═6bb69d21-4470-4cb3-88c7-556e3ff2afc7
# ╠═2e595ec1-de80-453a-8f7e-0789bc6c1231
# ╠═98db5124-17c7-433a-9c0d-acea8942cff1
# ╠═3f1e3cb1-6a71-4820-b8cc-d4e7b8ba218f
# ╠═3fe6c7f3-737d-411c-ac9d-3099df298089
# ╟─87c5322e-68b3-4fbf-9550-bdd7674bdcff
# ╠═34576736-848a-4cae-8181-6a659f172f4a
# ╠═769b7491-4b29-4d42-b9e6-ad5cd0bb72e9
