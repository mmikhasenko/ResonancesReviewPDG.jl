### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ d4a12a4c-480b-4323-be97-49ba4d40d817
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	# 
	Pkg.add("Plots")
	Pkg.add("PlotlyJS")
	Pkg.add("Parameters")
	# 
	using Plots
	using Parameters
	# 
	plotlyjs()
end

# ╔═╡ 5d1282ae-70f1-4fcd-b0cc-f97afbb2fc2f
md"""
## First and second sheet
"""

# ╔═╡ a3768c8d-3a16-4673-a8db-a817183c83f5
plotlyjs()

# ╔═╡ a9bdb1b6-9a0a-4102-a1e0-d959b7a5f398
theme(:boxed, lab="", fontfamily="Computer Modern")

# ╔═╡ 51440b8d-7a1e-4fd8-9b8e-1fdd0e13ff24
const iϵ = 1e-6im

# ╔═╡ 8a3ff5f9-dcea-4c7d-9099-c0b5c19fdc53
begin
	@with_kw struct myBW{T}
		m::T
		Γ::T
		mth::T
	end
	sheetI(bw::myBW, s::Number) =
		1/(bw.m^2-s-1im*bw.m*bw.Γ*kofs(s, bw.mth)/kofs(bw.m^2, bw.mth))
	sheetII(bw::myBW, s::Number) =
		1/(bw.m^2-s+1im*bw.m*bw.Γ*kofs(s, bw.mth)/kofs(bw.m^2, bw.mth))
	# 
	ofk(bw::myBW, k) =
		1/(bw.m^2-sofk(k, bw.mth)-1im*bw.m*bw.Γ*k/kofs(bw.m^2, bw.mth))
	sofk(k, mth) = k^2 + mth^2
	kofs(s, mth) = sqrt(s- mth^2)
end

# ╔═╡ dc07efbe-18cb-4c76-b94d-b2c5465a0e63
bw = myBW(; m=0.5, Γ=0.2, mth=2*0.13)

# ╔═╡ b818c3d0-30ed-4306-a654-940cca662a9f
k_pole = let
	# k^2+k*1im*bw.m*bw.Γ/sqrt(bw.m^2-bw.mth^2)+(mth^2-bw.m^2)
	a = 1
	b = 1im*bw.m*bw.Γ/sqrt(bw.m^2-bw.mth^2)
	c = (bw.mth^2-bw.m^2)
	D = b^2-4a*c
	(-b .+ [-1,1] .* sqrt(D)) ./ (2a)
end

# ╔═╡ cdcfbffc-248d-4529-b5c8-ab8cbb1201f2
const branch_point = (bw.mth^2, 0.0, imag(sheetI(bw, bw.mth^2+iϵ))+0.7);

# ╔═╡ 9b80d058-dd66-4084-91dd-2f1305fe010e
begin
	xv = range(-0.05, 0.45, 151)
	yv = range(-0.25, 0.25, 161)
	sv = xv' .+ 1im .* yv
end;

# ╔═╡ ad238351-b312-488e-83a8-fff1cd246d83
begin
	yv_up = yv[div(length(yv),2):end]
	yv_dn = yv[1:div(length(yv),2)]
end;

# ╔═╡ 453279b4-8855-40fc-b5a9-61510ce73f9c
md"""
## k-plane
"""

# ╔═╡ e7f32f79-e857-4b12-96da-e9b308191ca0
begin
	xkv = range(-0.7, 0.7, 100)
	ykv = range(-0.7, 0.7, 100)
end;

# ╔═╡ bc542219-f869-4379-9734-f6ebadcf6b57
md"""
## Surface
"""

# ╔═╡ 8aefc2d3-7fa3-413e-8020-0539e1cc663b
begin
	palette = :wong
	# 
	sheetI_color = theme_palette(palette)[5]
	sheetI_grad = cgrad([sheetI_color,sheetI_color*1.1])
	# 
	sheetII_color = theme_palette(palette)[1]
	sheetII_grad = cgrad([sheetII_color,sheetII_color*1.1])
	# 
	unitarity_color = theme_palette(palette)[10]
end;

# ╔═╡ 0a4d5bb9-bed9-498b-aeb1-7aa5d38e4729
let
	plot(colorbar=false)
	f(bw, z) = sheetI(bw, z)
	zv = f.(bw |> Ref, xv' .+ 1im .* yv) .|> abs .|> log10 .|> log10
	contour!(xv, yv, zv, c=cgrad(:roma, rev=true), levels=20)
	# 
	plot!([bw.mth^2, xv[1]] .+ 0.0im, lw=3, lc=:black)
	scatter!([bw.mth^2] .+ 0.0im, ms=6, mc=unitarity_color)
	# 
	plot!(aspect_ratio=1)
end

# ╔═╡ 679c5c98-d055-4817-a548-75fb4d8eb158
let
	plot(colorbar=false)
	f(bw, z) = imag(z) > 0 ?  sheetI(bw, z) : sheetII(bw, z)
	zv = f.(bw |> Ref, xv' .+ 1im .* yv) .|> abs .|> log10 .|> log10
	contour!(xv, yv, zv, c=cgrad(:roma, rev=true), levels=20)
	# 
	plot!([bw.mth^2, xv[end]] .+ 0.0im, lw=3, lc=unitarity_color)
	scatter!([bw.mth^2] .+ 0.0im, ms=6, mc=unitarity_color)
	plot!(aspect_ratio=1)
end

# ╔═╡ 407291e4-23df-4c9b-b64b-c9e4a9af9720
let
	plot(colorbar=false)
	f(bw, k) = ofk(bw, k)
	zv = f.(bw |> Ref, xkv' .+ 1im .* ykv) .|> abs .|> log10
	plot!([xkv[1],xkv[end]], [ykv[end],ykv[end]], fill=0,c=sheetI_color)
	plot!([xkv[1],xkv[end]], [ykv[1],ykv[1]], fill=0, c=sheetII_color)
	# 
	#
	for x in xv[1:10:end]
		plot!(kofs.(x .+ yv_up[2:end] .* 1im, bw.mth), lc=:black, lw=0.4)
		plot!(kofs.(x .+ yv_dn .* 1im, bw.mth), lc=:black, lw=0.4)
		plot!(-kofs.(x .+ yv_up[2:end] .* 1im, bw.mth), lc=:black, lw=0.4)
		plot!(-kofs.(x .+ yv_dn .* 1im, bw.mth), lc=:black, lw=0.4)
	end
	for y in yv[1:10:end]
		plot!(kofs.(y .* 1im .+ xv[2:end], bw.mth), lc=:black, lw=0.4)
		plot!(-kofs.(y .* 1im .+ xv[2:end], bw.mth), lc=:black, lw=0.4)
	end
	# 
	scatter!(k_pole, m=(:x, 6, :red))
	plot!([xkv[end],0.0im, ykv[end]*1im], lw=5, lc=unitarity_color)
	plot!([0, ykv[end]*1im], lw=5, lc=:black)
	scatter!([0.0im], ms=8, mc=unitarity_color)
	# 
	plot!(aspect_ratio=1, xlab="Re k", ylab="Im k")
	# 
	savefig(joinpath("k-plane.png"))
	savefig(joinpath("k-plane.pdf"))
	# 
	plot!()
end

# ╔═╡ 6d14b882-9469-4ed8-84b0-25d72df3c02d
theme_palette(palette)[10];

# ╔═╡ 3d659a6b-790a-42b7-9430-caa20f004700
clip(z, vclip=40) = abs(z)<vclip ? z : vclip*sign(z)

# ╔═╡ 35408a25-28b4-4279-99e1-041a82148a59
"""
    camera_eye(az_deg, el_deg; zoom=1.0)

Convert azimuth and elevation in degrees to Plotly `camera.eye` coordinates.
Zoom=1 means unit sphere, zoom >1 means zoom in (closer camera).
"""
function camera_eye(az_deg, el_deg; zoom=1.0)
    r = 2.6 / zoom
    az = deg2rad(az_deg-90)
    el = deg2rad(el_deg)

    x = r * cos(el) * cos(az)
    y = r * cos(el) * sin(az)
    z = r * sin(el)

    return Dict(:x => x, :y => y, :z => z)
end

# ╔═╡ 92c6239c-2017-4f7c-8543-cb0553a3235f
extra_kwargs = Dict(
			:series => Dict(
				:contours => Dict(
					:x=>Dict(:show => true, :color=>"black",
							 :start=>xv[1], :end=>xv[end], :size=>0.03),
					:y=>Dict(:show => true, :color=>"black",
							 :start=>yv_dn[1], :end=>yv_up[end], :size=>0.03),
			)),
			:plot => Dict(
				:annotations => [
					Dict(
						:showarrow=>:false,
						:x => 0.75,
					    :y => 0.95,
	      				:text => "Re s",
						:font => Dict(
							:family => "Computer Modern",
    						:color => "black",
	          				:size => 25),
					),
					Dict(
						:showarrow=>:false,
						:x => 0.2,
					    :y => 0.88,
	      				:text => "Im s",
						:font => Dict(
							:family => "Computer Modern",
    						:color => "black",
	          				:size => 25),
					)],
				:scene => Dict(
					:camera => Dict(
						:eye => camera_eye(30, 35; zoom=2.0),
						:center => Dict(:x=>-0.05,:y=>0.0,:z=>-0.3)),
					:aspectratio => Dict(:x=>1.1,:y=>1,:z=>0.5))));
		

# ╔═╡ 27f9523b-b15b-4200-8829-c6060f8ab2df
let
	plot(colorbar=false, size=(700,500), axis=nothing)
	zv = imag.(sheetI.(bw |> Ref, xv' .+ 1im .* yv_up))
	surface!(xv, yv_up, zv; c=sheetI_grad,
		frame=:none, extra_kwargs)
	# 
	zv = imag.(sheetII.(bw |> Ref, xv' .+ 1im .* yv_dn))
	surface!(xv, yv_dn, zv; c=sheetI_grad)
	# 
	path3d!(map(range(bw.mth^2, xv[end], 30)) do x
		z = imag(sheetI(bw, x+iϵ))
		δ=0.1
		(x, -0.003, z+0.7)
	end, l=(13,unitarity_color))
	scatter3d!([branch_point], m=(6,unitarity_color))
	plot!(zlim=(-30,30), clim=(-30,30))
	# 
	savefig(joinpath("sheet-I.png"))
	savefig(joinpath("sheet-I.pdf"))
	plot!()
end

# ╔═╡ 09ed5cb0-ad38-4e68-91f9-3ef4015b8bb2
let
	plot(colorbar=false, size=(700,500), axis=nothing)
	zv = imag.(sheetI.(bw |> Ref, xv' .+ 1im .* yv_up))
	zv = map(clip, zv)
	surface!(xv, yv_up, zv; c=sheetI_grad,
		frame=:none, extra_kwargs)
	# 
	zv = imag.(sheetI.(bw |> Ref, xv' .+ 1im .* yv_dn))
	zv = map(clip, zv)
	surface!(xv, yv_dn, zv; c=sheetII_grad)
	# 
	path3d!(map(range(bw.mth^2, xv[1], 30)) do x
		z = imag(sheetI(bw, x+iϵ))
		δ=0.1
		(x, 0, z+0.7)
	end, l=(13,:black))
	scatter3d!([branch_point], m=(6,unitarity_color))
	plot!(zlim=(-30,30), clim=(-30,30))
	# 
	savefig(joinpath("sheet-I_II.png"))
	savefig(joinpath("sheet-I_II.pdf"))
	plot!()
end

# ╔═╡ 4a4a1e27-42f3-4d24-9736-c7224fa056c3
let
	plot(colorbar=false, size=(700,500), axis=nothing)
	zv = imag.(sheetII.(bw |> Ref, xv' .+ 1im .* yv_up))
	zv = map(clip, zv)
	surface!(xv, yv_up, zv; c=sheetI_grad,
		frame=:none, extra_kwargs)
	# 
	zv = imag.(sheetI.(bw |> Ref, xv' .+ 1im .* yv_dn))
	zv = map(clip, zv)
	surface!(xv, yv_dn, zv; c=sheetII_grad)
	# 
	path3d!(map(range(bw.mth^2, xv[end], 30)) do x
		z = imag(sheetI(bw, x+iϵ))
		δ=0.1
		(x, -0.003, z+0.7)
	end, l=(13,unitarity_color))
	scatter3d!([branch_point], m=(6,unitarity_color))
	plot!(zlim=(-30,30), clim=(-30,30))
	# 
	savefig(joinpath("sheet-II.png"))
	savefig(joinpath("sheet-II.pdf"))
	plot!()
end

# ╔═╡ 8cbd4580-cd18-4bb6-be93-9e0f43542591
md"""
## Test ground
"""

# ╔═╡ 2017bece-338a-48d6-bf14-4909739d2eb5
surface(range(-1,1,100), range(-1,1,100), (x,y)->abs2(x+1im*y),
		colorbar=false,
		c=cgrad(:roma),
		frame=:none,
		extra_kwargs = Dict(
			:series => Dict(
				# :hidesurface => true,
				:contours => Dict(
					:z=>Dict(:show => true, :color=>"black",
							 :width => 10,
							 :start=>0.01, :end=>3.2, :size=>0.1)
			)),
			:plot => Dict(
				:yaxis => Dict(
					:color => "orange",
			        :showline => true,     # show main axis line
			        :showgrid => true,     # show grid lines
			        :zeroline => true,    # show x=0 line
					:visible => true,
			        :showbackground => true,
    			),
				:scene => Dict(
					:camera => Dict(
						:eye => camera_eye(60, 30; zoom=1.5),
						:center => Dict(:x=>0,:y=>0,:z=>-0.3)),
					:aspectratio => Dict(:x=>1,:y=>1,:z=>1))))
)

# ╔═╡ Cell order:
# ╟─5d1282ae-70f1-4fcd-b0cc-f97afbb2fc2f
# ╠═d4a12a4c-480b-4323-be97-49ba4d40d817
# ╠═a3768c8d-3a16-4673-a8db-a817183c83f5
# ╠═a9bdb1b6-9a0a-4102-a1e0-d959b7a5f398
# ╠═51440b8d-7a1e-4fd8-9b8e-1fdd0e13ff24
# ╠═8a3ff5f9-dcea-4c7d-9099-c0b5c19fdc53
# ╠═b818c3d0-30ed-4306-a654-940cca662a9f
# ╠═dc07efbe-18cb-4c76-b94d-b2c5465a0e63
# ╠═cdcfbffc-248d-4529-b5c8-ab8cbb1201f2
# ╠═9b80d058-dd66-4084-91dd-2f1305fe010e
# ╠═ad238351-b312-488e-83a8-fff1cd246d83
# ╠═0a4d5bb9-bed9-498b-aeb1-7aa5d38e4729
# ╠═679c5c98-d055-4817-a548-75fb4d8eb158
# ╟─453279b4-8855-40fc-b5a9-61510ce73f9c
# ╠═e7f32f79-e857-4b12-96da-e9b308191ca0
# ╠═407291e4-23df-4c9b-b64b-c9e4a9af9720
# ╟─bc542219-f869-4379-9734-f6ebadcf6b57
# ╠═8aefc2d3-7fa3-413e-8020-0539e1cc663b
# ╠═6d14b882-9469-4ed8-84b0-25d72df3c02d
# ╠═92c6239c-2017-4f7c-8543-cb0553a3235f
# ╠═27f9523b-b15b-4200-8829-c6060f8ab2df
# ╠═3d659a6b-790a-42b7-9430-caa20f004700
# ╠═09ed5cb0-ad38-4e68-91f9-3ef4015b8bb2
# ╠═4a4a1e27-42f3-4d24-9736-c7224fa056c3
# ╟─35408a25-28b4-4279-99e1-041a82148a59
# ╟─8cbd4580-cd18-4bb6-be93-9e0f43542591
# ╠═2017bece-338a-48d6-bf14-4909739d2eb5
