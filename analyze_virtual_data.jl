### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 71ffaefe-2441-11ed-24a9-5d47d1677675
begin
    using PlutoPlotly
    using FFTW
    using PlutoUI
    using CSV
	using PlutoTeachingTools
    using DataFrames
    using DSP
    using Statistics
    using JLD2
    using HDF5
end

# ╔═╡ 0ae12f40-67d7-4507-bedd-d9083edb9127
ChooseDisplayMode()

# ╔═╡ 6e3e2f23-6223-4b5d-ae1b-ba0d924ecda9
TableOfContents()

# ╔═╡ 01a81689-bcd2-400e-9495-fc6814517ff2
md"""
# Four earthquake-receiver pairs
"""

# ╔═╡ c3b13854-8379-4a0c-b51d-ec9bf9f91a43
md"""
$(
@bind plot_seis MultiCheckBox(["Observed", "Virtual", "Source Signatures (Dereverberated Virtual Seismograms)"]
, default=["Observed", "Virtual", "Source Signatures (Dereverberated Virtual Seismograms)"]))
"""

# ╔═╡ 57f27002-f865-4634-b474-d61e8dc3642a
md"""
- The aim is to study earthquake sources using seismograms plotted in red. Chile is a shallow source; so its seismogram is _heavily_ affected by scattering. The scattering is associated with the free surface, e.g., pP, PP, sP arrivals. On the other hand, Bonin is a deep source with less scattering immediately after P. We wish scattering for Chile was as low as Bonin. Alas!
- Symmetric autoencoders have learned a useful representation of seismograms. We can redatum scattering effects from Bonin to Chile to generate virtual seismograms! Virtual Chile seismogram has Bonin-like scattering.
- We can even ask the network to remove scattering altogether to dereverberate seismograms! In this animation, virtual seismograms are generated for different choices of source and scattering effects.
"""

# ╔═╡ f5a4ed58-708b-4bf7-a8d7-0b93ca9594cb
ThreeColumn(md"[^reference]: Redatuming physical systems using symmetric autoencoders, Phys. Rev. Research, 2022", md"[^credits]: Geophysical Inversion Group, Indian Institute of Science (IISc)", md"[^animation]: Pluto.jl (plutojl.org)")

# ╔═╡ 9933486a-51dd-4852-97f7-7a3fe88ae013
md"""# Virtual seismograms!
#### Dotted Seismograms => Real; Solid Seismograms => Virtual
"""

# ╔═╡ 224bf128-c221-4e13-a498-c91a0f2a2549
md"# Summary"

# ╔═╡ b6e8fdd0-4105-4942-9609-2a863ef0895a
PlutoUI.LocalResource("eq_symae_summary.png")

# ╔═╡ 6ce3fbb4-2e7e-41b6-951a-1b4e3e40fe59
md"""
# Appendix
"""

# ╔═╡ 77ea52ba-b35b-49d5-b10d-4605505aa804
fnames=["data7", "data59", "data94", "data176","datasynt"]

# ╔═╡ a476f8c5-86f3-4095-ad91-f05e5fd56bcb
data = broadcast(fnames) do d
	DataFrame(CSV.File(string("data/vdata_v2/",d,".csv")))
end;

# ╔═╡ c2008296-86c0-4270-b455-c54b834a9c15
path_names = ["〽️From Chile<br>To 1A.CORRE", "〽️From Okhotsk<br>To CN.CLRN", "〽️From Fiji<br>To AU.KELC", "〽️From Bonin<br>To IC.HIA", "No 〽️<br>Dereverberation"]

# ╔═╡ 7d29f6a0-a641-4514-966d-5fec9cb8f231
ThreeColumn(md"""##### Scattering Paths:""", md"""
$(@bind selected_path1 Select(path_names, default=rand(path_names)))""", 
	md"""$(@bind selected_path2 Select(path_names, default=rand(path_names)))
""")

# ╔═╡ 3e06fc95-5a5b-456a-8c12-eba22a57133e
eq_full_names=Dict(["'chl1_8'"=>"💥Chile<br>2010<br>Mw 8.8", " 'kam_4'"=>"💥Kamchatka<br> 2018<br>Mw 6.1", " 'okt5_24'"=>"💥Okhotsk<br>2013<br>Mw 6.7"," 'hnd1_12'"=>"💥Hindukush<br>2015<br>Mw 7.5"," 'fij1_34'"=>"💥Fiji<br>2018<br>Mw 8.2"," 'nzd_0'"=>"💥New Zealand<br>2016<br>Mw 7.8", " 'bon1_70'"=>"💥Bonin<br>2015<br>Mw 7.8"," 'rat_4'"=>"💥Rat Island<br>2014<br>Mw 7.9", " 'pb2_23'"=>"💥Peru<br>2015<br>Mw 7.6"," 'dnl_16'"=>"💥Denali<br>2002<br>Mw 7.9"])

# ╔═╡ c050b9aa-11b8-4b9d-a181-1dac83496b91
ThreeColumn(md"""##### Earthquakes:""",md"""
$(@bind selected_eq1 Select(collect(eq_full_names), default=rand(eq_full_names)))""", 
	md"""##### 
$(@bind selected_eq2 Select(collect(eq_full_names), default=rand(eq_full_names)) )
""")

# ╔═╡ bd626e67-9b83-4b15-a30f-e1360404b3cc
TwoColumn(md"""
#### Select paths
$(@bind selected_paths
MultiCheckBox(path_names, default=["From Chile To 1A.CORRE", "Synthetic PREM P Arrival"], select_all=true, orientation=:column)
)""",  md"""
#### Select quakes
$(@bind selected_eq_names
MultiCheckBox(collect(eq_full_names), default=["Chl1_8", "Okt5_21"], select_all=true)
)""")

# ╔═╡ bdf059d4-da94-4799-ac2d-cafe847f854a
eqdata = DataFrame(CSV.File("data/events_list_12_oct.csv"));

# ╔═╡ ffbbebe4-85bb-4637-a247-64fcca6af185
path_lat_long = [Dict("Lats"=>[-36.122, -67.5828], "Longs"=>[ -72.898, 144.275]), Dict("Lats"=>[-18.1125, -35.9791], "Longs"=>[-178.153, 136.9078]), Dict("Lats"=>[53.1995, -4.6737], "Longs"=>[152.7864, 55.4792]), Dict("Lats"=>[27.8386, 49.2704], "Longs"=>[140.4931, 119.7414])];

# ╔═╡ 1be07d95-c5f2-46fa-8d8a-f87571cb7b24
tgrid = range(-200, stop=200, length=801);

# ╔═╡ e7c83f23-9596-422a-91fc-652be0cfd73f
freqgrid = collect(rfftfreq(length(tgrid), inv(step(tgrid))));

# ╔═╡ e708d208-1f32-4c91-91e1-779af992b443
function normalize(s)
    sstd = std(s)
    return s ./ sstd
end

# ╔═╡ 7d84a028-743a-4782-836f-db813b09a2a3
begin
	
fig = Plot(Layout(title=#"Real/Virtual Seismograms"
"<span style='color: black'>Symmetric Autoencoders<br>"*"<span style='color: black;font-size: 20px;'>Redatum Earthquake (💥) and Wave Scattering (〽️) Effects<br>"*"<span style='color: black'> In "*"<span style='color: red;font-size: 20px;'>Observed" *"<span style='color: black'> Seismograms"*"<span style='color: black'> To Generate "* "<span style='color: blue'>Virtual</span>" * "<span style='color: black'> Seismograms", 
	margin=attr(t=250, b=75, l=75, r=75),
		showlegend=false,
		legend=attr(
		orientation="v"
),
	row_title=attr(font=attr(color="red")),
		xaxis_range=[-100, 200],
		yaxis_range=[-4, 4],
        width=690, height=500,
        # xaxis_title2=,
        template="none",
        yaxis_showgrid=false,
	hovermode=false,
	yaxis=attr(
		automargin=true,
            tickfont=attr(size=10),),
        xaxis=attr(
			# title="e",
			automargin=true,
            tickfont=attr(size=10),
			nticks=5,
            gridwidth=1, gridcolor="Gray"
            # scaleanchor = "x",
            # scaleratio = 1,
        ),
        # yaxis=attr(showticklabels=false,automargin=true),
	font=attr(size=20), Subplots(column_titles = "<span style='font-weight:bold;font-size: 14px;'>".*[selected_path1, selected_path2, "Remove 〽️<br>(Dereverberation)"], row_titles = "<span style='font-weight:bold;font-size: 14px;'>" .* [eq_full_names[selected_eq1], eq_full_names[selected_eq2]], shared_xaxes=true, shared_yaxes=true, 
		x_title="<span style='font-size: 12px;'>"*"Time (s) Relative to P-wave Arrival",
		y_title="<span style='font-size: 12px;'>"*"Normalized Displacement Wavefield",
		horizontal_spacing=0.05, rows=2, cols=3, )))
	# ith eq and jth path
	# plots only if plot_code is in plot_seis
	function add_seismogram(eq, path, i, j, color="black", plot_flag=true)
		d=data[findall(x->x==path, path_names)[1]][!, eq]
		virtual_flag = iszero(d[1])
		virtual_flag && ("Virtual" ∉ plot_seis) && (return nothing)
		!virtual_flag && ("Observed" ∉ plot_seis) && (return nothing)
		if(plot_flag)
		add_trace!(fig, scatter(x=tgrid, y=normalize(d[2:end]),
	   legendgroup="wfr",
    legendgrouptitle_text=eq_full_names[eq],
	   name=path, mode="lines", line=virtual_flag ? attr(color=color) : attr(color="red")), row=i, col=j)
		end
	end
	add_seismogram(selected_eq1, selected_path1, 1, 1, "blue", )
	add_seismogram(selected_eq1, selected_path2, 1, 2, "blue", )
	add_seismogram(selected_eq2, selected_path1, 2, 1, "blue", )
	add_seismogram(selected_eq2, selected_path2, 2, 2, "blue", )
	# if("Virtual (Dereverberated)" in plot_seis)
	add_seismogram(selected_eq1, "No 〽️<br>Dereverberation", 1, 3, "blue", "Source Signatures (Dereverberated Virtual Seismograms)" in plot_seis)
	add_seismogram(selected_eq2, "No 〽️<br>Dereverberation", 2, 3, "blue", "Source Signatures (Dereverberated Virtual Seismograms)" in plot_seis)
	# end
	# add_seismogram(1,2)
	# add_seismogram(2,1)
	# add_seismogram(2,2)
#     add_trace!(fig, scatter(
#             x=randn(5),
#             y=randn(5)), row=1, col=2)
# add_trace!(fig, scatter(
#             x=randn(5),
#             y=randn(5)), row=2, col=1)
    PlutoPlotly.plot(fig)

end

# ╔═╡ 19053683-3a0e-4b80-911d-06a241636e22
data[1][!, selected_eq_names[1]]

# ╔═╡ 9f893b60-98c5-4559-9de8-ead952c25092
function plot_seismograms()
   scatter_plots = vec([(d=data[findall(x->x==pathn, path_names)[1]][!, eq]; scatter(x=tgrid, y=normalize(d[2:end]),
	   legendgroup=eq,
    legendgrouptitle_text=eq_full_names[eq],
	   name=pathn, mode="lines", line=(iszero(d[1])) ? nothing : attr(dash="dot", width=3))) for eq in selected_eq_names, pathn in selected_paths])
     
    p = plot(scatter_plots, Layout(
		title="Real/Virtual Seismograms",
		showlegend=true,
		legend=attr(
		orientation="v"
),
		xaxis_range=[-100, 200],
        width=700, height=400,
        xaxis_title="Time Relative to PREM P (s)",
        template="none",
        yaxis_showgrid=false,
        xaxis=attr(
			automargin=true,
            tickfont=attr(size=15), nticks=20,
            gridwidth=1, gridcolor="Gray"
            # scaleanchor = "x",
            # scaleratio = 1,
        ),
        yaxis=attr(showticklabels=false,automargin=true),
    ))
end

# ╔═╡ a4cf325a-41da-4462-b231-4931d2a8d453
plot_seismograms()

# ╔═╡ 02ed416b-3f1a-45e0-82a7-84bcd52461ab
function plot_scattergeo()
	
	p=plot(
    [scattergeo(
        lat=pll["Lats"],
        lon=pll["Longs"],
        mode="markers+lines",
        hoverinfo="text",
        showlegend=true,
		name=path_names[i],
		marker=attr(symbol=["star", "triangle-up"], size=10),
        projection_type="orthographic",
        landcolor="white",
        oceancolor="MidnightBlue",
        showocean=true,
        lakecolor="LightBlue"
	) for (i,pll) in enumerate(path_lat_long)],Layout(title="Earthquake-Receiver Paths"),
)
end

# ╔═╡ ed2bf741-2064-4341-8885-b9343e6b4bf4
plot_scattergeo()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.9"
DSP = "~0.7.8"
DataFrames = "~1.4.4"
FFTW = "~1.5.0"
HDF5 = "~0.16.13"
JLD2 = "~0.4.30"
PlutoPlotly = "~0.3.6"
PlutoTeachingTools = "~0.2.5"
PlutoUI = "~0.7.49"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "90ff9dd0728f507fa00affa7cfb8f92e65b92e4e"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "844b061c104c408b24537482469400af6075aae4"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.5"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "0e5c14c3bb8a61b3d53b2c0620570c332c8d0663"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.2.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "da8b06f89fce9996443010ef92572b193f8dca1f"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.8"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d4f69885afa5e6149d0cab3818491565cf41446d"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "b5df7c3cab3a00c33c2e09c6bd23982a75e2fbb2"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.16.13"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "4cc2bb72df6ff40b055295fdef6d92955f9dede8"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.2+2"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "c3244ef42b7d4508c638339df1bdbf4353e144db"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.30"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "b289a36229c94e326282f36b3e24416a08dc7bd9"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.21"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "680e733c3a0a9cea9e935c8c2184aea6a63fa0b5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.21"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "18f84637e00b72ba6769034a4b50d79ee40c84a9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.5"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Colors", "Dates", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "PlotlyBase", "PlutoUI", "Reexport"]
git-tree-sha1 = "dec81dcd52748ffc59ce3582e709414ff78d947f"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.3.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "ea3e4ac2e49e3438815f8946fa7673b658e35bdb"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "a14a99e430e42a105c898fcc7f212334bc7be887"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.4"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "90cb983381a9dc7d3dff5fb2d1ee52cd59877412"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "c02bd3c9c3fc8463d3591a62a378f90d2d8ab0f3"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.17"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═0ae12f40-67d7-4507-bedd-d9083edb9127
# ╠═6e3e2f23-6223-4b5d-ae1b-ba0d924ecda9
# ╟─01a81689-bcd2-400e-9495-fc6814517ff2
# ╠═ed2bf741-2064-4341-8885-b9343e6b4bf4
# ╟─7d84a028-743a-4782-836f-db813b09a2a3
# ╟─c050b9aa-11b8-4b9d-a181-1dac83496b91
# ╟─7d29f6a0-a641-4514-966d-5fec9cb8f231
# ╟─c3b13854-8379-4a0c-b51d-ec9bf9f91a43
# ╟─57f27002-f865-4634-b474-d61e8dc3642a
# ╟─f5a4ed58-708b-4bf7-a8d7-0b93ca9594cb
# ╟─9933486a-51dd-4852-97f7-7a3fe88ae013
# ╠═a4cf325a-41da-4462-b231-4931d2a8d453
# ╠═bd626e67-9b83-4b15-a30f-e1360404b3cc
# ╟─224bf128-c221-4e13-a498-c91a0f2a2549
# ╟─b6e8fdd0-4105-4942-9609-2a863ef0895a
# ╟─6ce3fbb4-2e7e-41b6-951a-1b4e3e40fe59
# ╠═a476f8c5-86f3-4095-ad91-f05e5fd56bcb
# ╠═77ea52ba-b35b-49d5-b10d-4605505aa804
# ╠═c2008296-86c0-4270-b455-c54b834a9c15
# ╠═3e06fc95-5a5b-456a-8c12-eba22a57133e
# ╠═bdf059d4-da94-4799-ac2d-cafe847f854a
# ╠═ffbbebe4-85bb-4637-a247-64fcca6af185
# ╠═1be07d95-c5f2-46fa-8d8a-f87571cb7b24
# ╠═e7c83f23-9596-422a-91fc-652be0cfd73f
# ╠═e708d208-1f32-4c91-91e1-779af992b443
# ╠═19053683-3a0e-4b80-911d-06a241636e22
# ╠═9f893b60-98c5-4559-9de8-ead952c25092
# ╠═02ed416b-3f1a-45e0-82a7-84bcd52461ab
# ╠═71ffaefe-2441-11ed-24a9-5d47d1677675
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
