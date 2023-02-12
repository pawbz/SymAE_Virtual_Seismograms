### A Pluto.jl notebook ###
# v0.19.14

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
    using DataFrames
    using DSP
    using Statistics
    using JLD2
    using Healpix
    using HDF5
	using StatsBase
    using FileIO
    using NumericalIntegration
end

# ╔═╡ bdf059d4-da94-4799-ac2d-cafe847f854a
eqdata = DataFrame(CSV.File("data/events_list_12_oct.csv"))

# ╔═╡ 332688ba-28d2-4f3d-b445-6d21fe22d8b5
stf = h5open("data/stf_virtual_syn_pb10s.hdf5");

# ╔═╡ ec520726-81b6-451b-a021-5413ef13453a
spec = h5open("data/all_spectras_2.0.hdf5");

# ╔═╡ eb6a8998-5de8-4023-9669-59807042f4cb
md"""
#### Select quakes
$(@bind eq_names
MultiCheckBox(eqdata[!, "Code"], default=["okt1"])
)"""

# ╔═╡ 9d5e12da-2f4b-4896-b8ac-ba9e9577f3b6
@bind bins PlutoUI.combine() do Child
    components = ["x", "y", "z"]
    s = [md"""
      earthquake: $eq $(
      	Child(MultiCheckBox(keys(stf[eq]), select_all=true, default=keys(stf[eq])))
      	)""" for eq in eq_names]

    md"""
    #### Select Bins
    $(s)
    """
end

# ╔═╡ b7b20eba-b956-42e6-9c65-77feed77b218
function get_err(a,b)
err = broadcast(a,b) do aa, bb
	msd(aa,bb)
end
i=sortperm(err)
err = err[i]
return err, eqdata[!, "Code"][i]
end

# ╔═╡ 1be07d95-c5f2-46fa-8d8a-f87571cb7b24
tgrid = range(-200, stop=200, length=801)

# ╔═╡ e7c83f23-9596-422a-91fc-652be0cfd73f
freqgrid = collect(rfftfreq(length(tgrid), inv(step(tgrid))))

# ╔═╡ 5b65dc66-7e18-4f0e-bd4b-8b94100b332a
@bind flow Slider(freqgrid, show_value=true, default=0.24)

# ╔═╡ 9b0ea5fc-87b8-4dca-af4a-aed203cc4e9c
function window(s)
    return DSP.Windows.tukey(length(s), 0.2) .* s
end

# ╔═╡ e708d208-1f32-4c91-91e1-779af992b443
function normalize(s)
    sstd = std(s)
    return s ./ sstd
end

# ╔═╡ 91653674-7447-49fb-a79a-84180e6cb03e
function get_all_spectra(name)
s=broadcast(eqdata[!, "Code"]) do eq
	broadcast(keys(spec[eq])) do bin
		collect(spec[eq][bin][name]) |> normalize
end
end
s = broadcast(s) do ss
	hcat(ss...)
end
	return s
end



# ╔═╡ 0d1a1ba1-1a56-4ddb-8a07-84343220d245
spec_onepath = get_all_spectra("spectra_syn_onepath")

# ╔═╡ 2173d092-c7b8-4b2e-a13e-873d4133ca24
spec_allpath = get_all_spectra("stacked_spectra_realpath")

# ╔═╡ 55ba3811-fd00-4110-9f4f-9b5913142820
spec_err1, eq1 = get_err(spec_onepath,spec_allpath)

# ╔═╡ f60d0fc3-9752-4bf1-ba2d-487455414969
plot(scatter(x=eq1,y=spec_err1), Layout(title="Mean-squared error b/w _onepath_ and _all path_ virtual spectrum"))

# ╔═╡ 62f3fab0-3b45-4fb4-9327-97a224a1d69b
spec_rawpath = get_all_spectra("stacked_spectra_raw")

# ╔═╡ 863aa7cb-9c03-4b66-ba23-5e070529ce59
spec_err2, eq2 = get_err(spec_rawpath,spec_allpath)

# ╔═╡ 1b7e27b8-f992-4d99-af66-732c42682158
plot(scatter(x=eq2,y=spec_err2), Layout(title="Mean-squared error b/w _raw path_ and _all path_ virtual spectrum"))

# ╔═╡ 47479c38-8ad4-4bdd-9064-6c3ddd401d50
spec_err3, eq3 = get_err(spec_rawpath,spec_onepath)

# ╔═╡ 1ca0ec42-ffb0-4880-89ea-267e9f2c2ad2
plot(scatter(x=eq3,y=spec_err3), Layout(title="Mean-squared error b/w _raw path_ and _one path_ virtual spectrum"))

# ╔═╡ b9c71123-0c0a-47d1-9dd3-3cda5fcdd7ff
function integrate(s)
    return cumul_integrate(tgrid, s)
end

# ╔═╡ 945ed2b3-eec4-4002-9c69-e581407a6d6d
function lowpass_filter(s)
    responsetype = Lowpass(flow; fs=2)
    designmethod = Butterworth(4)
    return filtfilt(digitalfilter(responsetype, designmethod), s)
end

# ╔═╡ f2e08b0e-d1ce-498e-8a0e-78f7c6857aff
function envelope(s)
    return abs.(hilbert(s))
end

# ╔═╡ b42e13af-309e-4938-a03f-fb913bb747c0
function get_stf(virtual_seismogram)
	virtual_seismogram |> window |> lowpass_filter |> window |> normalize |> integrate |> window |> envelope
end
	

# ╔═╡ e7f6daa9-9d66-493f-8d6f-d7d8cccc5dd2
@warn "Note that there is a difference between python and julia bin healpix indexing."

# ╔═╡ 8aa9387e-c36e-41b3-a45d-e5baeb2b1578
function plot_spec(ieq, name="spectra_syn_onepath")
    bins1 = bins[ieq][sortperm(parse.(Int, bins[ieq]))]
    bin_indices = broadcast(x -> parse(Int, x), bins1)
    scatter_plots = map(bin_indices) do bin
        angles = broadcast(pix2angRing(Resolution(4), bin + 1)) do x
            floor(Int, rad2deg(x))
        end
        s = collect(spec[eq_names[ieq]][string(bin)][name]) |> normalize
        
        return scatter(x=freqgrid, y=s, name=bin,
            legendgrouptitle_text=eq_names[ieq])
    end
    p = plot(scatter_plots, Layout(
        width=700, height=400,
        xaxis_title="Time Relative to PREM P (s)",
        template="none",
		xaxis_range=[-2, 0],
		yaxis_range=[-2, 1.5],
		yaxis_type="log",
		xaxis_type="log",
        xaxis=attr(
            tickfont=attr(size=15), nticks=10,
            gridwidth=1, gridcolor="Gray"
            # scaleanchor = "x",
            # scaleratio = 1,
        ),
        yaxis=attr(),
        legend_title="Second Group Title",
    ))
end

# ╔═╡ 962bc72e-ad41-4abe-91eb-4055d5e2dec7
plot_spec(1, "spectra_syn_onepath" )

# ╔═╡ ffc6b01d-820d-44b4-b260-5a37be0c0cc1
plot_spec(1, "stacked_spectra_raw")

# ╔═╡ 8f7ea00c-4f8c-464d-b576-d60987ddff62
plot_spec(1, "stacked_spectra_realpath")

# ╔═╡ 5788e16a-3bd4-4b9f-8929-b9e06a30b6f2
bins[1][sortperm(parse.(Int, bins[1]))]

# ╔═╡ 33c409e5-df66-45bf-8eb4-5d3a52f99be9
bins

# ╔═╡ a6f8585e-b3f7-42f4-b638-5821350c387e
begin
    fig = Figure()
    ga = GeoAxis(
        fig[1, 1]; # any cell of the figure's layout
        dest="+proj=wintri", # the CRS in which you want to plot
        coastlines=true # plot coastlines from Natural Earth, as a reference.
    )
    scatter!(ga, -120:15:120, -60:7.5:60; color=-60:7.5:60, strokecolor=(:black, 0.2))
    fig
end

# ╔═╡ 6ce3fbb4-2e7e-41b6-951a-1b4e3e40fe59
md"""
## Appendix
"""

# ╔═╡ 74f22f41-769d-4eb2-98b0-e3891dda4701
pix2angRing(Resolution(4), 192)

# ╔═╡ 43ef5993-f80b-4d60-bd4c-c2060f481684
pix2angRing(Resolution(4), 192)

# ╔═╡ 317001fb-12b0-4d61-b417-3bea0d1eb7a8
function plot_spectra()
    f = Figure()
    tgrid = range(-200, stop=200, length=801)

    Axis(f[1, 1], yscale=log10, xscale=log10)

    xs = rfftfreq(length(tgrid), inv(step(tgrid)))[2:end]
    for eq in eq_names, ii in i
        lines!(xs, spec[eq][ii][2:end], label=string(eq, " bin ", ii), linewidth=3)
    end
    axislegend(; position=:rt)
    limits!(0.01, 1, 0.1, 1000) # x1, x2, y1, y2
    f
end

# ╔═╡ f7f756b6-ffa1-4c76-b82a-dd44134a0596
begin
    # isolate binindices from fnames as Int
    binindices(fnames) =
        broadcast(fnames) do f
            parse(Int, split(splitext(f)[1], "_")[end])
        end
    # get binindices of a specific eq
    binindices(fnames, eqname) = binindices(eqfiles(fnames, eqname))
    # filter jld2 files of an eq
    function eqfiles(fnames, eqname)
        return filter(fnames) do f
            occursin(string("/", eqname, "_"), f) && (splitext(f)[end] == ".jld2")
        end
    end
    # give filename of a particular eq and index
    function eqfiles(fnames, eqname, binindex)
        f = filter(eqfiles(fnames, eqname)) do f
            occursin(string("_", binindex, "."), f)
        end
        return first(f)
    end
    function load(eqname, binindex)
        return vec(JLD2.load(joinpath("..", "data", eqfiles(data_files, eqname, binindex)), "data"))

    end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
Healpix = "9f4e344d-96bc-545a-84a3-ae6b9e1b672b"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
NumericalIntegration = "e7bfaba1-d571-5449-8927-abc22e82249b"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.10.4"
DSP = "~0.7.7"
DataFrames = "~1.3.4"
FFTW = "~1.5.0"
FileIO = "~1.15.0"
HDF5 = "~0.16.12"
Healpix = "~4.1.2"
JLD2 = "~0.4.22"
NumericalIntegration = "~0.3.3"
PlutoPlotly = "~0.3.6"
PlutoUI = "~0.7.39"
StatsBase = "~0.33.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "73ed459cfc2a48f4e2c2d4ae494553ec00cfae7e"

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

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CFITSIO]]
deps = ["CFITSIO_jll"]
git-tree-sha1 = "8425c47db102577eefb93cb37b4480e750116b0d"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.4.1"

[[deps.CFITSIO_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "9c91a9358de42043c3101e3a29e60883345b0b39"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "4.0.0+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "78bee250c6826e1cf805a88b7f1e86025275d208"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.46.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "4ba2a190a9d05a36e8c26182eb1ba06cd12c1051"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.7"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

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
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

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
git-tree-sha1 = "94f5101b96d2d968ace56f7f2db19d0a5f592e28"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.15.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "316daa94fad0b7a008ebd573e002efd6609d85ac"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.19"

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
deps = ["Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires"]
git-tree-sha1 = "19effd6b5af759c8aaeb9c77f89422d3f975ab65"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.16.12"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "4cc2bb72df6ff40b055295fdef6d92955f9dede8"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.2+2"

[[deps.Healpix]]
deps = ["CFITSIO", "LazyArtifacts", "Libsharp", "LinearAlgebra", "Pkg", "Printf", "RecipesBase", "StaticArrays", "Test"]
git-tree-sha1 = "ca9690ae0c58fa8627a01049a949867ffb8e7d15"
uuid = "9f4e344d-96bc-545a-84a3-ae6b9e1b672b"
version = "4.1.2"

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
git-tree-sha1 = "d19f9edd8c34760dca2de2b503f969d8700ed288"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

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
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

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

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[deps.Libsharp]]
deps = ["Libdl", "libsharp2_jll"]
git-tree-sha1 = "e09051a3f95b83091fc9b7a26e6c585c6691b5bc"
uuid = "ac8d63fe-4615-43ae-9860-9cd4a3820532"
version = "0.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "41d162ae9c868218b1f3fe78cba878aa348c2d26"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NumericalIntegration]]
deps = ["Interpolations", "LinearAlgebra", "Logging"]
git-tree-sha1 = "2a4ef5fc235053f9747d59cfdee19bcb8ba1e833"
uuid = "e7bfaba1-d571-5449-8927-abc22e82249b"
version = "0.3.3"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

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
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

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
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Colors", "Dates", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "PlotlyBase", "PlutoUI", "Reexport"]
git-tree-sha1 = "dec81dcd52748ffc59ce3582e709414ff78d947f"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.3.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "3010a6dd6ad4c7384d2f38c58fa8172797d879c1"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.0"

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
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "f86b3a049e5d05227b10e15dbb315c5b90f14988"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.9"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

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
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

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
git-tree-sha1 = "ed5d390c7addb70e90fd1eb783dcb9897922cbfa"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.8"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

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

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libsharp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9a7b8663b021ba776c89d9de7cdb2069cc27c00a"
uuid = "180207a7-b08e-5162-af94-7d62a04fe081"
version = "1.0.2+1"

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
# ╠═bdf059d4-da94-4799-ac2d-cafe847f854a
# ╠═332688ba-28d2-4f3d-b445-6d21fe22d8b5
# ╠═ec520726-81b6-451b-a021-5413ef13453a
# ╠═eb6a8998-5de8-4023-9669-59807042f4cb
# ╠═9d5e12da-2f4b-4896-b8ac-ba9e9577f3b6
# ╠═5b65dc66-7e18-4f0e-bd4b-8b94100b332a
# ╠═91653674-7447-49fb-a79a-84180e6cb03e
# ╠═0d1a1ba1-1a56-4ddb-8a07-84343220d245
# ╠═2173d092-c7b8-4b2e-a13e-873d4133ca24
# ╠═62f3fab0-3b45-4fb4-9327-97a224a1d69b
# ╠═b7b20eba-b956-42e6-9c65-77feed77b218
# ╠═55ba3811-fd00-4110-9f4f-9b5913142820
# ╠═863aa7cb-9c03-4b66-ba23-5e070529ce59
# ╠═47479c38-8ad4-4bdd-9064-6c3ddd401d50
# ╠═f60d0fc3-9752-4bf1-ba2d-487455414969
# ╠═1b7e27b8-f992-4d99-af66-732c42682158
# ╠═1ca0ec42-ffb0-4880-89ea-267e9f2c2ad2
# ╠═962bc72e-ad41-4abe-91eb-4055d5e2dec7
# ╠═ffc6b01d-820d-44b4-b260-5a37be0c0cc1
# ╠═8f7ea00c-4f8c-464d-b576-d60987ddff62
# ╠═1be07d95-c5f2-46fa-8d8a-f87571cb7b24
# ╠═e7c83f23-9596-422a-91fc-652be0cfd73f
# ╠═9b0ea5fc-87b8-4dca-af4a-aed203cc4e9c
# ╠═e708d208-1f32-4c91-91e1-779af992b443
# ╠═b9c71123-0c0a-47d1-9dd3-3cda5fcdd7ff
# ╠═945ed2b3-eec4-4002-9c69-e581407a6d6d
# ╠═f2e08b0e-d1ce-498e-8a0e-78f7c6857aff
# ╠═b42e13af-309e-4938-a03f-fb913bb747c0
# ╟─e7f6daa9-9d66-493f-8d6f-d7d8cccc5dd2
# ╠═8aa9387e-c36e-41b3-a45d-e5baeb2b1578
# ╠═5788e16a-3bd4-4b9f-8929-b9e06a30b6f2
# ╠═33c409e5-df66-45bf-8eb4-5d3a52f99be9
# ╠═a6f8585e-b3f7-42f4-b638-5821350c387e
# ╟─6ce3fbb4-2e7e-41b6-951a-1b4e3e40fe59
# ╠═71ffaefe-2441-11ed-24a9-5d47d1677675
# ╠═74f22f41-769d-4eb2-98b0-e3891dda4701
# ╠═43ef5993-f80b-4d60-bd4c-c2060f481684
# ╠═317001fb-12b0-4d61-b417-3bea0d1eb7a8
# ╠═f7f756b6-ffa1-4c76-b82a-dd44134a0596
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
