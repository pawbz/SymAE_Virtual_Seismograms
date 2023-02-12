### A Pluto.jl notebook ###
# v0.19.19

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
    using PlotlyKaleido
    using FFTW
    using PlutoUI
    using CSV
    using StatsBase
    using DataFrames
    using DSP
    using Statistics
    using JLD2
    using Healpix
    using HDF5
    using FileIO
    using SignalAnalysis
    using NumericalIntegration
    using PlutoTeachingTools
    using ContinuousWavelets, Wavelets
    using LinearAlgebra
end

# ╔═╡ 9e65a55e-d65b-4a70-bf5f-2bc2d2cf1481
ChooseDisplayMode()

# ╔═╡ adba76b2-d470-4824-a289-21cf8d2b5813
TableOfContents()

# ╔═╡ 2c4afc82-5889-452a-b051-9b832ee23e47
md"# Learning earthquake source time functions via redatuming"

# ╔═╡ eb6a8998-5de8-4023-9669-59807042f4cb
md"""
#### Select quakes
"""

# ╔═╡ c57a0b22-6102-4a24-93af-379ed11f77be
md"""# CWT
"""

# ╔═╡ d66c3edd-5a9d-4acf-ac8d-de51483f24df
md"# Envelope"

# ╔═╡ 5a187589-73e9-4e9f-b814-065a0db430b5


# ╔═╡ 25c4170d-c52d-429c-8e78-442027e35781
md"""
# Energy Rise
Would you like to take mean over all bins?

$(@bind avg_over_bins Select([true => "Yes", false => "No"]))
"""

# ╔═╡ 8c3230fc-85eb-4c1a-8c26-9d88b82d1b95
md"# Appendix"

# ╔═╡ bdf059d4-da94-4799-ac2d-cafe847f854a
eqdata = DataFrame(CSV.File("data/events_list_12_oct.csv"))

# ╔═╡ d57677dc-19b6-4621-96ee-123cceb1496e
@bind eq_names confirm(MultiCheckBox(eqdata[!, "Code"], default=["okt1"], select_all=true,))

# ╔═╡ 332688ba-28d2-4f3d-b445-6d21fe22d8b5
stf = h5open("data/stf_virtual_syn_pb10s.hdf5");

# ╔═╡ c7af0751-9004-40e5-b0c8-4ae7c55117e3
md"""
Select spectrogram bin for $(eq_names[1]): $(@bind spectrogram_bin Slider(keys(stf[eq_names[1]])[sortperm(parse.(Int, keys(stf[eq_names[1]])))], show_value=true))
"""


# ╔═╡ 9d5e12da-2f4b-4896-b8ac-ba9e9577f3b6
@bind bins confirm(PlutoUI.combine() do Child
    components = ["x", "y", "z"]
    s = [md"""
      earthquake: $eq $(
      	Child(MultiCheckBox(keys(stf[eq]), select_all=true, default=rand(keys(stf[eq]),1)))
      	)""" for eq in eq_names]

    md"""
    #### Select Bins
    $(s)
    """
end)

# ╔═╡ 5b65dc66-7e18-4f0e-bd4b-8b94100b332a
@bind flow Slider(range(0, stop=1, length=100), show_value=true, default=0.24)

# ╔═╡ 1be07d95-c5f2-46fa-8d8a-f87571cb7b24
tgrid = range(-200, stop=200, length=801)

# ╔═╡ e7c83f23-9596-422a-91fc-652be0cfd73f
freqgrid = collect(rfftfreq(length(tgrid), inv(step(tgrid))))

# ╔═╡ 9b0ea5fc-87b8-4dca-af4a-aed203cc4e9c
function window(s)
    return DSP.Windows.tukey(length(s), 0.2) .* s
end

# ╔═╡ e708d208-1f32-4c91-91e1-779af992b443
function normalize(s)
    sstd = std(s)
    return s ./ sstd
end

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


# ╔═╡ 496ad4f7-54fa-4183-9151-6c7546451990
stf_energy_cumsum = map(stf) do s
    mean(map(s) do sb
        ss = collect(sb["virtual_seis"])
        return cumsum(abs2.(ss ./ norm(ss, 2)))
    end)
end

# ╔═╡ e7f6daa9-9d66-493f-8d6f-d7d8cccc5dd2
@warn "Note that there is a difference between python and julia bin healpix indexing."

# ╔═╡ d70a130a-5ed0-4c8d-bd5f-34a5ef7dca32
all_bins = [keys(stf[eq]) for eq in eq_names];

# ╔═╡ 3db67b6e-e1be-4309-8ccf-7f024cf7a5d3
# collect sorted bins of all earthquakes
all_sorted_bins = [all_bins[ieq][sortperm(parse.(Int, all_bins[ieq]))] for ieq in 1:length(eq_names)];

# ╔═╡ 2cc04b6e-8e9b-4990-a6bd-3575b81d5f9f
# collect sorted bins of all earthquakes
sorted_bins = [bins[ieq][sortperm(parse.(Int, bins[ieq]))] for ieq in 1:length(eq_names)];

# ╔═╡ 0b3a2252-6043-46e3-9f1b-6f071b756256
sorted_bin_indices = [broadcast(x -> parse(Int, x), sorted_bins[ieq]) for ieq in 1:length(eq_names)];

# ╔═╡ 943e472b-a1d4-47fb-82ef-b4bc14aa5136
all_sorted_bin_indices = [broadcast(x -> parse(Int, x), all_sorted_bins[ieq]) for ieq in 1:length(eq_names)];

# ╔═╡ f56a56ff-f415-4843-a08a-18ac9ecee1fc
function plot_energyrise()
    if (avg_over_bins)

        scatter_plots = map(enumerate(eq_names)) do (ieq, eq_name)
            # take mean of stf across all bins
            splot = mean(map(enumerate(keys(stf[eq_name]))) do (ibin, bin)
                s = window(collect(stf[eq_name][string(bin)]["virtual_seis"]))
                return cumsum(abs2.(s ./ norm(s, 2)))
            end)
            return scatter(x=tgrid, y=splot, name=eq_names[ieq])
        end
    else
        scatter_plots = map(enumerate(eq_names)) do (ieq, eq_name)
            map(enumerate(sorted_bin_indices[ieq])) do (ibin, bin)
                s = window(collect(stf[eq_name][string(bin)]["virtual_seis"]))
                splot = cumsum(abs2.(s ./ norm(s, 2)))
                return scatter(x=tgrid, y=splot, name=bin,
                    legendgroup=eq_names[ieq],
                    legendgrouptitle_text=eq_names[ieq])
            end
        end
    end

    p = plot(vcat(scatter_plots...), Layout(
        legend=attr(
            x=1,
            y=1.02,
            yanchor="bottom",
            xanchor="right",
            orientation="h"
        ),
        width=650, height=500,
        xaxis_title="Time Relative to PREM P (s)",
        template="none",
        yaxis_showgrid=false,
        xaxis=attr(
            tickfont=attr(size=15), nticks=20,
            gridwidth=1, gridcolor="Gray"
            # scaleanchor = "x",
            # scaleratio = 1,
        ),
        yaxis=attr(showticklabels=false,),
        legend_title="Second Group Title",
    ))
end

# ╔═╡ e290dbf1-c0da-4619-9590-94271ed15aa4
penergy = plot_energyrise()

# ╔═╡ 9f893b60-98c5-4559-9de8-ead952c25092
function plot_stf()
    scatter_plots = map(enumerate(eq_names)) do (ieq, eq_name)

        map(enumerate(sorted_bin_indices[ieq])) do (ibin, bin)
            angles = broadcast(pix2angRing(Resolution(4), bin + 1)) do x
                floor(Int, rad2deg(x))
            end
            s = collect(stf[eq_names[ieq]][string(bin)]["virtual_seis"])
            sout = get_stf(s)

            # add 10 per each stf, up to previous eq
            if (ieq == 1)
                previous_yoffset = 0
            else
                previous_yoffset = sum(map(1:ieq-1) do iieq
                    return (10 * (length(sorted_bin_indices[iieq]) + 2))
                end)
            end
            # add 10 per each stf, for current eq
            current_yoffset = (10 * (length(sorted_bin_indices[ieq][1:ibin]) - 1))

            return scatter(x=tgrid, y=sout .+ current_yoffset .+ previous_yoffset, fill="tonexty", name=bin, text=vcat(fill(nothing, 50), string("    ", angles), fill(nothing, 800)), mode="lines+text", textposition="top",
                legendgroup=eq_names[ieq],
                legendgrouptitle_text=eq_names[ieq])
        end
    end

    p = plot(vcat(scatter_plots...), Layout(
        width=650, height=400 + 20 * sum(length.(bins)),
        xaxis_title="Time Relative to PREM P (s)",
        template="none",
        yaxis_showgrid=false,
        xaxis=attr(
            tickfont=attr(size=15), nticks=20,
            gridwidth=1, gridcolor="Gray"
            # scaleanchor = "x",
            # scaleratio = 1,
        ),
        yaxis=attr(showticklabels=false,),
        legend_title="Second Group Title",
    ))
end

# ╔═╡ 677ab4cf-b461-4974-9f14-fc5fe748341f
pstf = plot_stf()

# ╔═╡ 5633826c-484c-4089-8175-449551ab3c4d
function plot_cwt_stf()

    f = collect(stf[eq_names[1]][string(spectrogram_bin)]["virtual_seis"])

    c = wavelet(paul2, averagingType=NoAve(), β=2)
    # c = wavelet(Dog{1}(), averagingType=NoAve(), β=1);
    cwt_stf = ContinuousWavelets.cwt(f, c)



    fig = PlutoPlotly.Plot(Layout(height=400, width=700, Subplots(shared_xaxes=true, rows=2, cols=1, subplot_titles=["Virtual Seismogram" "CWT"])))
    add_trace!(fig, heatmap(x=tgrid, z=abs.(cwt_stf)', colorscale="seismic"), row=2, col=1)
    add_trace!(fig, scatter(y=f, x=tgrid), row=1, col=1)

    return PlutoPlotly.plot(fig)
end

# ╔═╡ f5ab5b34-eeb6-4978-9420-39e1f5a2eee2
pcwt = plot_cwt_stf()

# ╔═╡ 1a0f5073-a974-4056-93ca-efbd2a9dbee7
savefig(pcwt.Plot, "test.pdf")

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

# ╔═╡ 7f6c6656-d6cd-4c6f-a329-4d966db67ac2
md"""

A source time function is virtual 
Model training: During training, the last bin was prepared using synthetic seismograms. These synthetic seismograms are downloaded from the Data Services Products: syngine (IRIS) web service using the prem earth model. One of the deep earthquakes from Peru-Brazil's synthetic data is downloaded and as all the seismograms are synthetic data they do not contain the directivity effect due to which we can put all seismograms in one bin ignoring the scale separations. This data is also pre-processed to filter and resample it so as to match with the other data. In the synthetic traces, the second arrival of pP is observed to be coming around 125 s. Still, we require a path devoid of any second arrivals so as to remove pP ambiguity so we mute the second arrivals by tapering the window and all these muted seismograms are also put in the same bin of synthetic data. So, the last bin contains all muted and unmuted synthetic data, and the network is trained with it. [put images of synthetic trace]

The source time functions are prepared for each and every bin of all the earthquakes. 
We prepared virtual seismograms using sources from each bin and using the synthetic path. 
Redatuming is done by combining source information from each and every bin of all earthquakes 
and the path information is taken from the synthetic muted trace. 
After fusing the PREM-model ideal path with the source information of the real earthquake
1. These virtual seismograms are then tapered at the edges, using a Tukey window.
2. Low-pass - A low pass filter is applied to the virtual data. ??what is the basis of choosing this lowpass value??
3. Window- virtual seismograms are then tapered at the edges, using a tukey-window.
4. Normalize- The virtual seismograms are normalized. 
5. Integrate - As these seismograms are far-field data = velocity and not displacement seismograms. 
So we integrate these seismograms to convert them to displacement seismograms.
6. Window
7. Envelope- We take Hilbert and then absolute of the final data to get our source time functions.
"""

# ╔═╡ 74f22f41-769d-4eb2-98b0-e3891dda4701
pix2angRing(Resolution(4), 192)

# ╔═╡ 43ef5993-f80b-4d60-bd4c-c2060f481684
pix2angRing(Resolution(4), 192)

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

# ╔═╡ c1305c38-884c-434b-8056-6cdff7433c2e
md"## Wideband Ambiguity"

# ╔═╡ 9fa7b7f1-33d7-44d5-9528-d367e4a5b658
"""
DRAFT: trying out wide-band ambiguity function using multirate filtering algo. from DSP.resample
"""
function get_wbxaf(dobs)
    # bsize = min(1000, size(dobs, 3))
    nt = size(dobs, 1)
    nr = size(dobs, 2)

    a = range(0, stop=0.4, length=100)
    αvec = vcat(reverse(1 .- a), (1 .+ a)[2:end])
    # convert to range, needed to Interpolations later
    αvec = range(αvec[1], stop=αvec[end], length=length(αvec))
    # @show length(αvec)

    # normalize
    dobs = dobs ./ std(dobs, dims=1)

    # @showprogress 
    dtest = dobs[:, 1]
    nt_new = length(DSP.resample(dtest, minimum(αvec)))
    nt_lag = 250
    @info "number of significant samples along the first dimension should be less than $nt_new"
    dxobs = zeros(Float32, 2 * nt_lag - 1, nr, nr, length(αvec))
    # dxobs = zeros(Float32, 1, nr, nr, length(αvec))
    tlagvec = range(start=-nt_lag + 1, stop=nt_lag - 1, length=2 * nt_lag - 1)
    for (i, α) in enumerate(αvec)
        dx = view(dxobs, :, :, :, i)
        dobs_scaled = DSP.resample(dobs, α, dims=1)
        # nt_new = min(size(dobs_scaled, 1), size(dobs, 1))
        d1 = view(dobs, 1:nt_new, :, 1)
        d2 = view(dobs_scaled, 1:nt_new, :, 1)

        StatsBase.crosscor!(dx, d1, d2, -nt_lag+1:nt_lag-1, demean=false)

    end

    # return rmul!(dxobs,nt_new), αvec
    return dxobs, αvec, tlagvec
    return dropdims(dxobs, dims=(1)), αvec
end


# ╔═╡ 924bf468-dbe7-4088-a9a7-2af3dfa0562d
function plot_wbxaf()
    # sort bins of first eq
    bins1 = bins[1][sortperm(parse.(Int, bins[1]))]
    # collect virtual seismograms of all bins
    f = mapreduce(hcat, bins1) do bin
        collect(stf[eq_names[1]][string(bin)]["virtual_seis"])
    end
    return get_wbxaf(f)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ContinuousWavelets = "96eb917e-2868-4417-9cb6-27e7ff17528f"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
Healpix = "9f4e344d-96bc-545a-84a3-ae6b9e1b672b"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NumericalIntegration = "e7bfaba1-d571-5449-8927-abc22e82249b"
PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SignalAnalysis = "df1fea92-c066-49dd-8b36-eace3378ea47"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Wavelets = "29a6e085-ba6d-5f35-a997-948ac2efa89a"

[compat]
CSV = "~0.10.8"
ContinuousWavelets = "~1.1.2"
DSP = "~0.7.8"
DataFrames = "~1.4.4"
FFTW = "~1.5.0"
FileIO = "~1.16.0"
HDF5 = "~0.16.13"
Healpix = "~4.1.2"
JLD2 = "~0.4.29"
NumericalIntegration = "~0.3.3"
PlotlyKaleido = "~1.0.0"
PlutoPlotly = "~0.3.6"
PlutoTeachingTools = "~0.2.5"
PlutoUI = "~0.7.49"
SignalAnalysis = "~0.4.2"
StatsBase = "~0.33.21"
Wavelets = "~0.9.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "57b2b23c55ef9141411247b7f610074126ca24fa"

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
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "8c73e96bd6817c2597cfd5615b91fca5deccf1af"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.8"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "3bf60ba2fae10e10f70d53c070424e40a820dac2"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.1.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

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
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.ContinuousWavelets]]
deps = ["AbstractFFTs", "FFTW", "Interpolations", "LinearAlgebra", "SpecialFunctions", "Wavelets"]
git-tree-sha1 = "c1808d31e107c23331c2068ae7f321f65b70a370"
uuid = "96eb917e-2868-4417-9cb6-27e7ff17528f"
version = "1.1.2"

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

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "a7756d098cbabec6b3ac44f369f74915e8cfd70a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.79"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

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

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "9a0472ec2f5409db243160a8b030f94c380167a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.6"

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

[[deps.Healpix]]
deps = ["CFITSIO", "LazyArtifacts", "Libsharp", "LinearAlgebra", "Pkg", "Printf", "RecipesBase", "StaticArrays", "Test"]
git-tree-sha1 = "ca9690ae0c58fa8627a01049a949867ffb8e7d15"
uuid = "9f4e344d-96bc-545a-84a3-ae6b9e1b672b"
version = "4.1.2"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

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
git-tree-sha1 = "0cf92ec945125946352f3d46c96976ab972bde6f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.3.2"

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
git-tree-sha1 = "ec8a9c9f0ecb1c687e34c1fda2699de4d054672a"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.29"

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
git-tree-sha1 = "3236cad4ac05408090221fc259f678d913176055"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.19"

[[deps.Kaleido_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43032da5832754f58d14a91ffbe86d5f176acda9"
uuid = "f7e6163d-2fa5-5f23-b69c-1db539e41963"
version = "0.2.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

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
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

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

[[deps.MetaArrays]]
deps = ["Requires"]
git-tree-sha1 = "6647f7d45a9153162d6561957405c12088caf537"
uuid = "36b8f3f0-b776-11e8-061f-1f20094e1fc8"
version = "0.2.10"

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

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

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

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "6466e524967496866901a78fca3f2e9ea445a559"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlotlyKaleido]]
deps = ["Base64", "JSON", "Kaleido_jll"]
git-tree-sha1 = "64b125713e6ec1b5fac6ae1f9624b8b408ec9cb8"
uuid = "f2990250-8cf9-495f-b13a-cce12b45703c"
version = "1.0.0"

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
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

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
deps = ["SnoopPrecompile"]
git-tree-sha1 = "18c35ed630d7229c5584b945641a73ca83fb5213"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.2"

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
git-tree-sha1 = "dad726963ecea2d8a81e26286f625aee09a91b7c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.4.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SignalAnalysis]]
deps = ["DSP", "Distributions", "DocStringExtensions", "FFTW", "LinearAlgebra", "MetaArrays", "PaddedViews", "Random", "Requires", "SignalBase", "Statistics", "WAV"]
git-tree-sha1 = "08d6c4d5290cd2cbc54993ac4c84ca15f1de3738"
uuid = "df1fea92-c066-49dd-8b36-eace3378ea47"
version = "0.4.2"

[[deps.SignalBase]]
deps = ["Unitful"]
git-tree-sha1 = "14cb05cba5cc89d15e6098e7bb41dcef2606a10a"
uuid = "00c44e92-20f5-44bc-8f45-a1dcef76ba38"
version = "0.1.2"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

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

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "ffc098086f35909741f71ce21d03dadf0d2bfa76"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.11"

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

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

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
git-tree-sha1 = "e4bdc63f5c6d62e80eb1c0043fcc0360d5950ff7"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.10"

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

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d670a70dd3cdbe1c1186f2f17c9a68a7ec24838c"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.2"

[[deps.WAV]]
deps = ["Base64", "FileIO", "Libdl", "Logging"]
git-tree-sha1 = "7e7e1b4686995aaf4ecaaf52f6cd824fa6bd6aa5"
uuid = "8149f6b0-98f6-5db9-b78f-408fbbb8ef88"
version = "1.2.0"

[[deps.Wavelets]]
deps = ["DSP", "FFTW", "LinearAlgebra", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "58f7491c21ab2b1d69368c7f7e8a6a93cbf8b7bf"
uuid = "29a6e085-ba6d-5f35-a997-948ac2efa89a"
version = "0.9.5"

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
# ╠═9e65a55e-d65b-4a70-bf5f-2bc2d2cf1481
# ╠═adba76b2-d470-4824-a289-21cf8d2b5813
# ╟─2c4afc82-5889-452a-b051-9b832ee23e47
# ╟─eb6a8998-5de8-4023-9669-59807042f4cb
# ╟─d57677dc-19b6-4621-96ee-123cceb1496e
# ╟─c57a0b22-6102-4a24-93af-379ed11f77be
# ╟─c7af0751-9004-40e5-b0c8-4ae7c55117e3
# ╠═f5ab5b34-eeb6-4978-9420-39e1f5a2eee2
# ╟─d66c3edd-5a9d-4acf-ac8d-de51483f24df
# ╟─9d5e12da-2f4b-4896-b8ac-ba9e9577f3b6
# ╠═677ab4cf-b461-4974-9f14-fc5fe748341f
# ╠═5a187589-73e9-4e9f-b814-065a0db430b5
# ╟─25c4170d-c52d-429c-8e78-442027e35781
# ╠═e290dbf1-c0da-4619-9590-94271ed15aa4
# ╠═1a0f5073-a974-4056-93ca-efbd2a9dbee7
# ╟─8c3230fc-85eb-4c1a-8c26-9d88b82d1b95
# ╠═bdf059d4-da94-4799-ac2d-cafe847f854a
# ╠═332688ba-28d2-4f3d-b445-6d21fe22d8b5
# ╠═5b65dc66-7e18-4f0e-bd4b-8b94100b332a
# ╠═1be07d95-c5f2-46fa-8d8a-f87571cb7b24
# ╠═e7c83f23-9596-422a-91fc-652be0cfd73f
# ╠═9b0ea5fc-87b8-4dca-af4a-aed203cc4e9c
# ╠═e708d208-1f32-4c91-91e1-779af992b443
# ╠═b9c71123-0c0a-47d1-9dd3-3cda5fcdd7ff
# ╠═945ed2b3-eec4-4002-9c69-e581407a6d6d
# ╠═f2e08b0e-d1ce-498e-8a0e-78f7c6857aff
# ╠═b42e13af-309e-4938-a03f-fb913bb747c0
# ╠═496ad4f7-54fa-4183-9151-6c7546451990
# ╟─e7f6daa9-9d66-493f-8d6f-d7d8cccc5dd2
# ╠═d70a130a-5ed0-4c8d-bd5f-34a5ef7dca32
# ╠═3db67b6e-e1be-4309-8ccf-7f024cf7a5d3
# ╠═2cc04b6e-8e9b-4990-a6bd-3575b81d5f9f
# ╠═0b3a2252-6043-46e3-9f1b-6f071b756256
# ╠═943e472b-a1d4-47fb-82ef-b4bc14aa5136
# ╠═f56a56ff-f415-4843-a08a-18ac9ecee1fc
# ╠═9f893b60-98c5-4559-9de8-ead952c25092
# ╠═5633826c-484c-4089-8175-449551ab3c4d
# ╠═924bf468-dbe7-4088-a9a7-2af3dfa0562d
# ╠═a6f8585e-b3f7-42f4-b638-5821350c387e
# ╠═7f6c6656-d6cd-4c6f-a329-4d966db67ac2
# ╠═71ffaefe-2441-11ed-24a9-5d47d1677675
# ╠═74f22f41-769d-4eb2-98b0-e3891dda4701
# ╠═43ef5993-f80b-4d60-bd4c-c2060f481684
# ╠═f7f756b6-ffa1-4c76-b82a-dd44134a0596
# ╟─c1305c38-884c-434b-8056-6cdff7433c2e
# ╠═9fa7b7f1-33d7-44d5-9528-d367e4a5b658
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
