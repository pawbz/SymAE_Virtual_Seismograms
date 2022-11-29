### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 8d521710-6f90-11ed-2bf4-d9fbc93ef2b0
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

    using PlutoUI, LinearAlgebra, FFTW, Plots
end

# ╔═╡ 85ae8480-7748-45d4-be3a-3199bdf04061
include("./test_module.jl")

# ╔═╡ 37448356-f86f-4389-9ce5-72f7800446e2
plot(randn(10))

# ╔═╡ 72f02db8-6834-4443-959e-540a6d6d2b62
plot(fft(randn(10)))

# ╔═╡ e0f3cc2c-cd53-4f10-b641-948b823a7856
T.test_plot()

# ╔═╡ Cell order:
# ╠═8d521710-6f90-11ed-2bf4-d9fbc93ef2b0
# ╠═85ae8480-7748-45d4-be3a-3199bdf04061
# ╠═37448356-f86f-4389-9ce5-72f7800446e2
# ╠═72f02db8-6834-4443-959e-540a6d6d2b62
# ╠═e0f3cc2c-cd53-4f10-b641-948b823a7856
