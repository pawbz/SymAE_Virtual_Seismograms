### A Pluto.jl notebook ###
# v0.19.11

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

# ╔═╡ 812796f8-e38b-40df-a992-d3a0860d076c
using DrWatson

# ╔═╡ 84d77632-f263-45b7-95be-f01f7614de24
@quickactivate "SymAE_Virtual_Seismograms"

# ╔═╡ 100ea3fc-4029-4b41-a7a4-f88d77cc5c55
using PlutoUI

# ╔═╡ a9fb1094-2506-11ed-1e18-f903543c0ec7
md"""
## Learning earthquake sources using symmetric autoencoders
"""

# ╔═╡ 4c158cb0-35da-489a-894d-a4bfdce45698
datadir("wrfweg")

# ╔═╡ 090d84cd-492c-4a1e-a2a0-cce52e0df6e4
@bind x Slider(range(0, stop=100, length=10))

# ╔═╡ b0e56346-9a1c-4e3a-97f4-e949a9113487
x

# ╔═╡ Cell order:
# ╠═a9fb1094-2506-11ed-1e18-f903543c0ec7
# ╠═812796f8-e38b-40df-a992-d3a0860d076c
# ╠═84d77632-f263-45b7-95be-f01f7614de24
# ╠═4c158cb0-35da-489a-894d-a4bfdce45698
# ╠═100ea3fc-4029-4b41-a7a4-f88d77cc5c55
# ╠═090d84cd-492c-4a1e-a2a0-cce52e0df6e4
# ╠═b0e56346-9a1c-4e3a-97f4-e949a9113487
