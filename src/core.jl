using PyCall
using LinearAlgebra
using Plots
using Statistics
using ProgressMeter
using DSP
using StatsBase
using FFTW


function window!(virt, window)
    @assert size(window) == (1,size(virt[1],2))
    broadcast!(virt, virt) do v
        broadcast!(*, v, v, window)
    end
end

function normalize!(virt)
    broadcast!(virt, virt) do v
        vstd=std(v, dims=2)
        broadcast!(inv, vstd, vstd)
        broadcast!(*, v, v, vstd)
    end
end

function get_mean_spectrum(virt, X, F)
    broadcast(virt) do v
        mul!(X, F, v)
        broadcast!(abs, X, X)
        s=mean(X, dims=1)
        return real.(dropdims(s, dims=1))
    end
end





function padfft(x::AbstractMatrix,n)
    a=fill(zero(eltype(x)), (size(x,1), n-size(x,2)))
    return cat(x, a, dims=2)
end

function padfft(virt,n)
    broadcast(virt) do v
        padfft(v,n)
    end
end

function _unwrap!(phase)
    samples = size(phase, 2)
	unwrap!(phase, dims=(2))
    center = div(samples + 1, 2)
    ndelay = round.(phase[:, center:center] ./ pi)
    phase .= phase .- (pi .* ndelay .* reshape(1:samples, 1, samples) ./ center)
    return phase
end

function get_mean_cepstrum(virt, X, Xabs, Xph, F, nt)
   broadcast(virt) do v
        mul!(X, F, v)
        broadcast!(Xabs, X) do x
            log(abs(x)) + eps(Float32)
        end
        broadcast!(angle, Xph, X) 
	    _unwrap!(Xph)
        broadcast!(complex, X, Xabs, Xph)
        s=mean(X, dims=1)
        broadcast!(exp, s, s)
        s=dropdims(s, dims=1)
        s=irfft(s,size(v,2))
        its=div(length(s)-nt, 2)
        return fftshift(s)[its:its+nt-1]
    end
end

function icepstrum(X, n)
    return real(ifft(exp.(fft(X, (1))), (1)))[1:n, :]
end