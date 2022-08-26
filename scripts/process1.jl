using DrWatson
@quickactivate "SymAE_Virtual_Seismograms"
include(srcdir("core.jl"))

@pyinclude(scriptsdir("prepare_32eq.py"))


nt=801
time_window=reshape(DSP.Windows.tukey(nt, 0.2), 1, :)


# get mean spectra
spec=Dict()
@showprogress for eq_name in py"eq_names"
    py"""
    virt=generate_virt($eq_name)
    """
    virt=py"virt"
    window!(virt, time_window)
    normalize!(virt)
    F = plan_rfft(zero(virt[1]), (2))
    buffer = F*zero(virt[1])
    spec[eq_name]=get_mean_spectrum(virt, buffer, F)
end

safesave(datadir("mean_spectra.jld2"), spec)





# need alpha=1 for cepstrum summation, need to bring poles into unit circle
time_window=reshape(DSP.Windows.tukey(nt, 1), 1, :)
# get mean spectra
stf=Dict()
nfft=4096
@showprogress for eq_name in py"eq_names"
    py"""
    virt=generate_virt($eq_name)
    """
    virt=py"virt"
    window!(virt, time_window)
    normalize!(virt)
    virt=padfft(virt, nfft)
    F = plan_rfft(zero(virt[1]), (2))
    X = F*zero(virt[1])
    Xabs = copy(real.(X))
    Xph = copy(real.(X))
    stf[eq_name]=get_mean_cepstrum(virt, X, Xabs, Xph, F, nt)
end

wsave(datadir("mean_cepstrum_stf.jld2"), stf)
