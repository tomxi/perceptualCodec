function ThrQuietBM = init_thr_in_quiet(p0_level_dB, Nbins, fs, band_binWidth)
%function init_thr_in_quiet(p0_level_dB, N, fs, band_binWidth)
%N is number of spectral coefs
%fs is sampling frequency
%p0_level_dB is 1/2 lsb noise level of input signal
%
global plt;
%
%Bosi, page 155
f = ([1:Nbins]'/Nbins)*fs/2;
ThrQuiet_dB = 3.64*((f/1000).^-0.8) - 6.5*exp(-0.6*(f/1000 - 3.3).^2) + ...
    (10^-3)*((f/1000).^4);
%clip to not more than 80 dB
ThrQuiet_dB(ThrQuiet_dB > 80) = 80;
%set min level to p0
min_thr_dB = min(ThrQuiet_dB);
ThrQuiet_dB = ThrQuiet_dB - min_thr_dB + p0_level_dB;
if plt
    xs = fs/2*[0:Nbins-1]/Nbins;
    semilogx(xs, ThrQuiet_dB)
    grid
    title('Threshold in Quiet')
    xlabel('Hz')
    ylabel('dB')
    pause
end
%convert to Bark scale
ThrQuietMagHz = 10.^(ThrQuiet_dB/20);
ThrQuietBM = hz2bark(ThrQuietMagHz, band_binWidth);

end
