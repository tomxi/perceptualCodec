function ThrBdB = perceptual_model(arg)
%function thr = perceptual_model(arg)
%
%if arg is scalar, then it is p0 for threshold of hearing
%   then initialize helper variables
%if are is vector, then it is windowed signal block
%   then compute masking threshold

global plt
global fs N SMR
global Zmax bands_per_bark Nbands N2
global band_binWidth spread_ftn ThrQuietBdB

%persistent Zmax bands_per_bark Nbands N2
%persistent band_binWidth spread_ftn ThrQuietBdB

%initialize and return
[r,c] = size(arg);
if (r == 1 && c ==1 )
    p0 = arg;
    N2 = N/2;
    
    bands_per_bark = 2;
    [band2bin, band_binWidth, Zmax] = init_bark(fs, N2, bands_per_bark);

    %e.g. 25 Bark, 2 bands per Bark
    Nbands = Zmax*bands_per_bark;
    spread_ftn = init_spread(Zmax, bands_per_bark);

    %thr in quiet
    ThrQuietBM = thr_in_quiet(p0, N2, fs, band_binWidth);
    ThrQuietBdB = 20*log10(ThrQuietBM + realmin);
    return
end

%otherwise process
if (c == 1)
    wx = arg;
else
    fprintf('Must be mono\n');
    return;
end

%magnitude spectrum
XM = mag_spec(wx);
if (plt)
    loglog(XM);
    grid
    title('Signal Hz Magnitude Spectrum');
    pause
end

%bark spectrum
XBM = hz2bark(XM, band_binWidth);
if (plt)
    semilogy(XBM);
    grid
    title('Signal Bark Magnitude Spectrum');
    pause
end

%spread spectrum
%pad XBM with zeros so spreading function can run to end
XSpread = filter(spread_ftn, 1, [XBM; zeros(Nbands,1)]);
%now compensate for width of spread_ftn
XSpread = XSpread(Nbands+1:end);
if (plt)
    semilogy([.5:.5:25], XSpread);
    grid
    title('Signal Bark Magnitude Spread Spectrum');
    pause
end

%convert to dB
XBdB = 20*log10(XBM + realmin);
XBSdB = 20*log10(XSpread + realmin);
%reduce spread spectrum by 20dB to produce masking threshold
ThrBdB = XBSdB - SMR;
%combine with thr in quiet
ThrBdB = max(ThrBdB, ThrQuietBdB);

if (plt)
    fseq = [1:Nbands];
    plot(fseq, XBdB, 'k', fseq, XBSdB, 'b', fseq, ThrBdB, 'r');
    title('Signal Spectrum');
    ylabel('dB');
    xlabel('Bark')
    grid
    legend('Spectrum', 'Spread Spec', 'Threshold');
    pause
end

%convert to magnitude and Hz
XHzM = bark2hz(XBM, N2, band_binWidth);
ThrBM = 10.^(ThrBdB/20);
ThrHzM = bark2hz(ThrBM, N2, band_binWidth);

if (plt)
    fseq = [0:N2-1]/N2*fs/2;
    loglog(fseq, XHzM, 'k', fseq, ThrHzM, 'r');
    title('Signal Spectrum');
    ylabel('Magnitude');
    xlabel('Hz')
    grid
    legend('Spectrum', 'Threshold');
    pause
end





