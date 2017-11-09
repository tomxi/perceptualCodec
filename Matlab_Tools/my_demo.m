%demo of perceptual model tools

%global plot variable
global plt

plt = 0;

%read signal
[sig, fs] = audioread('horn23_2.wav');

%must be mono
sig = sig(:,1);

%coder window length and block length
N = 2048;
N2 = 1024;

%step in to 25th block
n = 25*N2-1;
x = sig(n+1:n+N);

%create window
win = sin(pi*([1:N]'-0.5)/N);
plot(win);
grid
title('FFT Window');
pause;

%magnitude spectrum
XM = mag_spec(x.*win);
loglog(XM);
grid
title('Signal Hz Magnitude Spectrum');
pause

%bark spectrum
bands_per_bark = 2;
[band2bin, band_binWidth, Zmax] = init_bark(fs, N2, bands_per_bark);
XBM = hz2bark(XM, band_binWidth);
semilogy(XBM);
grid
title('Signal Bark Magnitude Spectrum');
pause

%25 Bark, 2 bands per Bark
Nbands = Zmax*bands_per_bark;
spread_ftn = init_spread(Zmax, bands_per_bark);

%spread spectrum
%pad XBM with zeros so spreading function can run to end
XSpread = filter(spread_ftn, 1, [XBM; zeros(Nbands,1)]);
%now compensate for width of spread_ftn
XSpread = XSpread(Nbands+1:end);
semilogy([.5:.5:25], XSpread);
grid
title('Signal Bark Magnitude Spread Spectrum');
pause

%convert to dB
XBdB = 20*log10(XBM);
XBSdB = 20*log10(XSpread);
%reduce spread spectrum by 20dB to produce masking threshold
ThrBdB = XBSdB - 20;
%thr in quiet
ThrQuietBM = thr_in_quiet(-96, N2, fs, band_binWidth);
ThrQuietBdB = 20*log10(ThrQuietBM);
%combine with thr in quiet
ThrBdB = max(ThrBdB, ThrQuietBdB);


%plot
fseq = [1:Nbands];
plot(fseq, XBdB, 'k', fseq, XBSdB, 'b', fseq, ThrBdB, 'r');
title('Signal Spectrum');
ylabel('dB');
xlabel('Bark')
grid
legend('Spectrum', 'Spread Spec', 'Threshold');
pause

%convert to magnitude and Hz
XHzM = bark2hz(XBM, N2, band_binWidth);
ThrBM = 10.^(ThrBdB/20);
ThrHzM = bark2hz(ThrBM, N2, band_binWidth);

%plot
fseq = [0:N2-1]/N2*fs/2;
loglog(fseq, XHzM, 'k', fseq, ThrHzM, 'r');
title('Signal Spectrum');
ylabel('Magnitude');
xlabel('Hz')
grid
legend('Spectrum', 'Threshold');
pause

%call ThrHzM "Scale Factors"
%quantize Scale Factors
Thr_min = min(ThrHzM);
Thr_max = max(ThrHzM);
fprintf('Thr min, max: %f %f\n', Thr_min, Thr_max);
Qsteps = 256;
[ThrHzMQ, Tidx, Tdelta] = quantThr(ThrHzM, 1/256);

loglog(fseq, ThrHzMQ, 'r');
title('Quantized Threshold');
ylabel('Magnitude');
xlabel('Hz')
grid
pause

XScaledHzM = XHzM./ThrHzMQ;
[XQ, Xidx, Xdelta] = quantMDCT(XScaledHzM,255);

loglog(fseq, XHzM, 'k', fseq, XQ, 'b', fseq, ThrHzMQ, 'r');
title('Signal Spectrum');
ylabel('Magnitude');
xlabel('Hz')
grid
legend('Spectrum', 'Quantized Scaled Spectrum', 'Threshold');
pause







