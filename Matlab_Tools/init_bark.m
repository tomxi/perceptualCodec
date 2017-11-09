function [band2bin, band_binWidth, Zmax] = init_bark(fs, N, bands_per_bark)
%function [band2bin, band_width] = init_bark(fs, N, bands_per_bark)
%
%initialize helper arrays for Hz, Bark conversions
%fs is sampling frequency
%N is bins from DC up to Nyquist in fft

global plt

%Zwicker Bark partitions, Bosi page 182
f = ([1:N]'/N)*(fs/2);
z = 13*atan(0.76*f/1000) + 3.5*atan((f/7500).^2);
if plt
    plot(f, z)
    grid
    title('Bark Scale')
    xlabel('Hz')
    ylabel('Bark')
    pause
end
%
Zmax = ceil(max(z));
Nbands = Zmax*bands_per_bark;
band2bin = zeros(Nbands, 1);
band_binWidth = zeros(Nbands, 1);
kprev = 0;
k = 1;
for i = 1:Nbands
    while ( z(k) < i/bands_per_bark )
        k = k+1; 
        if (k == N)
            break
        end
    end
    band2bin(i) = k;
    band_binWidth(i) = k - kprev;
    kprev = k;
end

end