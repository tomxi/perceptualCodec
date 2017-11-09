function spread_ftn = init_spread(Zmax, bands_per_bark)
%function spread_ftn = init_spread(Zmax, bands_per_bark)
%initialize spreading function

global plt

%Schroeder spreading function, Bosi page 185
Nbands = Zmax*bands_per_bark;
dz = -Zmax:1/bands_per_bark:Zmax;
spread_ftn_dB = 15.81 + 7.5*(dz + 0.474) - 17.5*(1 + (dz + 0.474).^2).^0.5;
%magnitude
spread_ftn = 10.^(spread_ftn_dB/20);
%normalize to unity gain
spread_ftn = spread_ftn/sum(spread_ftn);
if plt
    plot([-Nbands:Nbands], spread_ftn_dB); grid
    title('Spreading Function')
    xlabel('Bark')
    ylabel('dB')
    axis([-Nbands, Nbands, -100, 0]);
    pause
end
return

