function XMagHz = bark2hz(XMagBark, Nbins, band_binWidth)
%function XMagHz = bark2hz(XMagBark, Nbins, band_binWidth)
%converts from Bark scale to Hz scale
%XMagBark is magnitude spectrum on bark scale
%Nbins is number of bins in Hz scale

XPowBark = XMagBark.*XMagBark;
XPowHz = zeros(Nbins, 1);
Nbands = length(band_binWidth);
k = 0;
for i = 1:Nbands
    n = band_binWidth(i);
    XPowHz(k+1:k+n) = XPowBark(i)/n;
    k = k+n;
end

XMagHz = sqrt(XPowHz);
