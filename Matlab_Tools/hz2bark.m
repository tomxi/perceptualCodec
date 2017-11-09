function XMagBark = hz2bark(XMagHz, band_binWidth)
%function XMagBark = hz2bark(XMagHz, band_binWidth)
%converts from Hz scale to Bark scale
%
%Xhz is magnitude spectrum on Hz scale
%band_binWidth is Bark helper varaible

XPowHz = XMagHz.*XMagHz; %convert to power spectrum
Nbands = length(band_binWidth);
XPowBark = zeros(Nbands, 1);
k = 0;
for i=1:Nbands
    XPowBark(i) = sum(XPowHz(k+1:k+band_binWidth(i)));
    k = k + band_binWidth(i);
end

XMagBark = sqrt(XPowBark);
end

