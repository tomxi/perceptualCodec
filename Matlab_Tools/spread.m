function XSMB = spread(XMagBark)
%function XBSmag = spread(XMagBark)
%XMagBark is bark-scale magnitude spectrum
%XSMB is spread bark-scale magnitude spectrum
%
%spread Bark spectrum by spreading function

%extend XMagBark to right so spreading function can "run off the end"
Nbands = length(XMagBark);
XMBx = [XMagBark; zeros(Nbands,1)];
XSMB = filter(spread_ftn, 1, XMBx);
%clip at left side
XSMB = XSMB(Nbands+1:end); 
return
