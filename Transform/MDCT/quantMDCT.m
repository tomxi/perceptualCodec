function [XQ, Qidx, Qdelta] = quantMDCT(XS, nbins)
%function [XQ, Xidx] = quantMDCT(XS, nbins)
%quantize MDCT Coefs using symmetric mid-tread quantizer

XSmax = max(XS);
XSmin = min(XS);
Qmax = max(XSmax, -XSmin);
if (Qmax < realmax)
    Qdelta = 2*Qmax/nbins;
else
    %if Qmax is inf
    Qdelta = 1;
end
if (Qdelta > 0)
    Qidx = round(XS/Qdelta);
else
    %to avoid divide by zero
    Qidx = zeros(length(XS), 1);
end
XQ = Qidx*Qdelta;