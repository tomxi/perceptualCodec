function [xq, xbins, xidx] = quantThr(xdB, xdelta)
%function [xq, xbins, xidx] = quantThr(xdB, xstep)

%positve value only mid-tread quantizer

xidx = round(xdB/xdelta);
xbins = max(xidx);
xq = xidx*xdelta;
