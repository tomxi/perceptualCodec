function nBits = qx244_SMRdB2nBits(SMRdB)
%function nBits = qx244_SMRdB2nBits(SMRdB)
%   convert SMR vectors to number of bits vector. Using 6dB down noise per
%   added bit rule.
len = length(SMRdB);
nBits = zeros(1, len);
for i = 1:len
    nBit = ceil(SMRdB(i) / 6); % number of bits should provide at least SMRdB down noise.
    nBits(i) = max(nBit, 1); % use at least 1 bit. This is to guard against when SMRdB is negative (past 15K Hz).
end
    