function [group_SMRdB, binEndIdxs] = qx244_bandSMR(X, mThrHzM, band_group)
%function group_SMR = qx244_bandSMR(X, mThrHzM, band_group)
%   Takes MDCT coeffs and masking threshold, and a vector discribing 
%   how bins are grouped into bands. Outputs SMR for each of the
%   band_groups. Output size is the same as band_group.
%   

global plt
global band2bin

[r,groups] = size(band_group);
group_SMR = zeros(r,groups);
binEndIdxs = zeros(r,groups);


binStartIdx = 1;
for gno = 1:groups
    binEndIdx = band2bin(band_group(gno));
    binEndIdxs(gno) = binEndIdx;
    binIdxs = binStartIdx:binEndIdx; % get the list of bins associated with the gno (group number)
    group_SMR(gno) = max( abs(X(binIdxs)) ) / min( mThrHzM(binIdxs) );
    binStartIdx = binEndIdx + 1; % increment bin idx.
end

group_SMRdB = 20*log10(group_SMR);

if plt
    stem(band_group, group_SMRdB);
    grid
    title('Group SMR in dB')
    xlabel('band number')
    ylabel('SMR in dB')
    pause
end

end

