function [S, X]  = analysis_filterbank(x1, C, M, X)
    
% Shift X
% Higest index is oldest, lowest is newest
X(512:-1:32+1) = X(512-32:-1:1);
% Insert new sample
X(32:-1:1) = x1;

% Window
Z = C .* X;

Y = zeros(64,1);
for i= 1:64
    for j = 0:7
        Y(i) = Y(i) + Z(i + 64*j);
    end
end

% Matrix
S = zeros(1, 32);
for i = 1:32
    for k = 1:64
        S(i) = S(i) + M(i, k) * Y(k);
    end
end

return



