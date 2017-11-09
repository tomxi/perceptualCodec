function [y1, V] = synthesis_filterbank(S, D, N, V)

% Shift V
V(1024:-1:64+1) = V(1024-64:-1:1);

% Matrixing
V(64:-1:1) = zeros(64,1);
for i = 1:64
    for k = 1:32
        V(i) = V(i) + N(i, k) * S(k);
    end
end

% Build U vector
U = zeros(512,1);
for i = 0:7
    for j = 1:32
        U(i*64 + j) = V(i*128 + j);
        U(i*64 + 32 + j) = V(i*128 + 96 + j);
    end
end

% Window
W = U .* D;

% Calculate output samples
% Standard uses S, but here we use y1
y1 = zeros(32,1);
for j = 1:32
    for i = 0:15
        y1(j) = y1(j) + W(j + 32*i);
    end
end

return
