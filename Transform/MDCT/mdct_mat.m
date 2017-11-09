function X = mdct_mat(x, n0)

[r,c] = size(x);
N = r; % N is window length, N/2 is block length

X = zeros(N/2,1);
for k = 0:N/2-1
    xcos = cos(2*pi/N*([0:N-1]'+n0)*(k+1/2));
    X(k+1) = sum(x.*xcos);
end
