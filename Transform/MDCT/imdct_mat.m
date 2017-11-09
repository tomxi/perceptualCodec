function x = imdct_mat(X, n0)

[r,c] = size(X);
N = 2*r; % N is window length, N/2 = is block length

x = zeros(N,1);
for n = 0:N-1
    xcos = cos(2*pi/N*(n+n0)*([0:N/2-1]'+1/2));
    x(n+1) = sum(X.*xcos);
end
x = x*(4/N);