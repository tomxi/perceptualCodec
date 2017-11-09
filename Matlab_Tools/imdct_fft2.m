function x = imdct_fft2(X)

[r,c] = size(X);
if c>1
    fprintf('Must be Mx1 column vector\n');
    return
end
if mod(r,2) ~= 0
    fprintf('Length of x must be even\n')
end

N = r;

tw_1 = exp(j*(pi/N)*[0:2*N-1]'*(N/2+1/2));
tw_2 = exp(j*pi/(2*N)*([0:2*N-1]'+ N/2 + 1/2 ));

X2 = [X;-flipud(X)];

X3 = X2.*tw_1;

x3 = ifft(X3);

x4 = x3.*tw_2;

x = real(x4);

