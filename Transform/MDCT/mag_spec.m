function X = mag_spec(x)
%calculate magnitude spectrum using odly-stacked fft
%x is real values

[r,c] = size(x);
N = r; %input length
if (c > 1)
    x = x(:,1); %use first channel only
end

%twiddle to shift from even to odd stacking
x = x.*exp(i*2*pi*0.5*[0:N-1]'/N);
X = abs(fft(x));

%return unique values, but skip Nyquist value
X = X(1:N/2);
