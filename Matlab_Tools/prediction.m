function prediction(ifile)
%linear prediction
%
%ifile in input *.wav file
%ofile is ifile_residual.wav

plt = 1;

n = length(ifile);
base = ifile(1:n-4); %clip off '.wav' extension 
ofile = [base,'_residual.wav'];

[s, fs, nbits] = wavread(ifile);
fprintf('Sampling rate: %d\n', fs);
[slen,c] = size(s);
if (c > 1)
    s = s(:,1); %use channel 1
end

%set parameters
%30 ms analysis window
N = round(.030 *fs); 
n = mod(N,4);
N = N-n; %make N divisible by 4
N2 = N/2;
N4 = N/4;
P = 16; %order 16 analysis

%padd out to integer number of analysis blocks 
%and then add one more analysis block
n = N - mod(slen, N);
s = [s; zeros(n+N, 1)];
len = length(s);
nblocks = (len - N)/N2;
fprintf('%d samples, %d blocks\n', slen, nblocks);

%initialize variables
zi = zeros(P,1);
residual = zeros(len, 1);
win = hann(N);
if (0)
    subplot(1,1,1); plot(win);
    pause;
end

%for each block
n = 0;
for bno = 1:nblocks
    %fprintf('%d\n', bno);
    %select block x 
    x = s(n+1: n+N);
    
    %design lpc filter
    %window data block with Hann window
    % A = lpc(X, N)
    % A=[ 1 A(2) ... A(N+1) ], of an Nth order predictor
    % Xp(n) = -A(2)*X(n-1) - A(3)*X(n-2) - ... - A(N+1)*X(n-N)
    a = lpc(hann(N).*x, P);
    
    %inverse filter middle half of this signal block
    % Y = filter(B,A,X) filters the data in vector X with the
    % a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
    %                      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
    %filter
    ns = n+N4+1:n+N4+N2;
    x = s(ns); %middle half of block
    [residual(ns), zf] = filter(a, 1, x, zi);
    zi = zf; %for next time
    
    if (plt)
        lower = n+1+N4-N2;
        if (lower < 1)
            lower = 1;
        end
        ns = lower:n+N4+N2; %plot prev and this block
        subplot(2,1,1); plot(ns, s(ns)); grid
        xlabel('sample'); ylabel('amplitude'); title('speech signal');
        subplot(2,1,2); plot(ns, residual(ns)); grid
        xlabel('sample'); ylabel('amplitude'); title('residual');
        reply = input('CR for more plots, else 0: ');
        if ~isempty(reply)
            plt = reply;
        end
    end
    %move ahead by N2
    n = n+N2;
end

wavwrite(residual, fs, nbits, ofile);

if (1)
    ns = 1:len;
    subplot(1,1,1); plot(ns, s, 'b', ns, residual, 'r'); grid
    legend('Speech', 'Residual');
    xlabel('sample'); ylabel('amplitude'); title('speech and resicual signals');
end

