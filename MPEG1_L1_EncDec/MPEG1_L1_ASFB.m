function MPEG1_L1_ASFB(ifile, ofile)
% function MPEG_LI_ASFB(ifile, ofile)
%
% ifile must be 1-channel *.wav audio file with 44.1 kHz sampling rate
%
% MPEG Layer I Analysis/Synthesis Filterbank
%

debug_plots = 1;

%read input audio file
[x, fs] = audioread(ifile);
[nsamp, nchn] = size(x);
if (nchn > 1)
    x = x(:,1); %select only left to make mono
end
%zero padd for end-of-file processing
x = [x; zeros(1024,1)];
%initialize y
y = zeros(length(x), 1);

% Load filterbank coefficient tables.
C = Table_analysis_window;
D = Table_synthesis_window;

% Initialize state variables
X = zeros(512, 1);
V = zeros(1024, 1);

% Compute matrixing functions
% Analysis
for i = 0:31
    for k = 0:63
        M(i+1, k+1) = cos( (2*i + 1)*(k - 16)*pi/64 );
    end
end
% Synthesis
for i = 0:63
    for k = 0:31
        N(i+1, k+1) = cos( (16 + i)*(2*k + 1)*pi/64 );
    end
end

% Debugging plots
if (debug_plots)
    close all
    figure(1); ns=1:512; plot(ns, C, 'b', ns, D/32, 'r');
    figure(2); mesh(M);
    figure(3); mesh(N);
    pause
    close all
    % Debugging variables
    Y = zeros(512, 1);
end

% Process input file
for i = 0:384:length(x)-384
    fprintf('%d\n', i); %sample
 
    % Next block of 384 = 12*32 input
    S = zeros(12, 32);
    for j = 0:(12-1)
        x1 = x(i+j*32+1:i+j*32+32);
        [S(j+1,:), X] = analysis_filterbank(x1, C, M, X);
    end
    
    for j = 0:(12-1)
        [y1, V] = synthesis_filterbank(S(j+1,:), D, N, V);
        y(i+j*32+1:i+j*32+32) = y1;
    end
    
    % Debugging plots
    if (debug_plots)
        subplot(3,1,1); plot(x1); grid;
        subplot(3,1,2); mesh(S); grid
        subplot(3,1,3); plot(y1); grid
        pause
    end
end

%write output audio file
audiowrite(ofile, y, fs);

close all
nsx = 1:nsamp;
nsy = nsx+512+1-32;
plot(nsx, x(nsx), 'b', nsx, y(nsy), 'r');
legend('Original', 'Decoded');
ylabel('Sample value');
xlabel('Sample number');

d = x(nsx) - y(nsy);
fprintf('SNR: %5.1f dB\n', snr(x(nsx), d));




