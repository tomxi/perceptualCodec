function mark_MPEG1(ifile, mfile, ofile)
% function mark_MPEG1(ifile, mfile, ofile)
%
% CDMA watermark using MPEG-1 Layer I filterbank and
% Psychoacoustic model
%

%read input audio file
[x, fs] = audioread(ifile);
[nsamp, nchn] = size(x);
if (nchn > 1)
    x = x(:,1); %select only left to make mono
end
%zero padd for end-of-file processing
x = [x; zeros(2048,1)];
%initialize y
y = zeros(length(x), 1);

%read mark file
[mark, fs] = audioread(mfile);
[mark_samp, nchn] = size(x);
if (nchn > 1)
    mark = mark(:,1); %select only left to make mono
end
%make same length as audio
if ( mark_samp < length(x) )
    fprintf('mark file too short\n');
    return;
else
    mark = mark(1:length(x));
end

% Initialize constants
Common;

% Load filterbank coefficient tables.
C = Table_analysis_window;
D = Table_synthesis_window;
% Load other tables
[TH, Map, LTq] = Table_absolute_threshold(1, fs, 128); % Threshold in quiet
CB = Table_critical_band_boundaries(1, fs);

% Initialize state variables
Xstate = zeros(512, 1);
Mstate = zeros(512, 1);
Vstate = zeros(1024, 1);

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

% Process blocks of input file
nblocks = ceil(nsamp/384); %process last block
bno = 0;
for i = 0:384:nblocks*384;
    if (mod(bno,100) == 0) %every 100 blocks
        fprintf('%d\n', i); %sample number
    end
    bno = bno+1;
 
    %%% Analysis Filterbank
    % Next block of 384 input samples generates
    % 12*32 subband samples
    S = zeros(12, 32);
    for j = 0:(12-1)
        x1 = x(i+j*32+1:i+j*32+32);
        [S(j+1,:), Xstate] = analysis_filterbank(x1, C, M, Xstate);
    end

    % Filter mark data
    S_mark = zeros(12, 32);
    for j = 0:(12-1)
        m1 = mark(i+j*32+1:i+j*32+32);
        [S_mark(j+1,:), Mstate] = analysis_filterbank(m1, C, M, Mstate);
    end

    
    %%% Scalefactor calculation [1, pp. 70].
    scf = Scale_factors(S);
   
    %%% Psychoacoustic analysis.
    %
        % Compute the FFT for time frequency conversion [1, pp. 110].
        X = FFT_Analysis(x, i+256); % 256 is filterbank delay

       % Determine the sound pressure level in each  subband [1, pp. 110].
       Lsb = Sound_pressure_level(X, scf);

       % Find the tonal (sine like) and non-tonal (noise like) components
       % of the signal [1, pp. 111--113]
       [Flags Tonal_list Non_tonal_list] = Find_tonal_components(X, TH, Map, CB);

       % Decimate the maskers: eliminate all irrelevant maskers [1, pp. 114]
       [Flags Tonal_list Non_tonal_list] = ...
          Decimation(X, Tonal_list, Non_tonal_list, Flags, TH, Map);

       % Compute the individual masking thresholds [1, pp. 113--114]
       [LTt, LTn] = ...
          Individual_masking_thresholds(X, Tonal_list, Non_tonal_list, TH, Map);

       % Compute the global masking threshold [1, pp. 114]
       LTg = Global_masking_threshold(LTq, LTt, LTn);

       if (DRAW)
          disp('Global masking threshold');
          hold on;
          plot(TH(:, INDEX), LTg, 'k--');
          hold off;
          title('Masking components and masking thresholds.');
       end

       % Determine the minimum masking threshold in each subband [1, pp. 114]
       LTmin = Minimum_masking_threshold(LTg, Map);
       if (DRAW)
          figure; plot(LTmin); title('Minimum masking threshold');
          xlabel('Subband number'); ylabel('dB'); pause;
       end

       % Compute the singal-to-mask ratio
       SMRsb = Lsb - LTmin;
    %   
    %%% End of psychoacoustic analysis.

    %%% insert watermark
    %
        %%% Add mark signal at masking level
        % Reduce SMRsb in order to make mark work
        SMRsb = SMRsb - 20; %dB
        
        SMRampl = 10.^(SMRsb/20);
        for k = 1:32
            sig_rms = rms(S(:,k));
            S_mark(:,k) = S_mark(:,k)*sig_rms/SMRampl(k); % decreased by SMR wrt sig_rms
            S(:,k) = S(:,k) + S_mark(:,k);
        end
    %
    %%% End mark   
    
    %%% Synthesis Filterbank
    %
    for j = 0:(12-1)
        [y1, Vstate] = synthesis_filterbank(S(j+1,:), D, N, Vstate);
        y(i+j*32+1:i+j*32+32) = y1;
    end
    
    % Debugging plots
    if (DRAW)
        subplot(3,1,1); plot(x1); grid;
        subplot(3,1,2); mesh(S); grid
        subplot(3,1,3); plot(y1); grid
        pause
        close all
    end
end

close all
ns = 1:nsamp;
nsy = ns+512+1-32; % Delay of y WRT x
x = x(ns);
y = y(nsy);

%write output audio file
audiowrite(ofile, y, fs);

plot(ns, x, 'b', ns, y, 'r');
legend('Original', 'Decoded');
ylabel('Sample value');
xlabel('Sample number');

d = x - y;
fprintf('SNR: %5.1f dB\n', snr(x, d));




