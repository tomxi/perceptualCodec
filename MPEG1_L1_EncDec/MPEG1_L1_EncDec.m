function MPEG1_L1_EncDec(ifile, ofile, bit_rate)
% function MPEG_Layer_I(ifile, ofile, bit_rate)
%
% ifile must be *.wav audio file with 44.1 kHz sampling rate
% ofile is reconstructed *.wav output file
% bit_rate is bits/second
%
% MPEG-1 Layer I Encoder/Decoder
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

%calculate bits per block
%bit_rate = 192000; %bits/second
b = 12*32*bit_rate/fs;
block_bits = floor(b);
%frac_bits are accumulated to signal padding bit
frac_bits = b - block_bits;
accum_frac_bits = 0;

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

% Debugging plots
if (DRAW)
    close all
    figure(1); ns=1:512; plot(ns, C, 'b', ns, D/32, 'r'); grid; title('Analysis window'); 
    figure(2); mesh(M); title('Analysis Matrixing');
    figure(3); mesh(N); title('Synthesis Matrixing');
    pause
    close all
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

    %%% Code
    %
    %coding_mode = 'noise';
    coding_mode = 'coded';
    switch lower(coding_mode)
        case 'noise'
            %%% Add noise at masking level
            SMRampl = 10.^(SMRsb/20);
            for k = 1:32
                sig_rms = rms(S(:,k));
                noise = (2*(rand(12,1)-0.5)); % [-1 1] uniform noise 
                noise = noise*sig_rms/SMRampl(k); % decreased by SMR wrt sig_rms
                S(:,k) = S(:,k) + noise;
            end
        case 'coded'
            %%% Quantizer per Layer I spec, but don't construct bitstream
            % Treat fractional bit count
            if (accum_frac_bits > 32)
                this_block_bits = block_bits + 32; 
                accum_frac_bits = accum_frac_bits - 32;
                % fprintf('Padding bit at %d\n', i);
            else
                this_block_bits = block_bits;
            end
            S = Layer1_code(S, scf, SMRsb, block_bits);
         otherwise
            fprintf('Unknown mode: %s\n', mode);
    end
    %
    %%% End code
    
    
    %%% Synthesis Filterbank
    %
    for j = 0:(12-1)
        [y1, Vstate] = synthesis_filterbank(S(j+1,:), D, N, Vstate);
        y(i+j*32+1:i+j*32+32) = y1;
    end
    
    %%% accumulate fractional bits
    accum_frac_bits = accum_frac_bits + frac_bits;
    
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




