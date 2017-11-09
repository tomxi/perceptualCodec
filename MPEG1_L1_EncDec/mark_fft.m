function mark_fft(ifile, mfile, ofile)
% function mark_fft(ifile, mfile, ofile)

%read input audio file
[x, fs] = audioread(ifile);
[sig_samp, nchn] = size(x);
if (nchn > 1)
    x = x(:,1); %select only left to make mono
end

%read mark file
[mark, fs] = audioread(mfile);
[mark_samp, nchn] = size(x);
if (nchn > 1)
    mark = mark(:,1); %select only left to make mono
end
%make same length as audio
if ( mark_samp < sig_samp )
    fprintf('mark file too short\n');
    return;
else
    mark = mark(1:sig_samp);
end

% Initialize constants
Common;

% Load tables
[TH, Map, LTq] = Table_absolute_threshold(1, fs, 128); % Threshold in quiet
CB = Table_critical_band_boundaries(1, fs);

%parameters
N = 512;  %transform length
olf = 0.5; %fft 50% window overlap
Nbins = N/2 + 1; %DC through Nyquist
N1 = N*olf;%overlap
N2 = N-N1; %new samples each block

%create window 
%window may have flat top if less than 50% overlap
if N1 == N2
    w = sin(pi*([1:N]-0.5)'/N); %50% overlap
else
    w1 = sin(pi/2*([1:N1]-0.5)'/N1);
    w2 = ones(N-2*N1, 1);
    w = [w1; w2; flipud(w1)];
end
%plot window
if (0)
    plot(w)
    grid
    title('Window')
    xlabel('Samples')
    ylabel('Amplitude')
    pause
end

nblocks = floor(sig_samp/N2);
yprev = zeros(N1,1);
z = zeros(nblocks*N2, 1);
n = 0;
for bno = 1:nblocks
    if (mod(bno, 100) == 0)
        fprintf('Block %d of %d\n', bno, nblocks);
    end
    
    %select signal block
    if (bno == 1)
        %first N2 has nothing to overlap with, so make first block zero
        %and start signal at second block
        x1 = [zeros(N2, 1); x(1:N2)]; %signal
        m1 = [zeros(N2, 1); mark(1:N2)]; %mark
    else
        x1 = x(n+1:n+N); %signal
        m1 = mark(n+1:n+N); %mark
    end
    
    %window
    wx = w.*x1; %signal
    wm = w.*m1; %mark
    
    %fft
    X = fft(wx);
    M = fft(wm);

    %save only unique part of spectrum
    X1 = X(1:N/2+1);
    M1 = M(1:N/2+1);
   
    %plot spectrum
    if (0)
        XdB = 20*log10(abs(X1));
        XdB = XdB - max(XdB);
        MdB = 10*log10(abs(M1));
        fseq = [0:N/2]'/(N/2)*fs/2;
        plot(fseq, XdB, 'b', fseq, MdB, 'r');
        grid
        xlabel('Frequency');
        ylabel('Magnitude, dB');
        title(['Spectrum, block ', num2str(bno)]);
        legend('Signal', 'Mark');
        reply = input('CR for more plots, else 0: ');
        %fprintf('%d\n', reply);
        if ~isempty(reply)
            plt = reply;
        end
    end

if (1)
    %%% Psychoacoustic analysis.
    %
        % Compute the FFT for time frequency conversion [1, pp. 110].
        % don't use this function, instead use code below
        %X = FFT_Analysis(s, n); % same as "x"
        % Prepare the Hanning window
        h = sqrt(8/3) * hanning(512, 'periodic');

        % Power density spectrum
        Xpm = max(20 * log10(abs(fft(x1 .* h)) / FFT_SIZE), MIN_POWER);

        % Normalization to the reference sound pressure level of 96 dB
        Delta = 96 - max(Xpm);
        Xpm = Xpm + Delta;
        
        % Determine the sound pressure level in each  subband [1, pp. 110].
        % don't use this function, instead use code below
        % Lsb = Sound_pressure_level(X, scf);
        % This is now SPL per FFT bin
        Lbin = Xpm(1:N/2+1);
        
        % Find the tonal (sine like) and non-tonal (noise like) components
        % of the signal [1, pp. 111--113]
        [Flags Tonal_list Non_tonal_list] = Find_tonal_components(Xpm, TH, Map, CB);

        % Decimate the maskers: eliminate all irrelevant maskers [1, pp. 114]
        [Flags Tonal_list Non_tonal_list] = ...
          Decimation(Xpm, Tonal_list, Non_tonal_list, Flags, TH, Map);

        % Compute the individual masking thresholds [1, pp. 113--114]
        [LTt, LTn] = ...
          Individual_masking_thresholds(Xpm, Tonal_list, Non_tonal_list, TH, Map);

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
        % don't use this function, instead use LTg directly
        %LTmin = Minimum_masking_threshold(LTg, Map);
        % map LTg back to FFT bins
        LTbin = zeros(N/2+1, 1);
         for k = 1:N/2
            LTbin(k) = LTg(Map(k));
         end
        LTbin(N/2+1) = LTbin(N/2); %Nyquist

        % Align masking threshold with signal
        X1dB = 20*log10(abs(X1));
        X1avg = mean(X1dB);
        LTavg = mean(LTbin);
        THRbin = LTbin - LTavg + X1avg - 20; 
       
        if (DRAW)
            ns = 1:N/2+1;
            figure(2)
            plot(ns, X1dB, 'b', ns, THRbin, 'r'); 
            legend('Signal', 'Threshold');
            title('Signal and Threshold');
            xlabel('FFT bin'); ylabel('dB'); 
            pause;
            close(gcf);
        end
        
    %   
    %%% End of psychoacoustic analysis.
    
    %%% insert watermark
    %
        % Convert THRbin to amplitude value
        THRampl = 10.^(THRbin/20);
        
        % Add mark at masking level
        X1 = X1 + M1.*THRampl;
    %
    %%% End mark  
end
        
    %re-construct from conjugate symmetric part, X1
    Xcs = flipud(conj(X1(2:N/2)));
    X = [X1; Xcs];
    
    %inverse fft
    y = ifft(X);
    y = real(y); %just in case ifft has some errors
    
    %window
    wy = w.*y;
    
    %overlap-add
    if (N1 == N2)
        z1 = yprev+wy(1:N1); %50% overlap
    else
        z1 = [(yprev + wy(1:N1)); wy(N1+1:N-N1)];
    end
    %save second portion of this block to use next block
    yprev = wy(N-N1+1:N);

    if (bno == 1)
        ; %don't save and don't advance input 
    else
        %save output of analysis/synthesis system
        z(n+1:n+N2) = z1;

        %advance input index by block length
        n = n+N2;
    end
end

%last block is stuck in overlap buffer, so don't use it
len = (nblocks-1)*N2;
x = x(1:len);
z = z(1:len);
%write output file
audiowrite(ofile, z, fs);

%calcualate SNR
e = x - z;
xpow = x'*x;
epow = e'*e;
snr = 10*log10(xpow/epow);
fprintf('SNR is %f dB for total file\n', snr);

%plot original and reconstructed waveform
ns = [1:len]';
plot(ns, x, 'b', ns, z, 'r');
grid
title('Analysis/Synthesis for FFT and 50% overlap')
ylabel('Amplitude')
xlabel('Sample')
grid
   