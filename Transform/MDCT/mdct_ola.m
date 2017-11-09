function mdct_ola(ifile)
%
%mdct analysis/synthesis with 50% overlap
%ifile in input *.wav file
%ofile is ifile_out.wav

global plt quantize verboseBlock
global fs N SMR X
global Zmax bands_per_bark Nbands N2 bark_group

global mThrHzM


%coder parameters
%N = 2048; %window length
N = 1024; %to match MPEG-1 Psychoacoustic Model 1
slow_mdct = 0; 

%dont need to change these
N2 = N/2; %block length
Nbins = N2;

n = length(ifile);
base = ifile(1:n-4); %clip off '.wav' extension 
ofile = [base,'_out.wav'];

[s, fs] = audioread(ifile);
fprintf('Sampling rate is %d\n', fs);
[slen,c] = size(s);
if (c > 1)
    s = s(:,1); %take only left channel
end
n = mod(slen, N2);
%padd with 1 block before
%padd after to fill to integer number of blocks plus one more
s1 = [zeros(N2,1);s;zeros(N2-n + N2,1)];
[slen1, c] = size(s1);

nblocks = (slen1-N)/N2;
w = sin(pi*([1:N]-0.5)'/N); %window
%plot window
if plt > 0
    plot(w)
    grid
    title('Window')
    xlabel('Samples')
    ylabel('Amplitude')
    pause
end

s_out = zeros(slen, 1); %decoded signal
yprev = zeros(N2,1);
n1 = 0;
n2 = 0;
fprintf('%d blocks\n', nblocks);


perceptual_model(-96);
thrMat = zeros(Zmax * bands_per_bark, nblocks+1);
bitsPerBlock = zeros(nblocks+1, 1);


tic;
for bno = 0:nblocks
    %fprintf('block %d\n', bno);
    ns = [n1+1:n1+N];
    x = s1(ns);
    
    %window
    wx = w.*x;
    
    %mdct
    if (slow_mdct)
        n0 = (N2/2 + 1/2);
        X = mdct_mat(wx, n0);
    else
        X = mdct_fft2(wx);
    end
    
    %diagnostic plots
    if plt > 0 && bno > plt
        fseq = [1:N/2]/(N/2)*(fs/2);
        plot(fseq(1:64),X(1:64));
        grid
        title(['MDCT Spectrum for block ', num2str(bno)])
        ylabel('Amplitude')
        xlabel('Frequency')
        reply = input('CR for more plots, else 0: ');
        if ~isempty(reply)
            plt = reply;
        end
    end    
    
    %
    %quantize X and transmit
    %
    if quantize
        mThrHzM = perceptual_model(wx); % length(mthrBdB) = bands_per_bark * Zmax
        band_group = bark_group .* bands_per_bark;
        [group_SMRdB, binEndIdxs] = qx244_bandSMR(X, mThrHzM, band_group);
        group_bits = qx244_SMRdB2nBits(group_SMRdB);
        % quantize each Bark group and save to trasmit package.
        binStartIdx = 1;
        numGroup = length(binEndIdxs);
        quantizedX = zeros(N2,1);
        info = zeros(numGroup,5);
        for i = 1:numGroup
            idx = binStartIdx:binEndIdxs(i);
            nBits = group_bits(i);
            nbins = 2^nBits - 1;
            [XQ, Qidx, Qdelta] = quantMDCT(X(idx), nbins);
            info(i,:) = [i, nBits,nbins, Qdelta, length(idx)];
            if verboseBlock && plt
                stem(idx, Qidx);
                title('Qidx to transmit');
                xlabel('bin number');
                ylabel('Qidx');
                grid;
                disp(info);
                pause;
                plot(idx, X(idx), 'k', idx, XQ, 'r');
                title('Original MDCT Coeffs vs. Quantized Coeffs');
                xlabel('bin number');
                ylabel('Magnitude');
                legend('MDCT', 'Quantized');
                pause;
            end
            quantizedX(idx) = XQ;
            binStartIdx = binEndIdxs(i) + 1;
        end
        bitsNeededForThisBlock = 14 * 8 + (info(:,2)' * info(:,5)); % 14 groups 8 bit each for SF, the inner product of bits/bin and number of bins in each group.
        bitsPerBlock(bno+1) = bitsNeededForThisBlock;
        if plt
            plot(1:N2, X, 'k', 1:N2, quantizedX, 'r');
            title('Original MDCT Coeffs vs. Quantized Coeffs');
            xlabel('bin number');
            ylabel('Magnitude');
            legend('MDCT', 'Quantized');
            pause;
            bitsNeededForThisBlock
        end
        X = quantizedX;
    end
    
    %inverse mdct
    if (slow_mdct)
        n0 = (N2/2 + 1/2);
        y = imdct_mat(X, n0);
    else
        y = imdct_fft2(X);
    end

    %window
    wy = 2*w.*y;
    
    %overlap-add
    z = yprev+wy(1:N2);
    %save second half of this block to use next block
    yprev = wy(N2+1:N);
    
    if bno > 0
        %save output
        s_out(n2+1:n2+N2) = z; %output overlap portion
        n2 = n2+N2; %advance output sample counter
    end

    %advance by 50% of block length
    n1 = n1+N2; %advance input sample counter
    %bno = bno+1; %advance block counter
end
toc;
%0.66 seconds with mdct_fft2
%168.55 seconds with mdct_mat (255 times slower)

%clip off last fraction of a block in coded output
s_out = s_out(1:slen);

%calcualate SNR
e = s - s_out;
spow = s'*s;
epow = e'*e;
snr = 10*log10(spow/epow);
fprintf('SNR is %f dB\n', snr);

%write out coded result
audiowrite(ofile, s_out, fs);

if (1)
    figure(1)
    %plot original and reconstructed waveform
    ns = 1:slen;
    plot(ns, s, 'b', ns, s_out, 'r');
    legend('Original', 'Reconstructed');
    title('Analysis/Synthesis for MDCT and 50% overlap')
    ylabel('Amplitude')
    xlabel('Sample')
    grid
    pause;
    plot(bitsPerBlock);
    title('bits per block');
    ylabel('number of bits');
    xlabel('block number');
    grid
    slenInSec = slen / fs;
    overAllBitRate = sum(bitsPerBlock)/ slenInSec / 1024 % kb/sec
end
    