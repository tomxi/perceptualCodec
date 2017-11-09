function mdct_ola(ifile, plt_value)
%
%mdct analysis/synthesis with 50% overlap
%ifile in input *.wav file
%ofile is ifile_out.wav

MDCT_Only = 0;

global plt
global fs N SMR
global Zmax bands_per_bark Nbands N2
global band_binWidth spread_ftn ThrQuietBdB

plt = plt_value;

%coder parameters
N = 2048; %window length
N2 = N/2; %block length
Nbins = N2;
SMR = 15;

Tdelta = 1.5; %stepsize of 1.5 dB is emperical
Xbins = 63;


n = length(ifile);
base = ifile(1:n-4); %clip off '.wav' extension 
ofile = [base,'_out.wav'];

[s, fs, nbits] = wavread(ifile);
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

if (~MDCT_Only)
    %initialize perceptual model
    perceptual_model(-96);
end

s_out = zeros(slen, 1); %decoded signal
yprev = zeros(N2,1);
n1 = 0;
n2 = 0;
fprintf('%d blocks\n', nblocks);
tic;
for bno = 0:nblocks
    %fprintf('block %d\n', bno);
    ns = [n1+1:n1+N];
    x = s1(ns);
    
    %window
    wx = w.*x;
    
    %mdct
    %using matrix multiply
    %n0 = (N2/2 + 1/2);
    %X = mdct_mat(wx, n0);
    %using fast fft kernel
    X = mdct_fft2(wx);
    
    if plt > 0 && bno > plt
        fseq = [1:N/2]/(N/2)*(fs/2);
        plot(fseq(1:64),X(1:64));
        grid
        title(['MDCT Spectrum for block ', num2str(bno)])
        ylabel('Amplitude')
        xlabel('Frequency')
        reply = input('CR for more plots, else 0: ');
        %fprintf('%d\n', reply);
        if ~isempty(reply)
            plt = reply;
        end
    end    
    
    if (MDCT_Only)
        %don't quantize
        XQ = X;
    else
        %compute Thr for this block
        ThrBdB = perceptual_model(wx);

        %quantize Thr using unity stepsize single-sided mid-tread quantizer
        [ThrBdBQ, Tbins, Tidx] = quantThr(ThrBdB, Tdelta);
        if (plt)
            stride = 1/bands_per_bark;
            ns = [stride:stride:Zmax]';
            plot(ns, ThrBdBQ, 'r');
            title('Quantized Threshold');
            ylabel('dB');
            xlabel('Bark')
            grid
            pause
        end
        %convert to magnitude and Hz
        ThrBMQ = 10.^(ThrBdBQ/20);
        ThrHzMQ = bark2hz(ThrBMQ, N2, band_binWidth);
      

        %quantize MDCT coefficients using quantized Thr
        %Thr = stepsize/2 = maximum allowable noise level
        XScaledHzM = X./(2*ThrHzMQ);
        %quantize with unity stepsize symmetric mid-tread quantizer
        [XQ, Xidx, Xdelta] = quantMDCT(XScaledHzM, Xbins);

        if (plt)
            fseq = [0:N2-1]/N2*fs/2;
            loglog(fseq, abs(X), 'k', fseq, abs(XQ), 'b', fseq, ThrHzMQ, 'r');
            title('Signal Spectrum');
            ylabel('Magnitude');
            xlabel('Hz')
            grid
            legend('Spectrum', 'Quantized Scaled Spectrum', 'Threshold');
            pause
        end
     
        %
        %need to transmit Tidx, Tdelta, Xidx, Xdelta
        %
        %count bits needed
        count1 = ceil(log2(Tbins)) * length(Tidx);
        count1 = count1 + 32; %just send Tdelta as a float
        count2 = ceil(log2(Xbins)) * length(Xidx);
        count2 = count2 + 32; %just send Xdelta as a float
        fprintf('block %d: %d bits (%d %d)\n', bno, count1+count2, count1, count2);
        
        %reconstruct spectrum from mid-tread quantizer
        XQ = Xidx*Xdelta;
        
        %use quantized scale factors
        %scale spectrum by Scale Factors
        XQ = XQ.*(2*ThrHzMQ);
        
        if (plt)
            loglog(fseq, abs(X), 'k', fseq, abs(XQ), 'b', fseq, ThrHzMQ, 'r');
            title('Signal Spectrum');
            ylabel('Magnitude');
            xlabel('Hz')
            grid
            legend('Spectrum', 'Quantized Spectrum', 'Threshold');
            pause
        end
    end
    
    %inverse mdct
    %using matrix multiply
    %n0 = (N2/2 + 1/2);
    %y = imdct_mat(XQ, n0);
    %using fast fft kernel
    y = imdct_fft2(XQ);

    %window
    wy = 2*w.*y;
    
    %overlap-add
    z = yprev+wy(1:N2);
    %save second half of this block to use next block
    yprev = wy(N2+1:N);
    
    if bno > 0
        s_out(n2+1:n2+N2) = z;
        n2 = n2+N2;
    end

    %advance by 50% of block length
    n1 = n1+N2;
    bno = bno+1;
end
toc
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
wavwrite(s_out, fs, ofile);

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
end
    