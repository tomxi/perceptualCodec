function fft_ola(olf, ifile)
%function fft_ola(olf, ifile)
%olf   overlap fraction, must be 1/(power of 2)
%ifile input file
%ofile is ifile_out.wav
%
%fft overlap analysis and overlap add synthesis

plt = 0;

%parameters
N = 1024;  %transform length

Nbins = N/2 + 1; %DC through Nyquist
N1 = N*olf;%overlap
N2 = N-N1; %new samples each block
fprintf('Overlap fraction is %6.3f, %d samples\n', olf, N*olf);

[s, fs] = audioread(ifile);
[slen,c] = size(s);
fprintf('Fs %d, %d seconds, %d channels\n', fs, slen/fs, c);
if (c > 1)
    s = s(:,1); %use only left channel
end
n = length(ifile);
base = ifile(1:n-4); %strip .wav extension
ofile = [base, '_out.wav'];

n = mod(slen, N2);
%padd with N1 zero samples before
%padd after to fill to integer number of blocks plus one more
%this accounts for analysis/synthesis processing pipeline
s1 = [zeros(N1,1);s;zeros(N2-n + N2,1)];
[slen1, ~] = size(s1);

%window may have flat top if less than 50% overlap
if N1 == N2
    w = sin(pi*([1:N]-0.5)'/N); %50% overlap
else
    w1 = sin(pi/2*([1:N1]-0.5)'/N1);
    w2 = ones(N-2*N1, 1);
    w = [w1; w2; flipud(w1)];
end
%plot window
if plt > 0
    plot(w)
    grid
    title('Window')
    xlabel('Samples')
    ylabel('Amplitude')
    pause
end

nblocks = (slen1-N)/N2;
fprintf('%d blocks\n', nblocks);
yprev = zeros(N1,1);
s_out = zeros(nblocks*N2, 1);
n = 0;
for bno = 0:nblocks
    %fprintf('Block %d\n', bno);
    x = s1(n+1:n+N);
    
    %window
    wx = w.*x;
    
    %fft
    X = fft(wx);
    
    %drop conjugate symmetric part
    %save only unique part
    X1 = X(1:N/2+1);
        
    %delete X (a decoder would not have X, only X1)
    clear X
    
    %plot spectrum
    if (plt > 0)
        XdB = 20*log10(abs(X1));
        XdB = XdB - max(XdB);
        fseq = [0:N/2]'/(N/2)*fs/2;
        plot(fseq, XdB);
        grid
        xlabel('Frequency');
        ylabel('Magnitude, dB');
        title(['Spectrum, block ', num2str(bno)]);
        reply = input('CR for more plots, else 0: ');
        %fprintf('%d\n', reply);
        if ~isempty(reply)
            plt = reply;
        end
    end
    
    %
    %code spectrum
    %
    
    %re-construct conjugate symmetric part
    Xcs = flipud(conj(X1(2:N/2)));
    X = [X1; Xcs];
    
    %inverse fft
    y = ifft(X);
    y = real(y); %just in case ifft has some errors
    
    %window
    wy = w.*y;
    
    %overlap-add
    if (N1 == N2)
        z = yprev+wy(1:N1); %50% overlap
    else
        z = [(yprev + wy(1:N1)); wy(N1+1:N-N1)];
    end
    %save second portion of this block to use next block
    yprev = wy(N-N1+1:N);

    %save output of analysis/synthesis system
    s_out(n+1:n+N2) = z;
    
    %advance input index by block length
    n = n+N2;
end

%clip off N1 at start and fractional block at end
s_out = s_out(N1+1:N1+slen);

%write output file
audiowrite(ofile, s_out, fs);

%calcualate SNR
len = nblocks*N2;
e = s - s_out;
spow = s'*s;
epow = e'*e;
snr = 10*log10(spow/epow);
fprintf('SNR is %f dB for total file\n', snr);

%plot original and reconstructed waveform
ns = [1:slen]';
plot(ns, s, 'b', ns, s_out, 'r');
grid
title('Analysis/Synthesis for FFT and 50% overlap')
ylabel('Amplitude')
xlabel('Sample')
grid
   