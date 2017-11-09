function S = Layer1_code(S, scf, SMR, block_bits)

NMR = SMR; %in dB

header_bits = 32;
alloc_bits = 4*32;
sfac_bits = zeros(32);
samp_bits = 0;
allocation = zeros(1,32);
while (1)
    if (block_bits - header_bits - alloc_bits - sum(sfac_bits) - samp_bits < 12)
        break;
    end
    %find largest NMR
    [NMRmax, k] = max(NMR);
    allocation(k) = allocation(k) + 1; %add a bit
    sfac_bits(k) = 6;
    samp_bits = samp_bits + 12; %one bit per sample in this subband
    NMR(k) = NMR(k) - 6; %6 dB per bit
end

for k = 1:32
    %quantize
    if (allocation(k) == 0)
        %inverse quantize
        S(:,k) = zeros(12,1);
    else
        %quantize
        qmax = (2^allocation(k) - 1)/2;
        S_index = round( qmax*S(:,k)/scf(k) );
        %inverse quantize
        qmax = (2^allocation(k) - 1)/2;
        S(:,k) = S_index/qmax*scf(k);
    end
end
    