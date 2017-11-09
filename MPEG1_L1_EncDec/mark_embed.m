function mark_embed(ofile, fs, file_len, message_str)
% function mark_embed(ofile, fs, file_len, message_str)
%
% ofile is *.wav output file. This is contains pseudo-random 
%   "chip" sequences that represent the message
% fs is sampling rate of ofile, in Hz
% file_len is number of samples in ofile
% message_str is quoted ASCII string that is embedded
%
% This is the pseudo-random sequence that is to be added 
% to an audio signal such that it is below the perceptual threshold
%

DEBUG = 0;

chip_len = 16*1024;
chip_len = 4*1024;

% same chip is always used
% chip is uniform [-1 1]
chip = 2*(rand(chip_len, 1)-0.5); %column vector

% message is quoted string
% convert into bits
% 8 bits per ASCII character
% MSB is always zero
n = length(message_str);
msg_len = n*8;

% Message Synch is 8 ones, which is never an ASCII character
synch_len = 8; 
message_bits = zeros(synch_len+msg_len,1); %allocate array
message_bits(1:synch_len) = ones(synch_len,1); %synch bits
k = synch_len+1;
for i = 1:n
    c = uint8(message_str(i));
    for j=7:-1:0 %message bits are msb (or leftmost) first
        v = bitand(c, 2^j);
        if (v == 0) 
            message_bits(k) = -1; %symbol for zero
        else
            message_bits(k) = 1; %symbol for one
        end
        k = k+1;
    end    
    fprintf('%s %3d %3k: %2d %2d %2d %2d %2d %2d %2d\n', c, uint8(c),...
        k-8, message_bits(k-8:k-1) );
end
fprintf('%d message bits\n', length(message_bits));

num_chips = ceil(file_len/chip_len);
y = zeros(num_chips*chip_len, 1);
n = 0;
% create mark signal by looping through message bits
% this will put message synch character prior to every coding of 
% message
for i = 0:num_chips-1
    y(n+1:n+chip_len) = message_bits(mod(i,msg_len) + 1)*chip;
    n = n+chip_len;
end

% write mark
audiowrite(ofile, y, fs);

% save chip
save('chip.mat', 'chip');


    



