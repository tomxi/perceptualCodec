function mark_extract(mfile)
% function mark_extract(mfile)

DEBUG=1;

% load chip
load('chip.mat')
chip_len = length(chip);

[x, fs] = audioread(mfile);
[nsamp, nchn] = size(x);
if (nchn > 1)
    x = x(:,1);
end

% synchronize
% cross-correlation of chip and signal
% look across two chip lengths
% but take first significant peak
% peak indicates synchronization point
% take max(abs())since chip boundardy could be for + or - symbol
xc = xcorr(chip, x(1:2*chip_len) );
[v,i] = max(abs(xc));
% is it first peak?
if (i > chip_len)
    if (max(abs(xc(i-chip_len))) > 0.80*v)
        % take earlier peak
        i = i-chip_len;
    end
end
% n is chip start point
n = i - chip_len;
if (DEBUG)
    plot(xc);
    fprintf('Chip start at sample %d\n', n);
    pause
end

% search for synch word (8 'one' bits)
synch_len = 8;
num_ones = 0;
for i = n:chip_len:nsamp
    xc = xcorr( chip, x(i+1:i+chip_len) );
    % index of peak
    [v, idx] = max(abs(xc));
    if (DEBUG)
        plot(xc); grid
        fprintf('%d %f\n', idx, xc(idx));
        pause
    end
    % sign of peak
    if (xc(idx) > 0)
        num_ones = num_ones+1;
        if num_ones == synch_len;
            %found synch
            if (DEBUG)
                fprintf('Found Message Synch\n');
            end
            n = i+chip_len; %next bit position
            break;
        end
    else
        %not "1", so reset synch bit counter
        num_ones = 0;
    end
end
if (num_ones ~= synch_len)
    fprintf('Message Synch not found\n');
    return;
end

% extract mark bits
% continue if there is one more character of audio data
while (n+(8*chip_len) < nsamp)
    c = 0;
    for j = 7:-1:0 %bits per character
        xc = xcorr( chip, x(n+1:n+chip_len) );
        % index of peak
        [v, idx] = max(abs(xc));
        % sign of peak
        if (xc(idx) > 0)
            c = c + 2^j; %one
        end
        if (DEBUG)
            subplot(2,1,1); plot(x(n+1:n+chip_len)); grid
            subplot(2,1,2); plot(xc); grid
            fprintf('%3d %d %6.2f\n', c, idx, xc(idx));
            pause
        end
        n = n+chip_len;
    end
    if (c == 255)
        c = 10; %ASCII New Line
    end
    if (DEBUG)
        fprintf('%s\n', c);
    else
        fprintf('%s', c);
    end
end

fprintf('\n');
