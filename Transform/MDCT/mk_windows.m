function mk_windows()

Nlong = 1024;
num_short = 8;
Nshort = Nlong/num_short;
NL2 = Nlong/2;
NS2 = Nshort/2;

wlong = sin(pi*([1:Nlong]-0.5)'/Nlong); %long window
wshort = sin(pi*([1:Nshort]-0.5)'/Nshort); %long window
wstart = [wlong(1:Nlong/2); ones(4*NS2+NS2/2,1); wshort(Nshort/2+1:Nshort)];
wstop = flipud(wstart);

plot(wlong); title('long'); pause;
plot(wstart); title('start'); pause
plot(wshort); title('short'); pause
plot(wstop); title('stop'); pause
