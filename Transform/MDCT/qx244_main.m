global plt quantize verboseBlock
global fs N SMR X
global Zmax bands_per_bark Nbands N2
global band_binWidth spread_ftn ThrQuietBdB band2bin
global bark_group

global mThrHzM

bark_group = [2,4,6,8,10,14,18,19,20,21,22,23,24,25];
SMR = 20;
plt = 0;
quantize = 1;
verboseBlock = 0;

mdct_ola('vega_orig.wav');
