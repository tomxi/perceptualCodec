% script mark_all

mark_embed('mark.wav', 44100, 2500000, 'Hello world ');

%mark_MPEG1('all_441.wav', 'mark.wav', 'foo.wav');
mark_fft('svega_orig_441.wav', 'mark.wav', 'foo.wav');

mark_extract('foo.wav');