function [PS,f]=Powspek(x,SR)

Fs=SR;
T =1/Fs;
L =length(x);

NFFT = 2^nextpow2(L); 
Y=fft(x,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2);
PS=2*abs(Y(1:NFFT/2));
%figure;plot(f,PS);
