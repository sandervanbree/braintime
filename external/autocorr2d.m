function A=autocorr2d(I)

% By Tristan Ursell. Modified by Sander van Bree
%
% Compute the 2D spatial autocorrelation of a matrix or image I using the 
% Wiener - Khintchine Theorem. The output is the normalized correlation
% coefficient -1 < C < 1.
%
% The center pixel of A will have C = 1. Using images with odd dimensions 
% will give results that are easier to interpret.
%
% ref: http://mathworld.wolfram.com/Wiener-KhinchinTheorem.html
%
% %EX 1: 2D zero mean Gaussian process with variance = 1 
% %The result is a delta-function at zero, with spurious correlations
% around zero.
% 
% I=randn(501,501);
% A=autocorr2d(I);
% 
% figure;
% subplot(1,2,1)
% imagesc(I)
% axis equal tight
% box on
% xlabel('X')
% ylabel('Y')
% title('Original Image')
% colorbar
% 
% subplot(1,2,2)
% imagesc(-250:1:250,-250:1:250,A,[-1,1])
% axis equal tight
% box on
% xlabel('X')
% ylabel('Y')
% title('Spatial Autocorrelation')
% colorbar

%EX 2: Using the 'rice.png' image
%The central bright spot is the correlation coming from the rice grains,
%the diffuse horizontal band is the correlation coming from the
%background.

% I=imread('rice.png');
% A=autocorr2d(I);
% 
% figure;
% subplot(1,2,1)
% imagesc(I)
% axis equal tight
% box on
% xlabel('X')
% ylabel('Y')
% title('Original Image')
% colorbar
% 
% subplot(1,2,2)

I=double(I); %convert to double
I=I-mean(I(:)); %subtract mean
I=I/sqrt(sum(I(:).^2)); %normalize magnitude
fft_I=fft2(I); %compute fft2
A=real(fftshift(ifft2(fft_I.*conj(fft_I)))); %compute autocorrelation

% d_max=size(I,1)/2 +1;
% imagesc(-size(I,1)/2:size(I,1)/2-1,-size(I,1)/2:size(I,1)/2-1,A,[-1,1])
% axis equal tight
% box on
% xlabel('X')
% ylabel('Y')
% title('Spatial Autocorrelation')
% colorbar