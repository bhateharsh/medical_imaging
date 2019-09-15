%% Problem 0.6
% Author: Harsh Bhate
% email: bhate@gatech.edu
%% Clearing all
clear all;
clc;
%% Subproblem (1)
colorImg = imread('testImage.png');
colorImg = imresize(colorImg, [400 400]);
img = rgb2gray(colorImg);
% Plotting the Image
subplot(2,4,1);
imshow(img)
title('(a) Original Image',...
    'Interpreter','latex');
%% Subproblem (2)
fftImg = fft2(img);
fftImgCenter = fftshift(fftImg);
F = abs(fftImgCenter);
F = 20*log(F+1);
F = mat2gray(F);
% Plotting the FFT
subplot(2,4,5);
imshow(F);
title('(d) FFT of Original Image',...
    'Interpreter','latex');
%% Subproblem (3)
% Creating Mask
radius = 14;
origin = size(img,2)/2;
[xgrid, ygrid] = meshgrid(1:size(img,2), 1:size(img,1));
mask = ((xgrid-origin).^2 + (ygrid-origin).^2) <= radius.^2;
% Creating Low pass filter
lpf = zeros(size(img));
lpf(mask) = 1;
% Creating Mask
radius = 2;
origin = size(img,2)/2;
[xgrid, ygrid] = meshgrid(1:size(img,2), 1:size(img,1));
mask = ((xgrid-origin).^2 + (ygrid-origin).^2) <= radius.^2;
% Creating the High Pass Filter
hpf = ones(size(img));
hpf(mask) = 0;
% Applying low pass filter
lpfImg = times(lpf, fftImgCenter);
LPF = abs(lpfImg);
LPF = 20*log(LPF+1);
LPF = mat2gray(LPF);
% Applying high pass filter
hpfImg = times(hpf, fftImgCenter);
HPF = abs(hpfImg);
HPF = 20*log(HPF+1);
HPF = mat2gray(HPF);
% Plotting the LPF
subplot (2,4,6);
imshow(LPF);
title('(e) Low Pass Filter',...
    'Interpreter','latex');
% Plotting the HPF
subplot (2,4,7);
imshow(HPF);
title('(f) High Pass Filter',...
    'Interpreter','latex');
%% Subproblem (4)
% Reconstructing LPF
ilpfImg = ifftshift(lpfImg);
reconsLPFImg = ifft2(ilpfImg);
reLPFImg = uint8(real(reconsLPFImg));
% Plotting the Reconstructed LPF
subplot (2,4,2);
imshow(reLPFImg);
title('(b) Reconstructed Low Spatial Frequency Image',...
    'Interpreter','latex');
% Reconstructing HPF
ihpfImg = ifftshift(hpfImg);
reconsHPFImg = ifft2(ihpfImg);
reHPFImg = uint8(real(reconsHPFImg));
% Plotting the Reconstructed LPF
subplot (2,4,3);
imshow(reHPFImg);
title('(c) Reconstructed High Spatial Frequency Image',...
    'Interpreter','latex');
%% Subproblem (5)
newFFT = padarray(fftImgCenter, [400 400], 0);
F = abs(newFFT);
F = 20*log(F+1);
F = mat2gray(F);
% Plotting the new FFT
subplot(2,4,8);
imshow(F);
title('(d2) Zero Padded FFT',...
    'Interpreter','latex');
% Reconstructing Image
newImg = ifftshift(newFFT);
recImg = ifft2(newImg);
reImg = uint8(real(recImg));
% Plotting the Reconstructed Padded FFT
subplot (2,4,4);
imshow(reHPFImg);
title('(a2) Reconstructed Padded Image',...
    'Interpreter','latex');


