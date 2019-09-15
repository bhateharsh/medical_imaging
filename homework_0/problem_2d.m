%% Problem 0.2
% Author: Harsh Bhate
% email: bhate@gatech.edu
%% Clearing all
clear all;
clc;
%% Subproblem (d)
x = linspace(-10,10);
y = sinc(x);
% Plotting Function
subplot(3,1,1);
plot(x,y, 'k');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('$y(x) = sinc(x)$', 'interpreter', 'latex');
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');
% Finding the FFT
Y = fft(y);
Y_w = fftshift(Y);
Ys = abs(Y_w);
Ya = angle(Y_w);
% Plot FFT Magnitude
subplot(3,1,2);
plot(Ys,'k');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Magnitude $Y(\omega)$', 'interpreter', 'latex');
xlabel('$\omega$', 'interpreter', 'latex');
ylabel('$|Y(\omega)|$', 'interpreter', 'latex');
% Plot FFT Angle
subplot(3,1,3);
plot(Ya, 'k');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Phase $Y(\omega)$', 'interpreter', 'latex');
xlabel('$\omega$', 'interpreter', 'latex');
ylabel('$\angle Y(\omega)$', 'interpreter', 'latex');