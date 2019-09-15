%% Problem 0.6
% Author: Harsh Bhate
% email: bhate@gatech.edu
%% Clearing all
clear all;
clc;
%% Generate cos Signal
sampleTime = 0.005;
sampleFreq = 1/sampleTime;
x = -1:sampleTime:1;
y = -1:sampleTime:1;
[X, Y] = meshgrid(x,y);
u_0 = 4;
v_0 = 2;
s = cos(2*pi*(u_0*X + v_0*Y));
% Plotting the signal
subplot (3,2,1);
mesh(X,Y,s);
colormap(gray);
title('$s(x,y) = \cos[2 \pi (u_0x + v_0y)]$','interpreter', 'latex');
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');
zlabel('$s(x,y)$', 'interpreter', 'latex');
%% Find the FFT of the signal
S = fft2(s);
Sc= fftshift(S);
Sc = abs(Sc);
Sc = 20*log(Sc+1);
Sa = angle(Sc);
% Plotting the FFT Magnitude
subplot (3,2,3);
imshow(Sc);
title('$|S(u,v)|$','Interpreter','latex');
xlabel('u', 'interpreter', 'latex');
ylabel('v', 'interpreter', 'latex');
% Plotting the FFT Angle
subplot (3,2,5);
imshow(Sa);
title('$\angle S(u,v)$','Interpreter','latex');
xlabel('u', 'interpreter', 'latex');
ylabel('v', 'interpreter', 'latex');
%% Generate sin signal
sampleTime = 0.005;
sampleFreq = 1/sampleTime;
x = -1:sampleTime:1;
y = -1:sampleTime:1;
[X, Y] = meshgrid(x,y);
u_0 = 4;
v_0 = 2;
s = sin(2*pi*(u_0*X + v_0*Y));
% Plotting the signal
subplot (3,2,2);
mesh(X,Y,s);
colormap(gray);
title('$s(x,y) = sin[2 \pi (u_0x + v_0y)]$','interpreter', 'latex');
xlabel('x', 'interpreter', 'latex');
ylabel('y', 'interpreter', 'latex');
zlabel('$s(x,y)$', 'interpreter', 'latex');
%% Find the FFT of the signal
S = fft2(s);
Sc= fftshift(S);
Sc = abs(Sc);
Sc = 20*log(Sc+1);
Sa = angle(Sc);
% Plotting the FFT Magnitude
subplot (3,2,4);
imshow(Sc);
title('$|S(u,v)|$','Interpreter','latex');
xlabel('u', 'interpreter', 'latex');
ylabel('v', 'interpreter', 'latex');
% Plotting the FFT Angle
subplot (3,2,6);
imshow(Sa);
title('$\angle S(u,v)$','Interpreter','latex');
xlabel('u', 'interpreter', 'latex');
ylabel('v', 'interpreter', 'latex');
