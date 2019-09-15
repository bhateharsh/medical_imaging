%% Problem 0.2
% Author: Harsh Bhate
% email: bhate@gatech.edu
%% Clearing all
clear all;
clc;
%% Sub-problem(a)
% Generating rect(x)
L_0 = 4;
A_0 = 1;
sample_time = 0.05;
sample_freq = 1/sample_time
x = -(L_0/2):sample_time:(L_0/2);
nos_sample = length(x);
rect_x = A_0*rectpuls(x,L_0);
% Finding the fourier transform
X_omega = fftshift(fft(rect_x)/nos_sample);
% Finding if FFT is real
is_real = isreal(X_omega);
fprintf("Is FFT Real? %s \n", mat2str(is_real));
% Neglecting the imaginary Component
X_omega_abs = abs(X_omega);
X_omega_phase = angle(X_omega);
% Plotting the rect(x)
subplot(3,2,1);
plot (x, rect_x, 'k');
xlim([-5 5])
ylim([0 1.2])
title('Signal ($rect(x)$)',...
    'Interpreter','latex');
xlabel('$x$', ...
    'Interpreter','latex');
ylabel('$rect(x)$', ...
    'Interpreter','latex');
% Plotting the FFT
subplot(3,2,3);
domain = linspace(-sample_freq/2, sample_freq/2, nos_sample);
plot (domain, X_omega_abs, 'k')
% plot (x, rect_x, 'k');
% xlim([-4 4])
title('Fourier transform (Magnitude) ($X(\omega)$)',...
    'Interpreter','latex');
xlabel('$\omega$', ...
    'Interpreter','latex');
ylabel('$|X(\omega)|$', ...
    'Interpreter','latex');
subplot(3,2,5);
domain = linspace(-sample_freq/2, sample_freq/2, nos_sample);
plot (domain, X_omega_phase, 'k')
title('Fourier transform (Phase) ($X(\omega)$)',...
    'Interpreter','latex');
xlabel('$\omega$', ...
    'Interpreter','latex');
ylabel('$\angle X(\omega)$', ...
    'Interpreter','latex');
%% Shifted rect
x = -(L_0/2)+1:sample_time:(L_0/2)+1;
nos_sample = length(x);
rect_x = A_0*rectpuls(x,L_0);
% Finding the fourier transform
X_omega = fftshift(fft(rect_x)/nos_sample);
% Neglecting the imaginary Component
X_omega_abs = abs(X_omega);
X_omega_phase = angle(X_omega);
% Plotting the rect(x)
subplot(3,2,2);
plot (x, rect_x, 'k');
xlim([-5 5])
ylim([0 1.2])
title('Shifted Signal ($rect(x)$)',...
    'Interpreter','latex');
xlabel('$x$', ...
    'Interpreter','latex');
ylabel('$rect(x)$', ...
    'Interpreter','latex');

% Plotting the FFT
subplot(3,2,4);
domain = linspace(-sample_freq/2, sample_freq/2, nos_sample);
plot (domain, X_omega_abs, 'k')
% plot (x, rect_x, 'k');
% xlim([-4 4])
title('Fourier transform (Magnitude) ($X(\omega)$)',...
    'Interpreter','latex');
xlabel('$\omega$', ...
    'Interpreter','latex');
ylabel('$|X(\omega)|$', ...
    'Interpreter','latex');
subplot(3,2,6);
domain = linspace(-sample_freq/2, sample_freq/2, nos_sample);
plot (domain, X_omega_phase, 'k')
title('Fourier transform (Phase) ($X(\omega)$)',...
    'Interpreter','latex');
xlabel('$\omega$', ...
    'Interpreter','latex');
ylabel('$\angle X(\omega)$', ...
    'Interpreter','latex');


