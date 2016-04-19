close all; 
clear all; 

% This script is an attempt at implementing Otis Lamont Frost III's 1972
% algorithm for adaptive LCMV beamforming. 

% Constants
J = 4; % Number of delay taps
K = 4; % Number of sensors
F_curly = [1, -2, 1.5, 2]'; % The look direction filter
mu = 0.074; % Step size
Fs = 1; % Sample rate
tau = 1; % The time taken for a sound wave to move directly from one sensor to the next

% Signals
sig_length = 10000; % Note that this is before sampling for observation
s1 = 0.4 * randn(sig_length,1); % Look direction signal
bpFilt = designfilt('bandpassfir', 'FilterOrder', 100, ...
             'CutoffFrequency1', 0.22, 'CutoffFrequency2', 0.28,...
             'SampleRate', 1);
s1 = filter(bpFilt,s1);
s2 = 0.9 * randn(sig_length,1); % Interferer #1
bpFilt = designfilt('bandpassfir', 'FilterOrder', 100, ...
             'CutoffFrequency1', 0.1, 'CutoffFrequency2', 0.12,...
             'SampleRate', 1);
s2 = filter(bpFilt,s2);
s3 = 0.5 * randn(sig_length,1); % Interferer #2
bpFilt = designfilt('bandpassfir', 'FilterOrder', 100, ...
             'CutoffFrequency1', 0.35, 'CutoffFrequency2', 0.4,...
             'SampleRate', 1);
s3 = filter(bpFilt,s3);

% Create observations with delays and noise. Note that the sources are
% currently circularly shifted, is that a problem? 
noise_weight = 0.02; 
tau_1 = 0; % tau_x defines the delay in samples between two neighbouring sensors for signals 1, 2, 3. 
tau_2 = 5; 
tau_3 = 7; 
x1 = s1 + s2 + s3 + noise_weight * randn(length(s1), 1); 
x2 = s1 + s2([end-tau_2+1:end, 1:end-tau_2]) + s3([end-tau_3+1:end, 1:end-tau_3]) + noise_weight * randn(length(s1), 1); 
x3 = s1 + s2([end-2*tau_2+1:end, 1:end-2*tau_2]) + s3([end-2*tau_3+1:end, 1:end-2*tau_3]) + noise_weight * randn(length(s1), 1); 
x4 = s1 + s2([end-3*tau_2+1:end, 1:end-3*tau_2]) + s3([end-3*tau_3+1:end, 1:end-3*tau_3]) + noise_weight * randn(length(s1), 1); 

% Algorithm
initialize W as F
sample the input including four delay taps



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot fft of the three source signals
NFFT = length(s1);
f_s1 = fft(s1,NFFT);
F = [0 : 0.5/(NFFT/2) : 0.5];
magnitude_f_s1 = abs(f_s1);        % Magnitude of the FFT
phase_f_s1 = unwrap(angle(f_s1));  % Phase of the FFT
figure; 
plot(F, magnitude_f_s1(NFFT/2:end))
hold on;
f_s2 = fft(s2,NFFT);
magnitude_f_s2 = abs(f_s2);        % Magnitude of the FFT
phase_f_s2 = unwrap(angle(f_s2));  % Phase of the FFT
plot(F, magnitude_f_s2(NFFT/2:end)); 
hold on;
f_s3 = fft(s3,NFFT);
magnitude_f_s3 = abs(f_s3);        % Magnitude of the FFT
phase_f_s3 = unwrap(angle(f_s3));  % Phase of the FFT
plot(F, magnitude_f_s3(NFFT/2:end))
legend('s1', 's2', 's3');
title('Source power vs frequency');
xlabel('Frequency');ylabel('Power');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%