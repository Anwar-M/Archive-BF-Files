clear all; clc;
source_info = [0.2 0 1 1000 1; 0.1 0.15 1 1000 1];
mic_info = [0 0 0;-.2 -.1 0;.15 -.35 0];

[p, fs] = simulate_arraydata(source_info, mic_info, 343);
N_samples = size(p, 2);

% Very important part, one side fft => mult by 2 and scale to samples
P = 2*fft(p.')/N_samples;
x_fr = fs / N_samples * (0:floor(N_samples/2)-1);
