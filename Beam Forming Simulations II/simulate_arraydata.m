function [p, fs] = simulate_arraydata(source_info, mic_info, c)

% source_info(l, :) is [x y z freq p_eff] for each row
% fs > 2*f_max, f_max = max( mic_info(:,4) )

N_source = size(source_info, 1);
N_mic = size(mic_info, 1);
f_min = min(source_info(:, 4));
f_max = max(source_info(:, 4));

% sample freq to be 20 times the highest frequency
fs = 20*f_max;
% duration of signal to be 5 periods of lowest frequency
t_end = 5/f_min;
% compute sample points
t = 0:1/fs:(t_end-1/fs);
N_samples = length(t);

p = zeros(N_mic, N_samples);

for I = 1:N_source
    for J = 1:N_mic

        r = sqrt( sum((mic_info(J, :) - source_info(I, 1:3)).^2) );
        ph = r/c;
        p(J, :) = p(J, :) + source_info(I, 5)*cos(2*pi*source_info(I, 4)*(t-ph))/(4 * pi * r);
    end
end    

end