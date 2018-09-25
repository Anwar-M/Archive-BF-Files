function simulate_csm

h5in_file='Benchmark0.h5';
mic_info = h5read(h5in_file,'/MetaData/ArrayAttributes/microphonePositionsM').';

c = 343.3;
bf_freq = 1e3;
N_grid1D = 100;
x_range = [-1 1];
y_range = [-1 1];
z_range = 1;

%%

source_info(:,:,1) = [.75 0 1 1e3 sqrt(2)*2;
                   -.75 0 1 1e3 sqrt(2)*2;];

[p, fs] = simulate_arraydata(source_info(:,:,1), mic_info, c);
N_samples = size(p, 2);

% Very important part, one side fft => mult by 2 and scale to samples
P = 2*fft(p.')/N_samples;
x_fr = fs / N_samples * (0:floor(N_samples/2)-1);

indi = round(bf_freq*N_samples/fs+1);
CSM(1,:,:) = 0.5*P(indi,:)'*P(indi,:);
%%

source_info(:,:,2) = [.1 0 1 1e3 sqrt(2)*2;
                   -.1 0 1 1e3 sqrt(2)*2;];

[p, fs] = simulate_arraydata(source_info(:,:,2), mic_info, c);
N_samples = size(p, 2);

% Very important part, one side fft => mult by 2 and scale to samples
P = 2*fft(p.')/N_samples;
x_fr = fs / N_samples * (0:floor(N_samples/2)-1);

indi = round(bf_freq*N_samples/fs+1);
CSM(2,:,:) = 0.5*P(indi,:)'*P(indi,:);

%%
source_info(:,:,3) = [.15 0 1 1e3 sqrt(2)*2;
                   -.15 0 1 1e3 sqrt(2)*2;];

[p, fs] = simulate_arraydata(source_info(:,:,3), mic_info, c);
N_samples = size(p, 2);

% Very important part, one side fft => mult by 2 and scale to samples
P = 2*fft(p.')/N_samples;
x_fr = fs / N_samples * (0:floor(N_samples/2)-1);

indi = round(bf_freq*N_samples/fs+1);
CSM(3,:,:) = 0.5*P(indi,:)'*P(indi,:);

%%
source_info(:,:,4) = [.2 0 1 1e3 sqrt(2)*2;
                   -.2 0 1 1e3 sqrt(2)*2;];

[p, fs] = simulate_arraydata(source_info(:,:,4), mic_info, c);
N_samples = size(p, 2);

% Very important part, one side fft => mult by 2 and scale to samples
P = 2*fft(p.')/N_samples;
x_fr = fs / N_samples * (0:floor(N_samples/2)-1);

indi = round(bf_freq*N_samples/fs+1);
CSM(4,:,:) = 0.5*P(indi,:)'*P(indi,:);

%%

save('simulated_csm_data.mat', 'CSM', 'c', 'bf_freq', 'mic_info', 'source_info');

% map_beamforming(N_grid1D, x_range, y_range, z_range, ...
%                 CSM, mic_info, bf_freq, c, source_info);

end

function A = CB(x, CSM, mic_positions, freq, speed_of_sound)

n_mics = length(mic_positions);

r2 = mic_positions.' - x*ones(1, n_mics);

R_2 = sqrt( dot(r2,r2,1) );
xi_timedelay = R_2/speed_of_sound; 
g = (-exp(-2.*pi.*1i.*freq.*xi_timedelay)) ./ ( 4*pi*R_2 );

A = conj(g)*(CSM.')*g.'/(dot(g,g)^2);

end

function map_beamforming(N_grid, x_range, y_range, z_range, CSM, ...
                         mic_positions, freq, speed_of_sound, source_info)

xi_scan = zeros(N_grid, N_grid, 3);
A = zeros(N_grid, N_grid);
xi_scan(:, :, 1) = ones(N_grid, 1) * linspace(x_range(1), x_range(2), N_grid);
xi_scan(:, :, 2) = linspace(y_range(1), y_range(2), N_grid)' * ones(1, N_grid);
xi_scan(:, :, 3) = ones(N_grid, N_grid) * z_range;

for k1 = 1:N_grid
    for k2 = 1:N_grid
        A(k1, k2) = CB(squeeze(xi_scan(k1, k2, :)), CSM, mic_positions, freq, speed_of_sound);
    end
end

A = real(A);

dynamic_range = 3; %dB

reso = get(0, 'screensize');
f = figure('Visible', 'on', ...
           'Position', [floor(reso(3)/2)-500, floor(reso(4)/2)-350, 700, 600], ...
           'Resize', 'off');
       
BeamformFigure = axes;
hold(BeamformFigure);
set(BeamformFigure, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'Position', [0.1 0.1 .8*600/700 .8], ...
    'XTick', linspace(x_range(1), x_range(2), 11), ...
    'YTick', linspace(y_range(1), y_range(2), 11), ...
    'XLim', [x_range(1) x_range(2)], ...
    'YLim', [y_range(1) y_range(2)]);

min_spl_value = floor(log10(min(A(:))/2e-5^2)*10);
max_spl_value = ceil(log10(max(A(:))/2e-5^2)*10);
imagesc(linspace(x_range(1), x_range(2), N_grid), ...
        linspace(y_range(1), y_range(2), N_grid), ...
        10*log10(A/(2e-5)^2), 'Parent', BeamformFigure, [max_spl_value-dynamic_range max_spl_value]);
cb = colorbar('peer', BeamformFigure, 'Position', [0.9 0.45 0.03 0.4]);
title(cb, 'SPL (dB)');

plot(BeamformFigure, source_info(:,1), source_info(:,2), 'wx', 'MarkerSize', 12, 'LineWidth', 1.5);

hold off

if 0
    set(gcf, 'PaperPosition', [0 0 12 10]);
    set(gcf, 'PaperSize', [12 10]);
    print('-dpng', '-r600', ['source_figs\' datestr(now,'yyyymmdd-HHMMSS')]);
end

end
