function [A, a] = simulateBeamforming(mic_array, source_positions, source_type, ...
    f_source, image_freqs, f_samp, c_sound, M_wind, scan_plane_dims, ...
    z_scandistance, N_gridpoints, dynamic_range)

n_parameters = length(mic_array);
n_mics = n_parameters/2;
xi = source_positions;
fs = f_samp;

x(:, 1) = mic_array(1:n_mics);
x(:, 2) = mic_array(n_mics+1:end);
x(:, 3) = 0;

x_min_scanplane = scan_plane_dims(1);
x_max_scanplane = scan_plane_dims(2);
y_min_scanplane = scan_plane_dims(3);
y_max_scanplane = scan_plane_dims(4);

xi_scan = zeros(N_gridpoints, N_gridpoints, 3);
xi_scan(:, :, 1) = ones(N_gridpoints, 1) * linspace(x_min_scanplane, x_max_scanplane, N_gridpoints);
xi_scan(:, :, 2) = linspace(y_min_scanplane, y_max_scanplane, N_gridpoints)' * ones(1, N_gridpoints);
xi_scan(:, :, 3) = ones(N_gridpoints, N_gridpoints) * z_scandistance;

%% Signal generation by monopole source at 100 dB
beta_sq = 1 - dot(M_wind, M_wind);
ti = -4e-3:1/fs:4e-3-1/fs;
if strcmp(source_type, 'gauss');
    p = zeros(n_mics, 32e-3*fs);
    yi = 2e3*sqrt(2)*gauspuls(ti);
    ti = 0:1/fs:32e-3-1/fs;
    
    for I = 1:size(xi,1)
        for J = 1:n_mics
            
            r = sqrt( dot(M_wind, x(J,:) - xi(I,:))^2 + beta_sq*dot(x(J,:)-xi(I,:), x(J,:)-xi(I,:)) );
            ph = ( -dot(M_wind, x(J,:) - xi(I,:)) + r ) / (c_sound * beta_sq);
            
%             r = sqrt(sum((x(J, :) - xi(I, :)).^2));
%             ph = r/c_sound;

            delayed_first_sample = round(ph*fs);
            p(J, delayed_first_sample:delayed_first_sample+length(yi)-1) = ...
                p(J, delayed_first_sample:delayed_first_sample+length(yi)-1) + yi/(4 * pi * r);
        end
    end
    yi(length(yi)+1:32e-3*fs) = 0;
    
elseif strcmp(source_type, 'sine');
    p = zeros(n_mics, 8e-3*fs);
    ti = ti + 4e-3;
    yi = 2*sqrt(2)*cos(2*pi*f_source*ti);
    
    for I = 1:size(xi,1)
        for J = 1:n_mics
            
            r = sqrt( dot(M_wind, x(J,:) - xi(I,:))^2 + beta_sq*dot(x(J,:)-xi(I,:), x(J,:)-xi(I,:)) );
            ph = ( -dot(M_wind, x(J,:) - xi(I,:)) + r ) / (c_sound * beta_sq);
            
%             r = sqrt(sum((x(J, :) - xi(I, :)).^2));
%             ph = r/c_sound;
            delayed_first_sample = round(ph*fs);
            p(J, :) = p(J, :) + 2*sqrt(2)*cos(2*pi*f_source*(ti-ph))/(4 * pi * r);
        end
    end    
    
else
    error('Wrong source type!');
end

[A,a] = ConventionalBeamforming(p, fs, x, xi_scan, image_freqs, c_sound, M_wind);

%% Figures

mic_xmax = ceil(max(x(:, 1)));
mic_xmin = floor(min(x(:, 1)));
mic_ymax = ceil(max(x(:, 2)));
mic_ymin = floor(min(x(:, 2)));

reso = get(0, 'screensize');
f1 = figure('Visible', 'on', ...
           'Position', [floor(reso(3)/2)-525, floor(reso(4)/2)-250, 500, 500], ...
           'Resize', 'off');
           
MicsPositionFigure = axes;
hold(MicsPositionFigure);
set(MicsPositionFigure, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'Position', [0.125+.05 0.15 .7 .7], ...
    'XTick', linspace(mic_xmin, mic_xmax, 5), ...
    'YTick', linspace(mic_ymin, mic_ymax, 5), ...
    'XLim', [mic_xmin mic_xmax], ...
    'YLim', [mic_ymin mic_ymax], ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'Fontsize'    , 12, ...
    'LineWidth'   , 1, ...
    'FontName'   , 'Helvetica' );

hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');

plot(MicsPositionFigure, x(:,1), x(:,2), 'ko');

set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 14          );

hold off

f2 = figure('Visible', 'on', ...
           'Position', [floor(reso(3)/2)+25, floor(reso(4)/2)-250, 500, 500], ...
           'Resize', 'off');
       
BeamformFigure = axes;
hold(BeamformFigure);
set(BeamformFigure, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'Position', [0.125+.0375 0.15 .7 .7], ...
    'XTick', linspace(x_min_scanplane, x_max_scanplane, 5), ...
    'YTick', linspace(y_min_scanplane, y_max_scanplane, 5), ...
    'XLim', [x_min_scanplane x_max_scanplane], ...
    'YLim', [y_min_scanplane y_max_scanplane], ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'Fontsize'    , 12, ...
    'LineWidth'   , 1, ...
    'FontName'   , 'Helvetica' );

min_spl_value = floor(log10(min(A(:))/2e-5^2)*10);
max_spl_value = ceil(log10(max(A(:))/2e-5^2)*10);
imagesc(linspace(x_min_scanplane, x_max_scanplane, N_gridpoints), ...
        linspace(y_min_scanplane, y_max_scanplane, N_gridpoints), ...
        10*log10(A/(2e-5)^2), 'Parent', BeamformFigure, [max_spl_value-dynamic_range max_spl_value]);
cb = colorbar('peer', BeamformFigure, 'Position', [0.875 0.45 0.03 0.4]);
cTit = title(cb, 'SPL (dB)');
hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');
set([hXLabel, hYLabel, cTit], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel, cTit]  , ...
    'FontSize'   , 12          );

plot(BeamformFigure, xi(:,1),xi(:,2), 'wx', 'MarkerSize', 12, 'LineWidth', 2.5);

hold off

end