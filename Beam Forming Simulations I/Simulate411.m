% clear all;
clc;
addpath('BA');

%% Parameters
% imaged_frequency = 1000;
% imaged_frequency = 975:50:1025;
imaged_frequency = 2000;
% imaged_frequency = 1450:10:1500;
% xi = [ 0.5 -0.1 1;
%       -0.2 -0.4 1];
xi = [0 0 1;];
source_signal = 'sine';
f_source = imaged_frequency;
fs = 50e3;
c = 343.2;
M = 0.5*[1 0 0]; % wind

dynamic_range = 12; %dB

select_mics = 66;

print_to_file = 0;

%% Scanning grid
scanplane_x_max = 1.0*1;
scanplane_y_max = 1.0*1;
N_grid = 5*100;
z_scan = xi(1,3);
xi_scan = zeros(N_grid, N_grid, 3);
xi_scan(:, :, 1) = ones(N_grid, 1) * linspace(-scanplane_x_max, scanplane_x_max, N_grid);
xi_scan(:, :, 2) = linspace(-scanplane_y_max, scanplane_y_max, N_grid)' * ones(1, N_grid);
xi_scan(:, :, 3) = ones(N_grid, N_grid) * z_scan;

%% Read in mic positions
mikz = 6;
if select_mics == 1
    reso = get(0,'screensize');
    f2 = figure('Visible','on', ...
               'Position',[floor(reso(3)/2)-250, floor(reso(4)/2)-250, 500, 500], ...
               'Resize','off');

    ax = axes('position',[0.05 0.05 0.9 0.9]);
    box on; grid on;
    axis([-scanplane_x_max ...
           scanplane_x_max ...
          -scanplane_y_max ...
           scanplane_y_max]);
    [x(:,1), x(:,2)] = getpts(ax);
%     x(:,1) = xx; x(:,2) = yy;
    hold on
    for i = 1:size(x, 1)
        plot(x(i,1),x(i,2),'r*');
    end
    hold off
    close(f2);
    
elseif select_mics == 6
    load('config_data.mat');
    x(:,1) = A1(1:32, 2);
    x(:,2) = A1(33:end, 2);
    
elseif select_mics ==  2
    for I = 1:mikz
        x( (I-1)*mikz+1:I*mikz, 1 ) = -scanplane_x_max:2*scanplane_x_max/(mikz-1):scanplane_x_max;
        x( (I-1)*mikz+1:I*mikz, 2 ) = -scanplane_y_max + 2*(I-1)*scanplane_y_max/(mikz-1);
    end
elseif select_mics ==  3
    load('spiralOptimized.mat');
    x = x/2;
elseif select_mics == 4
    clear x;
    x(:,1) = 2*scanplane_x_max*rand(32,1) - scanplane_x_max;
    x(:,2) = 2*scanplane_y_max*rand(32,1) - scanplane_y_max;
elseif select_mics == 5
    clear x;
%     x = x_ga;
    x(:,1) = A1(1:32,4);
    x(:,2) = A1(33:64,4);
else
    x = dlmread('Array.txt');
%     x = z(:,5:6);
end

x(:, 3) = 0;
n_mics = size(x, 1);


%% Signal generation by monopole source at 100 dB
beta_sq = 1 - dot(M, M);

ti = -4e-3:1/fs:4e-3-1/fs;
if strcmp(source_signal, 'gauss');
    p = zeros(n_mics, 32e-3*fs);
    yi = gauspuls(ti);
    ti = 0:1/fs:32e-3-1/fs;
    
    for I = 1:size(xi,1)
        for J = 1:n_mics
            
            r = sqrt( dot(M, x(J,:) - xi(I,:))^2 + beta_sq*dot(x(J,:)-xi(I,:), x(J,:)-xi(I,:)) );
            ph = ( -dot(M, x(J,:) - xi(I,:)) + r ) / (c * beta_sq);
            
%             r = sqrt(sum((x(J, :) - xi(I, :)).^2));
%             ph = r/c;

            delayed_first_sample = round(ph*fs);
            p(J, delayed_first_sample:delayed_first_sample+length(yi)-1) = ...
                p(J, delayed_first_sample:delayed_first_sample+length(yi)-1) + yi/(4 * pi * r);
        end
    end
    yi(length(yi)+1:32e-3*fs) = 0;
    
elseif strcmp(source_signal, 'sine');
    p = zeros(n_mics, 8e-3*fs);
    ti = ti + 4e-3;
    yi = 2*sqrt(2)*cos(2*pi*f_source*ti);
    
    for I = 1:size(xi,1)
        for J = 1:n_mics
            
            r = sqrt( dot(M, x(J,:) - xi(I,:))^2 + beta_sq*dot(x(J,:)-xi(I,:), x(J,:)-xi(I,:)) );
            ph = ( -dot(M, x(J,:) - xi(I,:)) + r ) / (c * beta_sq);
            
%             r = sqrt(sum((x(J, :) - xi(I, :)).^2));
%             ph = r/c;
            delayed_first_sample = round(ph*fs);
            p(J, :) = p(J, :) + 2*sqrt(2)*cos(2*pi*f_source*(ti-ph))/(4 * pi * r);
        end
    end    
    
else
    error('wrong');
end
%%

[A,a] = ConventionalBeamforming(p, fs, x, xi_scan, imaged_frequency, c, M);
% A_CLEAN = BF_CLEANPST2(p', fs, c, linspace(-scanplane_x_max, scanplane_x_max, N_grid), ...
%     linspace(-scanplane_y_max, scanplane_y_max, N_grid), z_scan, imaged_frequency, x(:,1:2));
% A_CLEANSC = BF_CLEANSC(p', fs, c, linspace(-scanplane_x_max, scanplane_x_max, N_grid), ...
%     linspace(-scanplane_y_max, scanplane_y_max, N_grid), z_scan, imaged_frequency, x(:,1:2));

%% Figures
% A = A_CLEANSC;

mic_xmax = ceil(max(x(:, 1)));
mic_xmin = floor(min(x(:, 1)));
mic_ymax = ceil(max(x(:, 2)));
mic_ymin = floor(min(x(:, 2)));

reso = get(0, 'screensize');
f = figure('Visible', 'on', ...
           'Position', [floor(reso(3)/2)-500, floor(reso(4)/2)-250, 1000, 500], ...
           'Resize', 'off');
           
MicsPositionFigure = axes;
hold(MicsPositionFigure);
set(MicsPositionFigure, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'Position', [0.075 0.15 0.35 0.7], ...
    'XTick', linspace(mic_xmin, mic_xmax, 11), ...
    'YTick', linspace(mic_ymin, mic_ymax, 11), ...
    'XLim', [mic_xmin mic_xmax], ...
    'YLim', [mic_ymin mic_ymax]);
plot(MicsPositionFigure, x(:,1), x(:,2), 'r*');
title(['J = ' num2str(ObjectFunction([x(:,1)' x(:,2)'], length(x(:,1)), 3.83, pi*2000/c))]);


hold off

BeamformFigure = axes;
hold(BeamformFigure);
set(BeamformFigure, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'Position', [0.525 0.15 0.35 0.7], ...
    'XTick', linspace(-scanplane_x_max, scanplane_x_max, 11), ...
    'YTick', linspace(-scanplane_y_max, scanplane_y_max, 11), ...
    'XLim', [-scanplane_x_max scanplane_x_max], ...
    'YLim', [-scanplane_y_max scanplane_y_max]);

min_spl_value = floor(log10(min(A(:))/2e-5^2)*10);
max_spl_value = ceil(log10(max(A(:))/2e-5^2)*10);
imagesc(linspace(-scanplane_x_max, scanplane_x_max, N_grid), ...
        linspace(-scanplane_y_max, scanplane_y_max, N_grid), ...
        10*log10(A/(2e-5)^2), 'Parent', BeamformFigure, [max_spl_value-dynamic_range max_spl_value]);
cb = colorbar('peer', BeamformFigure, 'Position', [0.9 0.45 0.03 0.4]);
title(cb, 'SPL (dB)');

plot(BeamformFigure, xi(:,1),xi(:,2), 'wx', 'MarkerSize', 12, 'LineWidth', 2.5);

hold off


if print_to_file
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf, 'PaperPosition', [-6.0 +1 14.0 10.5]);
    set(gcf, 'PaperSize', 1.1*[14.0/2 10.5]);
    set(f,'Units','inches');
    screenposition = get(f,'Position');
    set(f,...
        'PaperPosition',[-screenposition(3)/2 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3)/2 screenposition(4)]);
    print('-dbmp','Figoer.bmp');
%     print -dpdf -painters epsFig
%     figure
%     bla = axes;
%     new_handuru = copyobj(BeamformFigure, bla);
%     plot(bla, xi(:,1),xi(:,2), 'kx', 'MarkerSize', 12, 'LineWidth', 2.5);
%     print(f,'-dpng','Figoer.png');
end