% clear all
close all;
%% Parameters
imaged_frequency = 2000;
xi = 35*[0 0 1;];
source_signal = 'sine';
f_source = imaged_frequency;
fs = 50e3;
c = 343.3;
M = 0*[1 0 0]; % wind

dynamic_range = 12*2; %dB

add_bounds = 1;
save_image = 1;
manip = 0;

f_optim = 2000;
%% Scanning grid
scanplane_x_max = 1.0*35;
scanplane_y_max = 1.0*35;
N_grid = 5*100;
z_scan = xi(1,3);
xi_scan = zeros(N_grid, N_grid, 3);
xi_scan(:, :, 1) = ones(N_grid, 1) * linspace(-scanplane_x_max, scanplane_x_max, N_grid);
xi_scan(:, :, 2) = linspace(-scanplane_y_max, scanplane_y_max, N_grid)' * ones(1, N_grid);
xi_scan(:, :, 3) = ones(N_grid, N_grid) * z_scan;

%% Read in mic positions
x = zeros(32,3);

load('best_config.mat');
x(:, 1) = A1(1:32);
x(:, 2) = A1(33:64);

% load('config_data.mat');
% x(:, 1) = A1(1:32,2);
% x(:, 2) = A1(33:64,2);

x(:, 3) = 0;

% x(:,1) = F1(1:32,2,1); x(:,2) = F1(33:64,2,1);
if manip
    x(33,:) = 1.0*[1 0 0];
    x(34,:) = 1.0*[0 1 0];
    x(35,:) = 1.0*[-1 0 0];
    x(36,:) = 1.0*[0 -1 0];
end
n_mics = size(x,1);

% x(:,1:2) = HexConfig(.2369);

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
     
            delayed_first_sample = round(ph*fs);
            p(J, :) = p(J, :) + 2*sqrt(2)*cos(2*pi*f_source*(ti-ph))/(4 * pi * r);
        end
    end    
    
else
    error('wrong');
end
%%
A = ConventionalBeamforming(p, fs, x, xi_scan, imaged_frequency, c, M);

%% Figures
set(0,'defaulttextinterpreter','latex');

figure; plot(x(:,1), x(:,2), 'ro', 'MarkerFaceColor', 'r'); axis equal; axis([-1 1 -1 1]);
hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');
% hTitle = title(['$J$ = ' num2str(round(ObjectFunction([x(:,1)' x(:,2)'], n_mics,2*3.83/2,pi*2000/c)*100)/100)]);
hTitle = title([]);
set(gca, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'Position', [0.125+2*.025 0.15 .75 .75], ...
    'XTick', linspace(-1, 1, 5), ...
    'YTick', linspace(-1, 1, 5), ...
    'XLim', [-1 1], ...
    'YLim', [-1 1], ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'Fontsize'    , 12, ...
    'LineWidth'   , 1, ...
    'FontName'   , 'Helvetica' );
    set([hXLabel, hYLabel, hTitle], ...
        'FontName'   , 'AvantGarde');
    set([hXLabel, hYLabel, hTitle]  , ...
        'FontSize'   , 14          );

if save_image
    set(gcf, 'PaperPosition', [0 0 10.5 10.5]);
    set(gcf, 'PaperSize', [10.5 10.5]);
    print(gcf,'-dpdf', 'mic_config_manual');
end

reso = get(0, 'screensize');
f = figure('Visible', 'on', ...
           'Position', [floor(reso(3)/2)-500, floor(reso(4)/2)-350, 700, 600], ...
           'Resize', 'off');

BeamformFigure = axes;
hold(BeamformFigure);

hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');
titStrng = ['$f$ = ' num2str(imaged_frequency) ' Hz'];
hTitle = title(titStrng);


set(BeamformFigure, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'Position', [0.2 0.15 .7*600/700 .7], ...
    'XTick', linspace(-scanplane_x_max, scanplane_x_max, 5), ...
    'YTick', linspace(-scanplane_y_max, scanplane_y_max, 5), ...
    'XLim', [-scanplane_x_max scanplane_x_max], ...
    'YLim', [-scanplane_y_max scanplane_y_max]);

min_spl_value = floor(log10(min(A(:))/2e-5^2)*10);
max_spl_value = ceil(log10(max(A(:))/2e-5^2)*10)+1;
imagesc(linspace(-scanplane_x_max, scanplane_x_max, N_grid), ...
        linspace(-scanplane_y_max, scanplane_y_max, N_grid), ...
        10*log10(A/(2e-5)^2), 'Parent', BeamformFigure, [max_spl_value-dynamic_range max_spl_value]);

cb = colorbar('peer', BeamformFigure, 'Position', [0.85 0.45 0.03 0.4], 'FontSize', 12);
cTit = title(cb, 'SPL [dB]');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel, cTit, hTitle], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel, cTit, hTitle]  , ...
    'FontSize'   , 14          );
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'Fontsize'    , 14, ...
  'LineWidth'   , 1         );

if add_bounds
    r_min = xi(3)*tan( 1.22*c/(2*f_optim) );
    r_max = xi(3)*tan( pi/6 );
%     viscircles([0,0], 0.1052+0.0252,'EdgeColor','w','LineStyle','--','LineWidth',1);
%     viscircles([0,0], 0.5774,'EdgeColor','w','LineStyle','--','LineWidth',1);
    viscircles([0,0], r_min+0*r_min/4,'EdgeColor','w','LineStyle','--','LineWidth',1);
    viscircles([0,0], r_max,'EdgeColor','w','LineStyle','--','LineWidth',1);
end

if save_image
    set(gcf, 'PaperPosition', [0 0 12 10]);
    set(gcf, 'PaperSize', [12 10]);
%     print -dpdf bf_minc.pdf
%     colormap gray;
    print(gcf,'-dpdf', 'bf_hexa');
end
