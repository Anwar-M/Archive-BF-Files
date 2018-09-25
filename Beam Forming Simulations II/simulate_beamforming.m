function simulate_beamforming

% h5in_file='Benchmark0.h5';
% mic_info = h5read(h5in_file,'/MetaData/ArrayAttributes/microphonePositionsM').';
load('optim_cfg.mat');
mic_info = mic_info';
% mic_info = rand(15,2); mic_info(:,3) = 0;

% load('optim_cfgf4000.mat');
% mic_info(:,1) = F1(1:32,2);
% mic_info(:,2) = F1(33:64,2);
% mic_info(:,3) = 0;

selection = 1;
save_img = 0;

c = 343.2;
bf_freq = 2000;
N_grid1D = 100;
x_range = 1*[-1 1];
y_range = 1*[-1 1];
z_range = 1;

r_ref = 1;

% load('.\proj\All_Data.mat');
% load('.\proj\Xpos.mat');
% p = data.file_full.';
% fs = 50e3;
% mic_info = [config.x; config.y; config.z].';
% source_info = [0 0 8.1 bf_freq 100];

if selection
    source_info = select_sources(x_range, y_range);
    source_info(:, 3) = z_range;
    source_info(:, 4) = bf_freq;
    source_info(:, 5) = sqrt(2)*2;
else
    source_info = [0 0 z_range bf_freq sqrt(2)*2];
%     source_info = [.75 1 5 3e3 sqrt(2)*2;
%                    -1.2 -.6 5 3e3 sqrt(2)*2;];
end
           
[p, fs] = simulate_arraydata(source_info, mic_info, c);
N_samples = size(p, 2);

% Very important part, one side fft => mult by 2 and scale to samples
P = 2*fft(p.')/N_samples;
x_fr = fs / N_samples * (0:floor(N_samples/2)-1);

indi = round(bf_freq*N_samples/fs+1);
CSM = 0.5*P(indi,:)'*P(indi,:);

tic;
map_beamforming(N_grid1D, x_range, y_range, z_range, r_ref, ...
                CSM, mic_info, bf_freq, c, source_info, save_img);
disp(toc);

% % Test
% X = source_info(:, [1 2 3 5]);
% n_src = size(X, 1);
% n_mic = size(mic_info, 1);
% A = zeros(n_mic, n_src);
% for i = 1:n_src
%     RR = sqrt((mic_info(:,1)-X(i,1)).^2 + (mic_info(:,2)-X(i,2)).^2 + (mic_info(:,3)-X(i,3)).^2);
%     A(:, i) = exp(-2*pi*1i*bf_freq.*RR/c)./(4*pi*RR);
% end
% D = diag(X(:,4),0).^2;
% C_mod = 0.5*conj(A)*D*A.';
% % C_mod = 0.5*A*D*A';
% map_beamforming(N_grid1D, x_range, y_range, z_range, r_ref, ...
%                 C_mod, mic_info, bf_freq, c, source_info, 0);

end

function A = CB(x, CSM, mic_positions, freq, speed_of_sound)

n_mics = length(mic_positions);

r2 = mic_positions.' - x*ones(1, n_mics);

R_2 = sqrt( dot(r2,r2,1) );
xi_timedelay = R_2/speed_of_sound; 
g = (-exp(-2.*pi.*1i.*freq.*xi_timedelay)) ./ ( 4*pi*R_2 );

A = conj(g)*(CSM.')*g.'/(dot(g,g)^2);

end

function map_beamforming(N_grid, x_range, y_range, z_range, r_ref, CSM, ...
                         mic_positions, freq, speed_of_sound, source_info, ...
                         save_img)

xi_scan = zeros(N_grid, N_grid, 3);
A = zeros(N_grid, N_grid);
xi_scan(:, :, 1) = ones(N_grid, 1) * linspace(x_range(1), x_range(2), N_grid);
xi_scan(:, :, 2) = linspace(y_range(1), y_range(2), N_grid)' * ones(1, N_grid);
xi_scan(:, :, 3) = ones(N_grid, N_grid) * z_range;

for k1 = 1:N_grid
    for k2 = 1:N_grid
        A(k1, k2) = CB(squeeze(xi_scan(k1, k2, :)), CSM, mic_positions, freq, speed_of_sound);
%         A(k1, k2) = SarradjBeamforming(squeeze(xi_scan(k1, k2, :)).', CSM, mic_positions.', freq, speed_of_sound, 'form3');
    end
end

A = real(A);

dynamic_range = 12; %dB

set(0,'defaulttextinterpreter','latex');
reso = get(0, 'screensize');
f = figure('Visible', 'on', ...
           'Position', [floor(reso(3)/2)-500, floor(reso(4)/2)-350, 700, 600], ...
           'Resize', 'off');
       
BeamformFigure = axes;
hold(BeamformFigure);
set(BeamformFigure, ...
    'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
    'Position', [0.2 0.15 .7*600/700 .7], ...
    'XTick', linspace(x_range(1), x_range(2), 5), ...
    'YTick', linspace(y_range(1), y_range(2), 5), ...
    'Fontsize'    , 14, ...
    'XLim', [x_range(1) x_range(2)], ...
    'YLim', [y_range(1) y_range(2)]);

max_spl_value = ceil(log10(sqrt(max(A(:)))/(4*pi*r_ref*2e-5))*20);
imagesc(linspace(x_range(1), x_range(2), N_grid), ...
        linspace(y_range(1), y_range(2), N_grid), ...
        2*10*log10(sqrt(A)/(4*pi*r_ref*2e-5)), 'Parent', BeamformFigure, [max_spl_value-dynamic_range max_spl_value]);
    
% max_spl_value = ceil(log10(sqrt(max(A(:)))/(2e-5))*20);
% imagesc(linspace(x_range(1), x_range(2), N_grid), ...
%         linspace(y_range(1), y_range(2), N_grid), ...
%         2*10*log10(sqrt(A)/(2e-5)), 'Parent', BeamformFigure, [max_spl_value-dynamic_range max_spl_value]);
cb = colorbar('peer', BeamformFigure, 'Position', [0.85 0.45 0.03 0.4], 'FontSize', 12);
title(cb, 'SPL (dB)');

% plot(BeamformFigure, source_info(:,1), source_info(:,2), 'wx', 'MarkerSize', 12, 'LineWidth', 1.5);
hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');
hTitle = title(['$f =$ ' num2str(source_info(1,4)) ' Hz']);
set([hXLabel, hYLabel, hTitle]  , ...
    'FontSize'   , 14          );

% viscircles([0,0], 1.1*z_range*0.1052,'EdgeColor','w','LineStyle','-','LineWidth',1.5);
% viscircles([0,0], z_range*0.5774,'EdgeColor','w','LineStyle','-','LineWidth',1.5);

hold off

if save_img
    set(gcf, 'PaperPosition', [0 0 12 10]);
    set(gcf, 'PaperSize', [12 10]);
    print('-dpng', '-r600', ['source_figs\' datestr(now,'yyyymmdd-HHMMSS')]);
    print('-dpng', '-r600', ['source_figs\' datestr(now,'yyyymmdd-HHMMSS')]);
%     print('-dpdf', ['source_figs\' datestr(now,'yyyymmdd-HHMMSS')]);
end

end

% function for manual selection of source position with the mouse
function src = select_sources(x, y)
    reso = get(0,'screensize');
    f2 = figure('Visible','on', ...
               'Position',[floor(reso(3)/2)-250, floor(reso(4)/2)-250, 500, 500], ...
               'Resize','off');

    ax = axes('position',[0.05 0.05 0.9 0.9]);
    box on; grid on;
    axis([x(1) x(2) y(1) y(2)]);
    [src(:,1), src(:,2)] = getpts(ax);
%     x(:,1) = xx; x(:,2) = yy;
    hold on
    for i = 1:size(src, 1)
        plot(src(i,1),src(i,2),'r*');
    end
    hold off
    close(f2);
end
