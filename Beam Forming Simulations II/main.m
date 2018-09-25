clear all; clc;
save_img = 0;
M = [0.4, 0, 0]; % Mach values for wind if used

%% Either one
tic;
load('SimDataDNWnb.mat');
mic_info = data.mic_config.';
c = data.c;
fs = data.sample_frequency;
frequencies = [7000];
N_samples = size(data.file_full, 1);
P = 2*fft(data.file_full)/N_samples;
x_fr = fs / N_samples * (0:floor(N_samples/2)-1);
indi = round(frequencies*N_samples/fs+1);
CSM(1,:,:) = 0.5*P(indi,:)'*P(indi,:);
CSM = squeeze(CSM(1,:,:));
disp(toc);

%% or two
% CsmReader;
% c = 343.3;
% CSM = cpreal + 1i*cpimag;

%% to from the beam
xmin = -1; xmax = 1;
ymin = -2; ymax = 0;
% Define scan grid
n_grid = 100;
z = 0.075;
x = linspace(xmin, xmax, n_grid);
y = linspace(ymin, ymax, n_grid);

A = zeros(n_grid, n_grid, numel(frequencies));
for i = 1:n_grid
    for j = 1:n_grid
        A(i, j, :) = ConventionalBeamforming([x(j), y(i), z], CSM, mic_info, frequencies, c, M);
%         A(i, j, :) = SarradjBeamforming([x(j), y(i), z], CSM, mic_info, frequencies, c, 'form1');
    end
end
%%
for i = 1:numel(frequencies)
%     SPL = 20*log10( sqrt(real(A(:,:,i))) / (4*pi*8.5*2e-5) );
    SPL = 20*log10( sqrt(0.5*real(A(:,:,i))) / 2e-5 );
    maxval = max(SPL(:));
    figure;
    imagesc(x, y, SPL, [round(maxval)-2*12 round(maxval)]); 
    title([num2str(frequencies(i)) ' Hz'], 'Fontsize', 14, 'FontName', 'Helvetica');
    hXLabel = xlabel('x [m]');
    hYLabel = ylabel('y [m]');
    set([hXLabel, hYLabel], ...
        'FontName', 'AvantGarde', 'FontSize', 14);
    
    set(gca, ...
    'YDir','normal', ...
    'XTick', linspace(xmin, xmax, 5), ...
    'YTick', linspace(ymin, ymax, 5), ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'Fontsize'    , 16, ...
    'FontName'   , 'Helvetica' );

    axis equal; axis([xmin xmax ymin ymax]);
    cb = colorbar;
    title(cb, 'SPL (dB)', 'Fontsize', 12, 'FontName', 'Helvetica');
    set(cb, 'Fontsize', 16, 'FontName', 'Helvetica');
    
%     set(gca,'units','centimeters')
%     pos = get(gca,'Position');
%     ti = get(gca,'TightInset');
% 
%     set(gcf, 'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3)+2.5 pos(4)+ti(2)+ti(4)]);
% %     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3)+2.5 pos(4)+ti(2)+ti(4)]);
    
    set(gcf, 'PaperPosition', [0 0 12 10]);
    set(gcf, 'PaperSize', [12 10]);
    if save_img
        print('-dpng', '-r600', ['CB_' num2str(frequencies(i)) '_Hz']);
        close;
    end
end