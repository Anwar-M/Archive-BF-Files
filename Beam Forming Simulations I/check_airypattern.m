% load('spiralOptimized.mat');
clc;
clear all; close all;
random = 0;
save_image = 1;
add_bounds = 1;
c = 343.3;
f = 2*1000;
xi = [0 0 35];
n = 32;
TH = rand(n,1)*2*pi;
% R = ones(n,1);
% R = rand(n,1);
R = [ones(n/2,1); rand(n/2,1)];
[x(:,1), x(:,2)] = pol2cart(TH,sqrt(R));

%%
cfg = 1;
set(0,'defaulttextinterpreter','latex');
close all;
% load('spiralOptimized.mat');
clear x;
if cfg == 1;
    load('int_config.mat');
elseif cfg == 2;
    load('int_random.mat');
elseif cfg == 3;
    load('spiralOptimized.mat');
else
    error('error is this');
end
% clear x; load('best_config.mat'); x(:,1) = A1(1:32); x(:,2) = A1(33:end);

Energy = 0;
if random
    while Energy < 42
        x = 2*rand(32,2)-1;
        a1 = 3.83;
        a2 = pi*f/c;
        Energy = ObjectFunction([x(:,1); x(:,2)], size(x,1), a1, a2);
    end
else
    a1 = 3.83;
    a2 = pi*f/c;
    Energy = ObjectFunction([x(:,1); x(:,2)], size(x,1), a1, a2);
end

dis1 = xi(3)*tan(a1*c/(2*pi*f));
dis2 = xi(3)*tan(asin(a2*c/(2*pi*f)));

n_samp = 250;
[A, a] = simulateBeamforming([x(:,1); x(:,2)], xi, 'sine', f, f, 50e3, c, [0,0,0], xi(3).*[-1 1 -1 1], xi(3), n_samp, 12*2);
hTitle = title(gca, ['$J$ = ' num2str(round(Energy*100)/100)]);
[X,Y] = meshgrid(linspace(-xi(3),xi(3),n_samp),linspace(-xi(3),xi(3),n_samp));
set(hTitle  , 'FontSize'   , 14          );

if add_bounds
    viscircles([0,0], dis1, 'EdgeColor','w','LineStyle','--','LineWidth',1);
    viscircles([0,0], dis2, 'EdgeColor','w','LineStyle','--','LineWidth',1);
end

if save_image
%     colormap gray;
    set(gcf, 'PaperPosition', [0 0 11.5 11.5]);
    set(gcf, 'PaperSize', [11.5 11.5]);
    if cfg==1
        print(gcf,'-dpdf','bf_cfg2'); close;
    elseif cfg==2
        print(gcf,'-dpdf','bf_cfg1'); close;
    end
%     colormap gray;
    set(gcf, 'PaperPosition', [0 0 11.5 11.5]);
    set(gcf, 'PaperSize', [11.5 11.5]);
    if cfg==1
        print(gcf,'-dpdf','mic_cfg2');
    elseif cfg==2
        print(gcf,'-dpdf','mic_cfg1');
    end
end

aa = abs(a)./max(abs(a(:)));

% figure; mesh(X, Y, aa);
figure; plot(X(1,:), aa(125,:),'Color','k');
hXLabel = xlabel('$x$ [m]');
% hYLabel = ylabel('p_0(x)/p_0(0)');
hYLabel = ylabel('Normalized $P_0$');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 14          );
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XTick'       , [-xi(3) -round(dis2*100)/100 -round(dis1*100)/100 round(dis1*100)/100 round(dis2*100)/100 xi(3)] , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'Fontsize'    , 12, ...
  'LineWidth'   , 1         );

if add_bounds
    line([dis1;dis1],[0;1],'LineWidth',1,'LineStyle','--','Color','k');
    line([-dis1;-dis1],[0;1],'LineWidth',1,'LineStyle','--','Color','k');
    line([dis2;dis2],[0;1],'LineWidth',1,'LineStyle','--','Color','k');
    line([-dis2;-dis2],[0;1],'LineWidth',1,'LineStyle','--','Color','k');
end
if save_image
    set(gcf, 'PaperPosition', [0 0 13 10.5]);
    set(gcf, 'PaperSize', [13 10.5]);
    if cfg==1
        print(gcf,'-dpdf','inter_cfg2');
    elseif cfg==2
        print(gcf,'-dpdf','inter_cfg1');
    end
end
