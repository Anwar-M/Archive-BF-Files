clear all;
f = 2000;
c = 340;
D = 2;

delta = 0.001;
dmin = 0.2; dmax = 0.35;
range = linspace(dmin, dmax, 1e2);

for i = 1:length(range)
    x = HexConfig(range(i));
    Dx = max(x(:,1))-min(x(:,1));
    Dy = max(x(:,2))-min(x(:,2));
    Deff = (Dx>=Dy)*Dx + (Dx<Dy)*Dy;
    E(i,1) = ObjectFunction([x(:,1) x(:,2)], 32, 2*3.83/D, pi*f/340);
    E(i,2) = ObjectFunction([x(:,1) x(:,2)], 32, 2*3.83/Deff, pi*f/340);
end

%%

figure; plot(range, E(:,1), 'Color', 'k'); axis([range(1) range(end) 6.2 8])
hXLabel = xlabel('d [m]');
hYLabel = ylabel('J');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 12          );
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'YGrid'       , 'off'      , ...
  'XTick'       , dmin:0.025:dmax, ...
  'YTick'       , 6.4:0.2:8, ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'Fontsize'    , 12, ...
  'LineWidth'   , 1         );

    set(gcf, 'PaperPosition', [0 0 14.0 10.5]);
    set(gcf, 'PaperSize', [14.0 10.5]);
%     print -dpdf bf_minc.pdf
    colormap gray;
    print(gcf,'-dpdf', 'a');