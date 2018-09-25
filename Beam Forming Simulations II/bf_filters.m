function bf_filters

CsmReader;
x_m = x_m.';

scanplane_x_max = 1.0*1;
scanplane_y_max = 1.0*1;
N_grid = 2*25;
z_scan = 1;

[X,Y] = meshgrid(linspace(-scanplane_x_max, scanplane_x_max, N_grid), ...
                 linspace(-scanplane_y_max, scanplane_y_max, N_grid));

x_s = [X(:) Y(:) ones(N_grid*N_grid, 1)*z_scan];

% x_s(:, :, 1) = ones(N_grid, 1) * linspace(-scanplane_x_max, scanplane_x_max, N_grid);
% x_s(:, :, 2) = linspace(-scanplane_y_max, scanplane_y_max, N_grid)' * ones(1, N_grid);
% x_s(:, :, 3) = ones(N_grid, N_grid) * z_scan;

ind = 4;
f = frequencies(ind);
C = cpreal(ind,:,:) + 1i*cpimag(ind,:,:);
 
B = beamform(squeeze(C), f, x_s, x_m, 343.3);
imagesc(x_s(:,1), x_s(:,2), abs(B));


end

function B = beamform(C, f, x_s, x_m, c, form)

N_scan = size(x_s, 1);
N_mic = size(x_m, 1);
B = zeros(N_scan, 1);

for I = 1:N_scan
    r_0 = norm(x_s(I,:));
    r_s = sqrt(sum((ones(N_mic,1)*x_s(I,:) - x_m).^2, 2));
    h = h_1(r_s, r_0, f, c, N_mic);
    B(I) = h'*C*h;
    keyboard;
end

end

function h = h_1(r_s, r_0, f, c, N_mic)

h = exp(-1i*2*pi*f*(r_s-r_0)/c)/N_mic;


% h = exp(-1i*2*pi*f*(r_s-r_0)/c)./(r_s*r_0*sum(1./(r_s.^2)));

end