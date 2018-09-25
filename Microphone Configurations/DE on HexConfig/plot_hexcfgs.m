clear all; clc; close all;
load('PerturbedHexCfg_F0.5_pc0.9.mat');

reso = get(0, 'screensize');
for i = 1:10
    figure('Visible', 'on', ...
           'Position', [floor(reso(3)/2)-525, floor(reso(4)/2)-250, 750, 375], ...
           'Resize', 'off');
    subplot(1,2,1); plot(F1(1:32, i, 1), F1(33:64, i, 1), 'b.');
    title(['J = ' num2str(G(i, 1))]); axis equal; axis([-1 1 -1 1]);
    subplot(1,2,2); plot(F1(1:32, i, end), F1(33:64, i, end), 'r.');
    title(['J = ' num2str(G(i, end))]); axis equal; axis([-1 1 -1 1]);
end

%% animate

x = 1:3999;

%// Plot starts here
figure, hold on

%// Set x and y limits of the plot
xlim([-1 1]);
ylim([-1 1]);

%// Plot point by point
for k = 1:numel(x)
    plot(F1(1:32, i, k), F1(33:64, i, k),'b.') %// Choose your own marker here
    axis equal; axis([-1 1 -1 1]);
    title(['Gen = ' num2str(k) ', J = ' num2str(G(i, k))]);

    %// MATLAB pauses for 0.001 sec before moving on to execue the next 
    %%// instruction and thus creating animation effect
    if k<50
        pause(0.25);
    elseif k<200
        pause(0.1);
    else
        pause(0.01);
    end
    cla(gca);
    if (k==1)
        set(gca, 'Box', 'on');
    end
end