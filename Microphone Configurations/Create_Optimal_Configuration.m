clear all;
%%
f = 2000;
c = 340;
xmin = -1; xmax = 1;
ymin = -1; ymax = 1;

d1 = 0.22;
% d1 = 1.38*c/f;
d2 = d1+0.04;
% d2 = 1.42*c/f+0.035;
d = [d1 d2];

x = zeros(32,2);
k = 1;

R1 = 0.15;
R2 = 0.45;
R3 = 0.75;
phi11 = 2*asin(d1/R1/2);
phi12 = 2*asin(d2/R1/2);

polars = zeros(33,2);
k = 2;
while (polars(k-1,1)+ rem(k,2)*phi11 + (~rem(k,2))*phi12) < 2*pi - rem(k,2)*phi11 - (~rem(k,2))*phi12
    polars(k,1) = polars(k-1,1) + rem(k,2)*phi11 + (~rem(k,2))*phi12;
    k = k + 1;
end
polars(1:k-1,2) = R1;

polars(k,2) = R2;
k = k+1;
phi11 = 2*asin(d1/R2/2);
phi12 = 2*asin(d2/R2/2);
while polars(k-1,1)< 2*pi - rem(k,2)*phi11 - (~rem(k,2))*phi12
    polars(k,1) = polars(k-1,1) + rem(k,2)*phi11 + (~rem(k,2))*phi12;
    polars(k,2) = R2;
    k = k + 1;
end

polars(k,2) = R3;
k = k+1;
phi11 = 2*asin(d1/R3/2);
phi12 = 2*asin(d2/R3/2);
while polars(k-1,1)< 2*pi - rem(k,2)*phi11 - (~rem(k,2))*phi12
    polars(k,1) = polars(k-1,1) + rem(k,2)*phi11 + (~rem(k,2))*phi12;
    polars(k,2) = R3;
    k = k + 1;
end

[x(:,1), x(:,2)] = pol2cart(polars(1:32,1),polars(1:32,2));
% polars(2:k-1,2) = R2;

figure; plot(x(:,1),x(:,2),'*'); axis equal; axis([-1 1 -1 1]);
%% hist
n_mic = 32;
A1 = [x(1:n_mic,1);x(1:n_mic,2)]; clear x;
neighbours = 2;
dist = zeros(n_mic, n_mic);

distances = zeros(n_mic, neighbours, 1);
mean_distances = zeros(neighbours, 1);

for j = 1:n_mic-1
    x = A1(j) - A1(j+1:n_mic);
    y = A1(j+n_mic) - A1(j+1+n_mic:end);
    dist(j, 1+j:n_mic) = sqrt(x.^2 + y.^2);
end
dist = dist + dist';
dist = sort(dist,1);

distances(:, :) = dist(2:neighbours+1,:)';
mean_distances(:) = mean(distances(:, :))';
% dist = zeros(n_mic, n_mic);