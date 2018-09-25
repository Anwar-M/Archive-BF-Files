clear all; clc;

pc = 0.55;
F = 0.7;
q = 64;
Nruns = 1;
Ng = 50;

alpha = [3.83/2 pi*2000/343];

[pop_total, energy_total, F1, G] = driver_array(pc, F, q, Ng, Nruns, alpha);

%%
for i = 1:N_G
z(i) = min(energy_total(q*(i-1)+1:q*i));
end
figure;
plot(z);