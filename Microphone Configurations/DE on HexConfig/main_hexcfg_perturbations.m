% SMALL PERTURBATIONS ON HEXAGONAL CONFIG FOR START POPULATION
clear all; clc;
addpath('DE files');
%%
matlabpool('open',feature('numCores'));

sigma = 0.0384;

q = 64; % population size

p_c_values = 0.9; % crossover probability
F_values = 0.5; % multiplication factor to 2

N_G = 4000; % number of generations
repeat = 10;

% # mics
n_mic = 32;

% array dimensions
x_min = -1;
x_max = 1;
y_min = -1;
y_max = 1;

% Objective function parameters
f_max = 2000;
D = (x_max - x_min);
alpha_max = pi*f_max/340;
alpha_min = 3.83/(D/2);

vlb = x_min*ones(1, n_mic*2);
vub = x_max*ones(1, n_mic*2);

for i = 1:length(p_c_values)
    p_c = p_c_values(i);
    for j = 1:length(F_values)
        F = F_values(j);
        
        fprintf('\n\t p_c=%.2f, F=%.2f \n', p_c, F)
        
        tic;
        load('config_data.mat');
        unperturbed_hexconfig_pop = ones(q, 1) * A1(:, 2)';
        parfor k = 1:repeat
%             seed_nr = abs( sum(100*clock) - cputime );
%             rand('state',seed_nr) % set the random generator to a different state each time
%             clear seed_nr

            % Perturbed initial population based on hex config
            startpop = unperturbed_hexconfig_pop + normrnd(0, sigma/2, q, n_mic*2);
            X1 = startpop;
            
            E1 = zeros(1, q);
            for k4 = 1:q
                E1(k4) = ObjectFunction(X1(k4, :), n_mic, alpha_min, alpha_max);
            end

            E1 = E1';

            pop_total = [];
            energy_total = [E1'];
            for k5 = 1:N_G-1
                [dummy1, Nt, X3] = create_descendants_DE_1(X1, F, p_c, vlb, vub, pop_total);
                [X4, dummy2, pop_total, energy_total] = new_generation_DE_2(X1, E1, X3, alpha_min, alpha_max, pop_total, energy_total);
                E1 = dummy2;
                % update initial population
                X1 = X4;

                [G(k,k5), Imin] = min(E1);
                F1(:, k, k5) = X1(Imin, :)';
            end
            
        end
        disp(toc);
        file_name = ['PerturbedHexCfg_F' num2str(F) '_pc' num2str(p_c) '.mat'];
        save(file_name, '-v7.3');
        disp('F: OK');
    end
    disp('p_c: OK');
end

matlabpool('close');
