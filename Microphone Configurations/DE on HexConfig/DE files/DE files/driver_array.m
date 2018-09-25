function [pop_total, energy_total, F1, G] = driver_array(pc, F, q, Ng, Nruns, alpha)

Nmics = 32;
bound = 1;

vlb = bound*ones(1, 2*Nmics);
vub = -1*bound*ones(1, 2*Nmics);

for i = 1 : Nruns
    tic;
%     rng('shuffle');
    seed_nr = abs( sum(100*clock) - cputime );
    rand('state',seed_nr) % set the random generator to a different state each time
    clear seed_nr
    
    % run DE
    pop_total = [];
    energy_total = [];
    clear E1
    startpop = lhd_new(q, vlb, vub);
    X1 = startpop;
    % energy function values for elements of X1
    for j = 1:q
        E1(j) = ObjectFunction(X1(j, :), Nmics, alpha(1), alpha(2));
    end
    
    % tranfer to column vector
    E1 = E1';
    pop_total = [pop_total; X1];
    energy_total = [energy_total E1'];
    for j = 1:Ng-1
        % create_descendants_DE_1
        % OUTPUT:
        % dummy1 = popinterm, i.e. the partner population
        % Nt = size of partner population
        % X3 = popnew, i.e., the new population obtained from
        % crossover between popold and popinterm
        [dummy1, Nt, X3] = create_descendants_DE_1(X1, F, pc, vlb, vub, pop_total);
        % new_generation_DE_1
        % OUTPUT:
        % X4 = popnew, i.e. the population for the next generation
        % dummy2 are the corresponding values for the energy
        % function
        % pop_total and energy_total are updated with respect to
        % input values
        [X4, dummy2, pop_total, energy_total] = new_generation_DE_2(X1, E1, X3, alpha(1), alpha(2), pop_total, energy_total);
        E1 = dummy2;
        % update initial population
        X1 = X4;
    end
    [dummy, Imin] = min(E1);
    F1(:, i) = X1(Imin, :)';
    G(i) = E1(Imin);
    disp(toc);
end

end