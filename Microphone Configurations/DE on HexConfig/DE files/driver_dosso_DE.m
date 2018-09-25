% Driver for 'dosso' minimisation

clear all

t1 = cputime;

% Search bounds
vlb = -3*ones(1,6);
vub =  3*ones(1,6);
% Energy function
energyfun = 'dosso';
% Crossover rate
pc = 0.55;
% Multiplication factor
F = 0.6;
% Number of runs
Nruns = 10;
% Population size
q = 16;
% Number of generations
Ng = 150;

for k3 = 1 : Nruns
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
    for k4 = 1:q
        E1(k4) = dosso(X1(k4,:),vlb,vub);
    end
    % tranfer to column vector
    E1 = E1';
    pop_total = [pop_total; X1];
    energy_total = [energy_total E1'];
    for k5 = 1:Ng-1
        % create_descendants_DE_1
        % OUTPUT:
        % dummy1 = popinterm, i.e. the partner population
        % n = size of partner population
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
        [X4, dummy2, pop_total, energy_total] = new_generation_DE_1('dosso', X1, E1, X3, vlb, vub, pop_total, energy_total);
        E1 = dummy2;
        % update initial population
        X1 = X4;
    end
    [dummy,Imin] = min(E1);
    F1(:,k3) = X1(Imin,:)';
    G(k3) = E1(Imin);
    t2 = cputime;
    t2 - t1
    t1 = t2;
end

% end of program