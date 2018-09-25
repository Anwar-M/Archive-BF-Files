clear all; clc;
F = 0.2; CR = 0.8;
NP = 64;      % number of population elements (q)
D  = 64;      % number of variables
n  = 10*NP;   

vlb = -1*ones(1, D);
vub = ones(1, D);
X1 = lhd_new(NP, vlb, vub);

popold = X1;

for i = 1:NP
    E1(i) = ObjectFunction(popold(i,:), D/2, 3.83, pi*2000/340);
end
    
for i = 1:40

    A = 1+floor(NP*0.9999999999999*rand(n,3));
    [I,J] = find(A(:,1)~=A(:,2) & A(:,1)~=A(:,3) & A(:,2)~=A(:,3));
    A = A(I,:); clear I J
    A = unique(A,'rows');
    n = length(A(:,1));
    M = zeros(n,NP);
    M(sub2ind(size(M),(1:n)',A(:,1))) = 1;
    M(sub2ind(size(M),(1:n)',A(:,2))) = F;
    M(sub2ind(size(M),(1:n)',A(:,3))) = -F;

    popinterm = M*popold;

    lowerbound = ones(n,1)*vlb;
    upperbound = ones(n,1)*vub;
    z1 = all(popinterm >= lowerbound,2);
    z2 = all(popinterm <= upperbound,2);
    z = z1 .* z2;
    I = find(z==1);
    if isempty(I)
        popinterm = popold;
    else
        popinterm = popinterm(I,:);
    end
    
    bounds_good(i) = sum(z);
%     if bounds_good(i)>0; keyboard; end
    
    n = length(popinterm(:,1));

    for p = 1 : NP
        % apply crossover from 1 vector of popold with all vectors of popinterm
        mui = rand(n,D) < CR;         % all random numbers < CR are 1, 0 otherwise
        mpo = mui < 0.5;              % inverse mask to mui: mui added to mpo gives a matrix consisting of ones
        % popinterm2 does not result in other parameter VALUES than those in pop_old and popinterm, but in other
        % combinations of the D parameter values
        popinterm2  = (ones(n,1)*popold(p,:)).*mpo + popinterm.*mui;
        % popinterm2 = setdiff(popinterm2, pop_total, 'rows'); % Camiel, 10 march 2004
        % select at random one element from the admissable set
        n2 = length(popinterm2(:,1));
        if (n2 == 0)
            disp('line 79 in create_descendants_1.m')
            keyboard
        end
        element = randperm(n2);
        popnew(p,:) = popinterm2(element(1),:);
        clear mui mpo popinterm2 n2 element

    end
    
    
    popnewer = popold;
    val = E1;
    for k=1:NP
%   tempval = ObjectFunction([fixed_parameters(:,1)' ui(k,1:n_mic-n_fixed_parameters) ...
%       fixed_parameters(:,2)' ui(k,n_mic-n_fixed_parameters+1:end)], ...
%       n_mic, alpha_min, alpha_max); % check cost of competitor  ui
        tempval = ObjectFunction(popnew(k,:), D/2, 3.83, pi*2000/340);

      if (abs(tempval) <= abs(E1(k)))    % if competitor is better than value in "cost array"
         popnewer(k,:) = popnew(k,:);               % replace old vector with new one (for new iteration)
         val(k)   = tempval;                  % save value in "cost array"
      end
    end
    
    popold = popnewer;
    E1 = val;
end
disp(sum(bounds_good))