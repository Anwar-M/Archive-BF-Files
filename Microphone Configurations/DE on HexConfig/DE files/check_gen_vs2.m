% check_gen checks whether a population is in the admissable region

function z = check_gen_vs2(popnew, vlb, vub)

z = 1;
NP = size(popnew,1);    % number of population elements
D = size(popnew,2);     % number of dimensions

mvlb = ones(NP,1)*vlb;
mvub = ones(NP,1)*vub;
z1 = all(popnew >= mvlb);
z2 = all(popnew <= mvub);
if all(z1)==0 | all(z2)==0
   z=0;
end
