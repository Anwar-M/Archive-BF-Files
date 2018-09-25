function [popnew, val, pop_total, energy_total] = new_generation_DE_2(popold, fvalold, ui, alpha_min, alpha_max, pop_total, energy_total)

NP = size(popold,1);    % number of population elements (q)
n_mic = size(popold, 2)/2;
popnew = popold;
val = fvalold;

%------------------------------------------------------------------------------------------------------
% Select which vectors are allowed to enter the new population (compare popold and ui)
%------------------------------------------------------------------------------------------------------
for k=1:NP
  tempval = ObjectFunction(ui(k,:), n_mic, alpha_min, alpha_max); % check cost of competitor  ui

  if (abs(tempval) <= abs(fvalold(k)))    % if competitor is better than value in "cost array"
     popnew(k,:) = ui(k,:);               % replace old vector with new one (for new iteration)
     val(k)   = tempval;                  % save value in "cost array"
  end
  
  % save all evaluated vectors with their energies (10 march 2004, Camiel)
  pop_total = [pop_total; ui(k,:)];
  energy_total = [energy_total tempval];
  clear tempval
end

% end of program