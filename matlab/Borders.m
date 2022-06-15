function [Lmax] = Borders(Grad,MaxStep,MinLimit)
% The function of determining the boundary value of the Lmax variable Lambda.

if abs(max(Grad)) >= abs(min(Grad))
    Max = abs(max(Grad));
else
    Max = abs(min(Grad));
end

if Max < MinLimit
    Lmax = 0;
else
    Lmax = MaxStep/Max;
end

end

