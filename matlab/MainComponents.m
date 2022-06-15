function [p] = MainComponents(C0Inv,VT,KC0)
% The function of calculating the principal components from the singular
% value decomposition of the Fisher matrix.

p = zeros(KC0,1);

for t = 1:KC0
    p(t) = VT(1,:)*C0Inv;
end

end

