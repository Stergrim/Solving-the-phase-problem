function [Likelihood] = LikelihoodFunction(C0,VectorPolinoms,R,Set,Norm,Noise,RealQ)
% The likelihood function.

N = Set(2);

% Calculation of intensity distribution
I = DirectTask(C0,VectorPolinoms,R,Set);
Q = Model(I,Norm,Noise);

% Calculation of the constant
const = Factorial(RealQ + 1,N);
Likelihood = const;

for k = 1:N
    for l = 1:N
       Likelihood = Likelihood + (Q(k,l) + 1) - (RealQ(k,l) + 1)*log(Q(k,l) + 1);
    end
end

end

