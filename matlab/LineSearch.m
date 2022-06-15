function [Lambda] = LineSearch(LimitIters,ErrorLambda,Lleft,Lright,C0,VectorPolinoms,R,Set,Norm,Noise,RealQ,B,Function)
% One-dimensional minimization function based on the golden section method.

% Calculation of the golden ratio
tau = (sqrt(5) - 1)/2;

Iters = 0;
k = 1;

Left = zeros(LimitIters,1);
Right = zeros(LimitIters,1);
L1 = zeros(LimitIters,1);
R1 = zeros(LimitIters,1);

% Defining the initial boundaries
Left(k) = Lleft;
Right(k) = Lright;

while abs(Right(k)-Left(k)) > ErrorLambda
    Iters = Iters + 1;
    
    if Iters > LimitIters
        break
    end
    
    L1(k) = Left(k) + (1-tau)*(Right(k)-Left(k));
    R1(k) = Left(k) + tau*(Right(k)-Left(k));
    
    if IncrementFunction(C0,VectorPolinoms,R,Set,Norm,Noise,RealQ,L1(k),B,Function) > IncrementFunction(C0,VectorPolinoms,R,Set,Norm,Noise,RealQ,R1(k),B,Function)
        Left(k+1) = L1(k);
        Right(k+1) = Right(k);
        L1(k+1) = R1(k);
        R1(k+1) = Left(k+1) + tau*(Right(k+1)-Left(k+1));
    else
        Left(k+1) = Left(k);
        Right(k+1) = R1(k);
        R1(k+1) = L1(k);
        L1(k+1) = Left(k+1) + (1-tau)*(Right(k+1)-Left(k+1));
    end
    k = k + 1;
end

Lambda = (Right(k) + Left(k))/2;

end

