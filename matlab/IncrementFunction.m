function [Increment] = IncrementFunction(C0,VectorPolinoms,R,Set,Norm,Noise,RealQ,Lambda,B,Function)
% The function of calculating the increment of all Zernike coefficients.

C0 = C0 + Lambda*B;
Increment = Function(C0,VectorPolinoms,R,Set,Norm,Noise,RealQ);

end

