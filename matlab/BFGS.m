function [C0Inv,Iters,MSDInv] = BFGS(C0Start,Function,SetInverse,VectorPolinoms,R,SetDirect,Norm,Noise,KC0,RealQ)
% The function of the implementation of the BFGS algorithm.

% Disclosure of the parameter vector
MaxIters = SetInverse(1);
MaxStep = SetInverse(2);
LeftBord = SetInverse(3);
RightBord = SetInverse(4);
MSD = SetInverse(5);
ErrorCoeff = SetInverse(6);

% Formation of a new vector of coefficients, so as not to re-assign the
% input coefficients
C0Inv = C0Start;

% Formation of a unit matrix that performs the role of the inverse Hesse matrix
H = eye(KC0);

Iters = 0;

% Calculation of the gradient of the likelihood function
Grad1 = Derivative(C0Inv,KC0,RealQ,VectorPolinoms,R,SetDirect,Norm,Noise);

flag1 = 0;

while flag1 == 0
    Iters = Iters + 1;
    if Iters > MaxIters
        break
    end
    
    % Formation of a vector of directions of movement of changes in coefficients
    B = -1*H*Grad1;
    
    % Determination of the limits of the step coefficient, in order not to
    % go beyond the range of acceptable values
    Lmax = Borders(Grad1,MaxStep,10^(-3));
    
    if Lmax == 0
        break
    end
    
    Lleft = -Lmax;
    Lright = Lmax;
    
    % Calculation of the step coefficient by the method of one-dimensional
    % minimization of the function
    Lambda = LineSearch(10,10^(-10),Lleft,Lright,C0Inv,VectorPolinoms,R,SetDirect,Norm,Noise,RealQ,B,Function);
    
    % Formation of a vector of increments of C0Inv coefficients
    d = Lambda*B;
    
    % Restriction on increment of coefficients
    for t = 1:KC0
        if abs(d(t)) > MaxStep && d(t) > 0
            d(t) = MaxStep;
        elseif abs(d(t)) > MaxStep && d(t) < 0
            d(t) = -MaxStep;
        end
    end
    
    % Formation of new coefficients
    C0Inv = C0Inv + d;
    
    flag2 = 0;
    
    % If at least one of the coefficients reaches the boundary, then we
    % finish the work
    for t = 1:KC0
        if C0Inv(t) <= LeftBord
            C0Inv(t) = LeftBord;
            flag2 = 1;
        elseif C0Inv(t) >= RightBord
            C0Inv(t) = RightBord;
            flag2 = 1;
        end
    end
    
    if flag2 == 1
        break
    end
    
    % Calculation of the mean square deviation (MSD) of two images: real
    % and calculated
    MSDInv = MSDFunction(C0Inv,VectorPolinoms,R,SetDirect,Norm,Noise,RealQ);
    
    % If the MSD is less than the specified one, then we stop working
    if MSDInv <= MSD
        break
    end
    
    % Calculation of the new gradient
    Grad2 = Derivative(C0Inv,KC0,RealQ,VectorPolinoms,R,SetDirect,Norm,Noise);
    
    DeltaGrad = Grad2 - Grad1;
    Grad1 = Grad2;
    Indicator = 0;
    
    % If the step of changing all coefficients is less than the specified
    % one, then we terminate the work
    for t = 1:KC0
        if abs(d(t)) <= ErrorCoeff
            Indicator = Indicator + 1;
        end
    end
    
    if Indicator == KC0
        break
    end
    
    % Calculate a new approximate inverse Hesse matrix
    ddT = d*d';
    dTD = d'*DeltaGrad;
    if dTD == 0
        break
    end
    
    DTHD = DeltaGrad'*H*DeltaGrad;
    dDTH = d*DeltaGrad'*H;
    HDdT = H*DeltaGrad*d';
    H = H + (1 + DTHD/dTD)*ddT/dTD - dDTH/dTD - HDdT/dTD;
    
    % At iterations 4, 8, 12 and 16, we update the inverse Hesse matrix,
    % that is, we form the unit matrix again. It is optimal with a maximum
    % number of iterations equal to 20.
    if Iters == 4 || Iters == 8 || Iters == 12 || Iters == 16
        H = eye(KC0);
    end
end

end

