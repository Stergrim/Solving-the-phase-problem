function [Front] = WaveFront(VectorPolinoms,R,M,C0)
% A function that calculates the wavefront.

% C0 - vector of Zernike coefficients.

Front = zeros(M);

% Determining the size of the coefficient vector, whether the vector is a
% column or a row vector.
Size_C = size(C0);

% If the vector is a column, then transpose.
if Size_C(1) > Size_C(2)
    C0 = C0';
end

% Summation of polynomials multiplied by the corresponding coefficients.
for a = 1:M
    for b = 1:M
        if R(a,b) <= 1
            Front(a,b) = C0*VectorPolinoms(:,a,b);
        end
    end
end

end
