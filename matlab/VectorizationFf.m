function [VectorM] = VectorizationFf(Matrix,N,KC0)
% The vectorization function of a three-dimensional matrix.

VectorM = zeros(N^(2),KC0);

for t = 1:KC0
    a = 1;
    for k = 1:N
        for l = 1:N
            VectorM(a,t) = Matrix(t,k,l);
            a = a + 1;
        end
    end
end

end

