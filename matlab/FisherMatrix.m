function [Ff] = FisherMatrix(C0,KC0,VectorPolinoms,R,Set,Norm,Noise)
% The function of calculating the Fisher matrix.

% Disclosure of the parameter vector
M = Set(1);
N = Set(2);
Dp = Set(3);
z = Set(4);
wave = Set(5);
step = Set(6);
F = Set(7);

Ff = zeros(KC0,N,N);

% Wavefront formation
Front = WaveFront(VectorPolinoms,R,M,C0);

% Calculation of intensity distribution
[I,L,beta] = DirectTask(C0,VectorPolinoms,R,Set);

% Image model formation
Q = Model(I,Norm,Noise);

% Calculation of parameters g and A
g = wave*z/(Dp*step);
A = -(z-F)*(Dp^2)/(8*wave*F^2);

for k = 1:N
    for l = 1:N
        for a = 1:M
            for b = 1:M
                for t = 1:KC0
                    if R(a,b) <= 1
                        Ff(t,k,l) = Ff(t,k,l) + exp(-sqrt(-1)*2*pi*(((k-(N-1)/2-1)*(a-(M-1)/2-1)+(l-(N-1)/2-1)*(b-(M-1)/2-1))/(M*g)-Front(a,b)-A*R(a,b)^2))*VectorPolinoms(t,a,b);
                    end
                end
            end
        end
    end
end

% Matrix conjugation
L = conj(L);

for t = 1:KC0
    for k = 1:N
        for l = 1:N
            Ff(t,k,l) = 4*pi*beta*Norm*imag(L(k,l)*Ff(t,k,l));
            Ff(t,k,l) = Ff(t,k,l)/sqrt(Q(k,l) + 1);
        end
    end
end

end
