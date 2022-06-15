function [I,L,beta] = DirectTask(C0,VectorPolinoms,R,Set)
% The function of calculating the intensity distribution in the image plane.

% Disclosure of the parameter vector
M = Set(1);
N = Set(2); % Image size in pixels
Dp = Set(3); % Exit pupil diameter in meters
z = Set(4); % Distance to the image plane in meters
wave = Set(5); % Wavelength in meters
step = Set(6); % Calculation step or receiver pixel size in meters
F = Set(7); % Focal length of the system in meters

% Wavefront formation
Front = WaveFront(VectorPolinoms,R,M,C0);

L = zeros(N);
I = zeros(N);

% Calculation of parameters g and A
g = wave*z/(Dp*step);
A = -(z-F)*(Dp^2)/(8*wave*F^2);

% Calculation of intensity distribution
for k = 1:N
    for l = 1:N
        for a = 1:M
            for b = 1:M
                if R(a,b) <= 1
                    L(k,l) = L(k,l) + exp(-sqrt(-1)*2*pi*(((k-(N-1)/2-1)*(a-(M-1)/2-1)+(l-(N-1)/2-1)*(b-(M-1)/2-1))/(M*g)-Front(a,b)-A*R(a,b)^2));
                end
            end
        end
        I(k,l) = (abs(L(k,l)))^2;
    end
end

% Calculation of the normalization coefficient
beta = 1/max(I(:));

% Normalization of the intensity distribution
I = beta*I;

end