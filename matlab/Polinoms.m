function [VectorPolinoms,R] = Polinoms(M,k)
% A function that forms a set of Zernike polynomials written into a vector.

% Ì - discretization of the exit pupil,  i.e. the number of partitions of
% the side of the square described around the exit pupil.
% k - number of polynomials.
% R - the matrix of values is the radius of the vector corresponding to the
% round pupil.
% Matrix R is necessary for subsequent functions to limit the scope of
% calculations.

VectorPolinoms = zeros(k,M,M);

% Cartesian coordinate system
x = zeros(1,M);
y = zeros(1,M);

% Polar coordinate system
R = zeros(M);
tt = zeros(M);

for a = 1:M
    for b = 1:M
        x(a) = -1 + 2*(a-1)/(M-1);
        y(b) = -1 + 2*(b-1)/(M-1);

        % Translation of the polar coordinate system to Cartesian
        R(a,b) = sqrt(x(a)^2 + y(b)^2); % Translation of the radius vector
        if R(a,b) <= 1 % The calculation is not carried out outside the
                       % circle of the unit radius
            % Translation of the polar angle
            if (x(a) > 0 && y(b) >= 0)||(x(a) > 0 && y(b) <= 0)
                tt(a,b) = atan(y(b)/x(a));
            elseif (x(a) < 0 && y(b) >= 0)||(x(a) < 0 && y(b) <= 0)
                tt(a,b) = pi + atan(y(b)/x(a));
            elseif (x(a) == 0 && x(a-1) > 0 && y(b) > 0)||...
(x(a) == 0 && x(a-1) < 0 && y(b) < 0)
                tt(a,b) = -pi/2;
            elseif (x(a) == 0 && x(a-1) < 0 && y(b) > 0)||...
(x(a) == 0 && x(a-1) > 0 && y(p) < 0)
                tt(a,b) = pi/2;
            end
            
            % The counter of the number of polynomials
            l = 1;
            
            % Formation of polynomials
            for n = 1:k
                for m = -n:2:n
                    if m < 0
                        for s = 0:(n+m)/2
                            VectorPolinoms(l,a,b) = VectorPolinoms(l,a,b)+sqrt(2*(n+1))*...
                                (((-1)^s)*factorial(n-s)/(factorial(s)*factorial((n-m)/2-s)*factorial((n+m)/2-s)))*...
                                sin(-m*tt(a,b))*R(a,b)^(n-2*s);
                        end
                    elseif m == 0
                        for s = 0:n/2
                            VectorPolinoms(l,a,b) = VectorPolinoms(l,a,b)+sqrt(n+1)*...
                                (((-1)^s)*factorial(n-s)/(factorial(s)*factorial(n/2-s)*factorial(n/2-s)))*...
                                R(a,b)^(n-2*s);
                        end
                    elseif m > 0
                        for s = 0:(n-m)/2
                            VectorPolinoms(l,a,b) = VectorPolinoms(l,a,b)+sqrt(2*(n+1))*...
                                (((-1)^s)*factorial(n-s)/(factorial(s)*factorial((n-m)/2-s)*factorial((n+m)/2-s)))*...
                                cos(m*tt(a,b))*R(a,b)^(n-2*s);
                        end
                    end
                    
                    l = l + 1;
                    
                    if l > k
                        break
                    end
                end
                
                if l > k
                    break
                end
            end
        end
    end
end

end
