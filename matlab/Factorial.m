function [const] = Factorial(Q,N)
% The factorial calculation function, because the standard function does
% not allow values greater than 170.

const = 0;
for k = 1:N
    for l = 1:N
        if Q(k,l) <= 100
            const = const + log(factorial(Q(k,l)));
        else
            for t = 1:(Q(k,l)-100)
                const = const + log(Q(k,l)-(t-1));
            end
            const = const + log(factorial(100));
        end
    end
end

end

