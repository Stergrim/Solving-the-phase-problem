function [MSD] = MSDFunction(C0,VectorPolinoms,R,Set,Norm,Noise,RealQ)
% The function of calculating the deviation of two images.

N = Set(2);

% Calculation of intensity distribution
I = DirectTask(C0,VectorPolinoms,R,Set);
Q = Model(I,Norm,Noise);

MSD = 0;

for k = 1:N
    for l = 1:N
       MSD = MSD + (Q(k,l) - RealQ(k,l))^(2);
    end
end

MSD = sqrt(MSD)/N;

end

