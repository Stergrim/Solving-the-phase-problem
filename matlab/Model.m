function [Q] = Model(I,Norm,Noise)
% The image generation function is based on the intensity distribution I,
% the maximum on the Norm image and the noise matrix Noise.

% Image model formation
Q = Norm*I + Noise;

% Rounding to integer values to match the digital receiver
Q = round(Q);

end

