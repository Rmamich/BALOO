function [A] = gravJac_CRTBP_page(x, gravParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Calculated the jacobian of the acceleration components wrt the
% position components for TBP, J2-J4
% 
% Inputs:
%   x           state vector
%                   x(1:3)     nondimensional position in rotating frame []
%   gravParam   Array of elements for TPB, J2-4
%                   (1)         mu []
%                   (2)         
% 
% Outputs:
%   A           jacobian of the acceleration components of CRTBP wrt the
%                   position components for CRTBP
% 
% Author: Rachel Mamich
% Date Modified: 07/31/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N, n_bank] = size(x);
x = reshape(x, N, 1, n_bank);
mu = gravParam(1);
rx = x(1, :, :);
ry = x(2, :, :);
rz = x(3, :, :);


r1Vec = [rx + mu; ry; rz];
r2Vec = [rx-1+mu; ry; rz];


r1 = pagenorm(r1Vec);
r1expneg3 = r1.^(-3);
r1expneg5 = r1.^(-5);

r2 = pagenorm(r2Vec);
r2expneg3 = r2.^(-3);
r2expneg5 = r2.^(-5);


xmu = rx+mu;
x1mu = rx-1+mu;


A = zeros(3, 3, n_bank);
A(1, 1, :) = 1 - (1-mu).*r1expneg3 + 3*(1-mu)*pagemtimes((xmu).^2,r1expneg5) - mu*r2expneg3 + 3*mu*pagemtimes((x1mu).^2,r2expneg5);
A(1, 2, :) = 3*(1-mu)*pagemtimes((xmu),pagemtimes(ry,r1expneg5)) + 3*mu*pagemtimes(ry,pagemtimes((x1mu),r2expneg5));
A(1, 3, :) = 3*(1-mu)*pagemtimes((xmu),pagemtimes(rz,r1expneg5)) + 3*mu*pagemtimes(rz,pagemtimes((x1mu),r2expneg5));
A(2, 1, :) = A(1, 2);
A(2, 2, :) = 1 - (1-mu)*r1expneg3 + 3*(1-mu)*pagemtimes(ry.^2,r1expneg5) - mu*r2expneg3 + 3*mu*pagemtimes(ry.^2,r2expneg5);
A(2, 3, :) = 3*(1-mu)*pagemtimes(ry,pagemtimes(rz,r1expneg5)) + 3*mu*pagemtimes(ry,pagemtimes(rz,r2expneg5));
A(3, 1, :) = A(1, 3);
A(3, 2, :) = A(2, 3);
A(3, 3, :) = -(1-mu)*r1expneg3 + 3*(1-mu)*pagemtimes(rz.^2,r1expneg5) - mu*r2expneg3 + 3*mu*pagemtimes(rz.^2,r2expneg5);







end