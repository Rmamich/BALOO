function [Xk] = Xk_relInertOnly_accelNoise(dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Takes in dt and returns the matrix necessary to
% transform the diagonal of the power pectral density matrix to the vector
% representation of the upper triangular portion of the discrete time
% matrix.
% 
% Assumption: Assumes that noise is white noise applied as jerk
% 
% Inputs:
%   dt      change in time for integration [s] (tkp1 - tk)
% 
% Outputs:
%   Xk      matrix to transfrom diahonal of psd to discrete time
%               (vectorized upper tiangular)
% 
% Author: Rachel Mamich 
% Date Created: 10/20/23
% Date Modified: 10/20/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xk = zeros(45, 3)*dt;

indexes_dt5div36 = [1, 10, 18];
indexes_dt4div12 = [4, 13, 21];
indexes_dt3div6 = [7, 16, 24];
indexes_dt3div4 = [25, 31, 36];
indexes_dt2div2 = [28, 34, 39];
indexes_dt = [40, 43, 45];


Xk(indexes_dt5div36, 1:3) = 1/3*dt^3*eye(3);
Xk(indexes_dt4div12, 1:3) = 1/2*dt^2*eye(3);
Xk(indexes_dt3div6, 1:3) = 1/2*dt*eye(3);
Xk(indexes_dt3div4, 1:3) = dt*eye(3);
Xk(indexes_dt2div2, 1:3) = eye(3);
Xk(indexes_dt, 1:3) = 1/dt*eye(3);





end