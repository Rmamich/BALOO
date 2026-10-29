function [A] = gravJac(x, gravParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Calculated the jacobian of the acceleration components wrt the
% position components for TBP, J2-J4
% 
% Inputs:
%   x           state vector
%                   x(1:3)     chief position in the inertial frame [m] 
%                   x(4:6)     chief velocity in the inertial frame [m/s]
%                   x(7:9)     chief unmodeled acceleration [m/s^2]
%                   x(10:12)   deputy position in the inertial frame [m]
%                   x(13:15)   deputy velocity in the inertial frame [m/s]
%                   x(16:18)   deputy unmodeled acceleration est [m/s^2]
%   gravParam   Array of elements for TPB, J2-4
%                   (1)         mu [m^3/s^2]
%                   (2)         R equatorial radius of central body [m]
%                   (3)         J2
%                   (4)         J3
%                   (5)         J4
% 
% Outputs:
%   A           jacobian of the acceleration components wrt the
%                   position components for TBP, J2-J4
% 
% Author: Rachel Mamich
% Date Modified: 10/18/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = gravParam(1);
R = gravParam(2);
J2 = gravParam(3);
J3 = gravParam(4);
J4 = gravParam(5);
rx = x(1);
ry = x(2);
rz = x(3);


r = [rx; ry; rz];

rexp3 = (r'*r)^(3/2);
rexp5 = (r'*r)^(5/2);
rexp7 = (r'*r)^(7/2);
rexp9 = (r'*r)^(9/2);
rexp11 = (r'*r)^(11/2);
rexp13 = (r'*r)^(13/2);

AJ4x = [15/8*J4*R^4*mu*((5*rx^4 + 6*rx^2*ry^2 - 36*rx^2*rz^2 + ry^4 -12*ry^2*rz^2 + 8*rz^4)/rexp11 -11*(rx^5 + 2*rx^3*ry^2 - 12*rx^3*rz^2 + rx*ry^4 -12*ry^2*rz^2*rx + 8*rx*rz^4)*rx/rexp13), 15/8*J4*R^4*mu*((4*rx^3*ry + 4*rx*ry^3 - 24*ry*rz^2*rx)/rexp11 -11*(rx^5 + 2*rx^3*ry^2 - 12*rx^3*rz^2 + rx*ry^4 -12*ry^2*rz^2*rx + 8*rx*rz^4)*ry/rexp13), 15/8*J4*R^4*mu*((-24*rx^3*rz - 24*ry^2*rz*rx + 32*rx*rz^3)/rexp11 -11*(rx^5 + 2*rx^3*ry^2 - 12*rx^3*rz^2 + rx*ry^4 -12*ry^2*rz^2*rx + 8*rx*rz^4)*rz/rexp13)];
AJ4y = [15/8*J4*R^4*mu*((4*rx^3*ry + 4*rx*ry^3 - 24*rx*ry*rz^2)/rexp11 - 11*(rx^4*ry + 2*rx^2*ry^3 - 12*rx^2*ry*rz^2 + ry^5 - 12*ry^3*rz^2 + 8*ry*rz^4)*rx/rexp13), 15/8*J4*R^4*mu*((rx^4 + 6*rx^2*ry^2 -12*rx^2*rz^2 + 5*ry^4 - 36*ry^2*rz^2+8*rz^4)/rexp11 - 11*(rx^4*ry + 2*rx^2*ry^3 - 12*rx^2*ry*rz^2 + ry^5 - 12*ry^3*rz^2 + 8*ry*rz^4)*ry/rexp13), 15/8*J4*R^4*mu*((-24*rx^2*ry*rz -24*ry^3*rz + 32*ry*rz^3)/rexp11 - 11*(rx^4*ry + 2*rx^2*ry^3 - 12*rx^2*ry*rz^2 + ry^5 - 12*ry^3*rz^2 + 8*ry*rz^4)*rz/rexp13)];
AJ4z = [5/8*J4*R^4*mu*((60*rx^3*rz + 60*rx*ry^2*rz - 80*rx*rz^3)/rexp11 - 11*(15*rx^4*rz+30*rx^2*ry^2*rz -40*rx^2*rz^3 + 15*ry^4*rz-40*ry^2*rz^3+8*rz^5)*rx/rexp13), 5/8*J4*R^4*mu*((60*rx^2*ry*rz + 60*ry^3*rz -80*ry*rz^3)/rexp11 - 11*(15*rx^4*rz+30*rx^2*ry^2*rz -40*rx^2*rz^3 + 15*ry^4*rz-40*ry^2*rz^3+8*rz^5)*ry/rexp13), 5/8*J4*R^4*mu*((15*rx^4 + 30*rx^2*ry^2 - 120*rx^2*rz^2+15*ry^4-120*ry^2*rz^2+40*rz^4)/rexp11 - 11*(15*rx^4*rz+30*rx^2*ry^2*rz -40*rx^2*rz^3 + 15*ry^4*rz-40*ry^2*rz^3+8*rz^5)*rz/rexp13)];


AJ3x = [-5/2*J3*R^3*mu*((9*rx^2*rz + 3*ry^2*rz - 4*rz^3)/rexp9 -9*(3*rx^3*rz + 3*ry^2*rx*rz - 4*rx*rz^3)*rx/rexp11), -5/2*J3*R^3*mu*((6*ry*rx*rz)/rexp9 -9*(3*rx^3*rz + 3*ry^2*rx*rz - 4*rx*rz^3)*ry/rexp11), -5/2*J3*R^3*mu*((3*rx^3 + 3*ry^2*rx - 12*rx*rz^2)/rexp9 -9*(3*rx^3*rz + 3*ry^2*rx*rz - 4*rx*rz^3)*rz/rexp11)];
AJ3y = [-5/2*J3*R^3*mu*(6*rx*ry*rz/rexp9 -9*(3*rx^2*ry*rz + 3*ry^3*rz - 4*ry*rz^3)*rx/rexp11), -5/2*J3*R^3*mu*((3*rx^2*rz + 9*ry^2*rz - 4*rz^3)/rexp9 -9*(3*rx^2*ry*rz + 3*ry^3*rz - 4*ry*rz^3)*ry/rexp11), -5/2*J3*R^3*mu*((3*rx^2*ry + 3*ry^3 - 12*ry*rz^2)/rexp9 -9*(3*rx^2*ry*rz + 3*ry^3*rz - 4*ry*rz^3)*rz/rexp11)];
AJ3z = [1/2*J3*R^3*mu*((12*rx^3 + 12*rx*ry^2 - 48*rx*rz^2)/rexp9 - 9*(3*rx^4 + 6*rx^2*ry^2 - 24*rx^2*rz^2 + 3*ry^4 - 24*ry^2*rz^2 + 8*rz^4)*rx/rexp11), 1/2*J3*R^3*mu*((12*rx^2*ry + 12*ry^3 - 48*ry*rz^2)/rexp9 - 9*(3*rx^4 + 6*rx^2*ry^2 - 24*rx^2*rz^2 + 3*ry^4 - 24*ry^2*rz^2 + 8*rz^4)*ry/rexp11), 1/2*J3*R^3*mu*((-48*rx^2*rz - 48*ry^2*rz + 32*rz^3)/rexp9 - 9*(3*rx^4 + 6*rx^2*ry^2 - 24*rx^2*rz^2 + 3*ry^4 - 24*ry^2*rz^2 + 8*rz^4)*rz/rexp11)];


AJ2x = 3/2 * J2 * mu * R^2*[5*rz^2/rexp7-35*rz^2*rx^2/rexp9-1/rexp5+5*rx^2/rexp7, -35*rz^2*rx*ry/rexp9 + 5*rx*ry/rexp7,                 5*rx*rz/rexp7 + 10*rx*rz/rexp7 - 35*rz^3*rx/rexp9];
AJ2y = 3/2 * J2 * mu * R^2*[5*ry*rx/rexp7 - 35*rz^2*ry*rx/rexp9,                  5*rz^2/rexp7-35*rz^2*ry^2/rexp9-1/rexp5+5*ry^2/rexp7, 5*ry*rz/rexp7 + 10*ry*rz/rexp7 - 35*rz^3*ry/rexp9];
AJ2z = 3/2 * J2 * mu * R^2*[15*rz*rx/rexp7 - 35*rz^3*rx/rexp9,                    15*rz*ry/rexp7 - 35*rz^3*ry/rexp9,                    15*rz^2/rexp7 - 35*rz^4/rexp9-3/rexp5+15*rz^2/rexp7];


Atbp_x = [-mu/rexp3 + 3*mu/rexp5*rx^2, 3*mu/rexp5*rx*ry, 3*mu/rexp5*rx*rz];
Atbp_y = [3*mu/rexp5*rx*ry, -mu/rexp3+3*mu/rexp5*ry^2,3*mu/rexp5*ry*rz];
Atbp_z = [3*mu/rexp5*rx*rz, 3*mu/rexp5*ry*rz, -mu/rexp3+3*mu/rexp5*rz^2];

AJ4 = [AJ4x; AJ4y; AJ4z];
AJ3 = [AJ3x; AJ3y; AJ3z];
AJ2 = [AJ2x; AJ2y; AJ2z];
ATBP = [Atbp_x; Atbp_y; Atbp_z];

A = AJ4 + AJ3 + AJ2 + ATBP;

end