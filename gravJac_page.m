function [A] = gravJac_page(x, gravParam)
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
[N, n_bank] = size(x);

x = reshape(x, N, 1, n_bank);

rx = x(1, :, :);
ry = x(2, :, :);
rz = x(3, :, :);

r = x(1:3, :, :);

rNorm =  pagenorm(r);
rexp3 = rNorm.^3;
rexp5 = rNorm.^5;
rexp7 = rNorm.^7;
rexp9 = rNorm.^9;
rexp11 = rNorm.^11;
rexp13 = rNorm.^13;

rx2 = rx.^2;
rx3 = rx.^3;
rx4 = rx.^4;
rx5 = rx.^5;


ry2 = ry.^2;
ry3 = ry.^3;
ry4 = ry.^4;
ry5 = ry.^5;


rz2 = rz.^2;
rz3 = rz.^3;
rz4 = rz.^4;
rz5 = rz.^5;



rxry = pagemtimes(rx, ry);
rxrz = pagemtimes(rx,rz);


rx2ry = pagemtimes(rx2, ry);
rxrz3 = pagemtimes(rx, rz3);

rx2rz = pagemtimes(rx2,rz);
rx2ry2 = pagemtimes(rx2, ry2);
rx2rz2 = pagemtimes(rx2, rz2);
rx2ry3 = pagemtimes(rx2,ry3);
rx2rz3 = pagemtimes(rx2,rz3);


ryrz = pagemtimes(ry, rz);
ryrz2 = pagemtimes(ry, rz2);
ryrz3 = pagemtimes(ry, rz3);
ryrz4 = pagemtimes(ry,rz4);


ry2rz = pagemtimes(ry2,rz);
ry2rz2 = pagemtimes(ry2, rz2);
ry2rz3 = pagemtimes(ry2,rz3);




rx3ry2 = pagemtimes(rx3, ry2);

rx3rz = pagemtimes(rx3, rz);

rx3rz2 = pagemtimes(rx3, rz2);

ry3rz = pagemtimes(ry3,rz);
ry3rz2 = pagemtimes(ry3,rz2);

rx4ry = pagemtimes(rx4,ry);
rx4rz = pagemtimes(rx4,rz);

ry4rz = pagemtimes(ry4,rz);



J4xy_const = 15/8*J4*R^4*mu;
J4z_const = 5/8*J4*R^4*mu;

AJ4x = J4xy_const*[((5*rx4 + 6*rx2ry2 - 36*rx2rz2 + ry4 -12*ry2rz2 + 8*rz4) ./ rexp11 - 11*pagemtimes((rx5 + 2*rx3ry2 - 12*rx3rz2 + pagemtimes(rx,ry4) -12*pagemtimes(ry2rz2,rx) + 8*pagemtimes(rx,rz4)), rx ./ rexp13)),...
    ((4*pagemtimes(rx3,ry) + 4*pagemtimes(rx,ry3) - 24*pagemtimes(ryrz2,rx)) ./rexp11 - 11*pagemtimes((rx5 + 2*rx3ry2 - 12*rx3rz2 + pagemtimes(rx,ry4) -12*pagemtimes(ry2rz2,rx) + 8*pagemtimes(rx,rz4)), ry./rexp13)),...
    ((-24*pagemtimes(rx3,rz) - 24*pagemtimes(ry2rz,rx) + 32*pagemtimes(rx,rz3)) ./rexp11 -11*pagemtimes((rx5 + 2*rx3ry2 - 12*rx3rz2 + pagemtimes(rx,ry4) -12*pagemtimes(ry2rz2,rx) + 8*pagemtimes(rx,rz4)), rz./rexp13))];


AJ4y = J4xy_const*[((4*pagemtimes(rx3,ry) + 4*pagemtimes(rx,ry3) - 24*pagemtimes(rxry,rz2)) ./rexp11 - 11*pagemtimes((pagemtimes(rx4,ry) + 2*rx2ry3 - 12*pagemtimes(rx2ry,rz2) + ry5 - 12*ry3rz2 + 8*pagemtimes(ry,rz4)), rx./rexp13)),...
    ((rx4 + 6*rx2ry2 -12*pagemtimes(rx2,rz2) + 5*ry4 - 36*ry2rz2 + 8*rz4) ./rexp11 - 11*pagemtimes((pagemtimes(rx4,ry) + 2*rx2ry3 - 12*pagemtimes(rx2ry,rz2) + ry5 - 12*ry3rz2 + 8*pagemtimes(ry,rz4)), ry./rexp13)),...
    ((-24*pagemtimes(rx2ry,rz) -24*ry3rz + 32*ryrz3)./rexp11 - 11*pagemtimes((rx4ry + 2*rx2ry3 - 12*pagemtimes(rx2ry,rz2) + ry5 - 12*ry3rz2 + 8*ryrz4), rz./rexp13))];


AJ4z = J4z_const*[((60*rx3rz + 60*pagemtimes(rx,ry2rz) - 80*rxrz3) ./rexp11 - 11*pagemtimes((15*rx4rz + 30*pagemtimes(rx2ry2,rz) - 40*rx2rz3 + 15*ry4rz - 40*ry2rz3 + 8*rz5), rx ./rexp13)),...
    ((60*pagemtimes(rx2ry,rz) + 60*ry3rz -80*ryrz3) ./rexp11 - 11*pagemtimes((15*rx4rz + 30*pagemtimes(rx2ry2,rz) -40*rx2rz3 + 15*ry4rz - 40*ry2rz3 + 8*rz5),ry ./rexp13)),...
    ((15*rx4 + 30*rx2ry2 - 120*pagemtimes(rx2,rz2) + 15*ry4 - 120*pagemtimes(ry2,rz2) + 40*rz4) ./ rexp11 - 11*pagemtimes((15*rx4rz + 30*pagemtimes(rx2ry2,rz) - 40*rx2rz3 + 15*ry4rz - 40*ry2rz3 + 8*rz5), rz ./rexp13))];






J3xy_const = -5/2*J3*R^3*mu;
J3z_const = 1/2*J3*R^3*mu;
AJ3x = J3xy_const*[((9*rx2rz + 3*ry2rz - 4*rz3) ./ rexp9 -9*pagemtimes((3*rx3rz + 3*pagemtimes(ry2,rxrz) - 4*rxrz3), rx ./ rexp11)),...
    ((6*pagemtimes(ry,rxrz)) ./rexp9 -9*pagemtimes((3*rx3rz + 3*pagemtimes(ry2,rxrz) - 4*rxrz3),ry ./rexp11)),...
    ((3*rx3 + 3*pagemtimes(ry2,rx) - 12*pagemtimes(rx,rz2)) ./ rexp9 -9*pagemtimes((3*rx3rz + 3*pagemtimes(ry2,rxrz) - 4*rxrz3), rz ./rexp11))];



AJ3y = J3xy_const*[(6*pagemtimes(rxry, rz) ./ rexp9 -9*pagemtimes((3*pagemtimes(rx2ry,rz) + 3*ry3rz - 4*ryrz3),rx ./rexp11)),...
    ((3*rx2rz + 9*ry2rz - 4*rz3) ./rexp9 -9*pagemtimes((3*pagemtimes(rx2ry,rz) + 3*ry3rz - 4*ryrz3),ry ./rexp11)),...
    ((3*rx2ry + 3*ry3 - 12*ryrz2) ./rexp9 -9*pagemtimes((3*pagemtimes(rx2ry,rz) + 3*ry3rz - 4*ryrz3),rz ./rexp11))];


AJ3z = J3z_const*[((12*rx3 + 12*pagemtimes(rx,ry2) - 48*pagemtimes(rx,rz2)) ./ rexp9 - 9*pagemtimes((3*rx4 + 6*rx2ry2 - 24*rx2rz2 + 3*ry4 - 24*ry2rz2 + 8*rz4), rx ./rexp11)),...
    ((12*rx2ry + 12*ry3 - 48*ryrz2) ./ rexp9 - 9*pagemtimes((3*rx4 + 6*rx2ry2 - 24*rx2rz2 + 3*ry4 - 24*ry2rz2 + 8*rz4), ry ./rexp11)),...
    ((-48*rx2rz - 48*ry2rz + 32*rz3) ./ rexp9 - 9*pagemtimes((3*rx4 + 6*rx2ry2 - 24*rx2rz2 + 3*ry4 - 24*ry2rz2 + 8*rz4), rz ./rexp11))];



J2_const = 3/2 * J2 * mu * R^2;
AJ2x = J2_const*[5*rz2 ./rexp7 - 35*pagemtimes(rz2,rx2)./rexp9 - 1 ./ rexp5 + 5*rx2./rexp7,...
    -35*pagemtimes(rz2,rxry) ./ rexp9 + 5*rxry ./rexp7,...
    5*rxrz ./ rexp7 + 10*rxrz ./ rexp7 - 35*rxrz3 ./ rexp9];
AJ2y = J2_const*[5*rxry ./rexp7 - 35*pagemtimes(rz2,rxry) ./rexp9,...
    5*rz2 ./ rexp7 - 35*ry2rz2 ./ rexp9 - 1 ./ rexp5 + 5*ry2./rexp7,...
    5*ryrz ./ rexp7 + 10*ryrz ./ rexp7 - 35*ryrz3 ./ rexp9];
AJ2z = J2_const*[15*rxrz ./ rexp7 - 35*rxrz3 ./rexp9,...
    15*ryrz ./ rexp7 - 35*ryrz3 ./ rexp9,...
    15*rz2 ./ rexp7 - 35*rz4 ./ rexp9 - 3 ./ rexp5 + 15*rz2 ./ rexp7];





Atbp_x = [-mu ./ rexp3 + 3*mu ./ rexp5.*rx2, 3*mu./ rexp5 .* rxry, 3*mu ./ rexp5 .* rxrz];
Atbp_y = [3*mu ./ rexp5 .* rxry, -mu ./ rexp3 + 3*mu ./ rexp5 .* ry2, 3*mu ./ rexp5 .* ryrz];
Atbp_z = [3*mu ./ rexp5 .* rxrz, 3*mu ./ rexp5 .* ryrz, -mu ./ rexp3 + 3*mu ./ rexp5  .* rz2];

% rx = rx(:, 1, 1);
% ry = ry(:, 1, 1);
% rz = rz(:, 1, 1);
% rexp5 = rexp5(:, 1, 1);
% rexp3 = rexp3(:, 1, 1);
% rexp7 = rexp7(:, 1, 1);
% rexp9 = rexp9(:, 1, 1);
% rexp11 = rexp11(:, 1, 1);
% rexp13 = rexp13(:, 1, 1);
% 
% 
% 
% AJ4x_check = [15/8*J4*R^4*mu*((5*rx^4 + 6*rx^2*ry^2 - 36*rx^2*rz^2 + ry^4 -12*ry^2*rz^2 + 8*rz^4)/rexp11 -11*(rx^5 + 2*rx^3*ry^2 - 12*rx^3*rz^2 + rx*ry^4 -12*ry^2*rz^2*rx + 8*rx*rz^4)*rx/rexp13), 15/8*J4*R^4*mu*((4*rx^3*ry + 4*rx*ry^3 - 24*ry*rz^2*rx)/rexp11 -11*(rx^5 + 2*rx^3*ry^2 - 12*rx^3*rz^2 + rx*ry^4 -12*ry^2*rz^2*rx + 8*rx*rz^4)*ry/rexp13),...
%     15/8*J4*R^4*mu*((-24*rx^3*rz - 24*ry^2*rz*rx + 32*rx*rz^3)/rexp11 -11*(rx^5 + 2*rx^3*ry^2 - 12*rx^3*rz^2 + rx*ry^4 -12*ry^2*rz^2*rx + 8*rx*rz^4)*rz/rexp13)];
% AJ4y_check = [15/8*J4*R^4*mu*((4*rx^3*ry + 4*rx*ry^3 - 24*rx*ry*rz^2)/rexp11 - 11*(rx^4*ry + 2*rx^2*ry^3 - 12*rx^2*ry*rz^2 + ry^5 - 12*ry^3*rz^2 + 8*ry*rz^4)*rx/rexp13), 15/8*J4*R^4*mu*((rx^4 + 6*rx^2*ry^2 -12*rx^2*rz^2 + 5*ry^4 - 36*ry^2*rz^2+8*rz^4)/rexp11 - 11*(rx^4*ry + 2*rx^2*ry^3 - 12*rx^2*ry*rz^2 + ry^5 - 12*ry^3*rz^2 + 8*ry*rz^4)*ry/rexp13), 15/8*J4*R^4*mu*((-24*rx^2*ry*rz -24*ry^3*rz + 32*ry*rz^3)/rexp11 - 11*(rx^4*ry + 2*rx^2*ry^3 - 12*rx^2*ry*rz^2 + ry^5 - 12*ry^3*rz^2 + 8*ry*rz^4)*rz/rexp13)];
% AJ4z_check = [5/8*J4*R^4*mu*((60*rx^3*rz + 60*rx*ry^2*rz - 80*rx*rz^3)/rexp11 - 11*(15*rx^4*rz+30*rx^2*ry^2*rz -40*rx^2*rz^3 + 15*ry^4*rz-40*ry^2*rz^3+8*rz^5)*rx/rexp13), 5/8*J4*R^4*mu*((60*rx^2*ry*rz + 60*ry^3*rz -80*ry*rz^3)/rexp11 - 11*(15*rx^4*rz+30*rx^2*ry^2*rz -40*rx^2*rz^3 + 15*ry^4*rz-40*ry^2*rz^3+8*rz^5)*ry/rexp13), 5/8*J4*R^4*mu*((15*rx^4 + 30*rx^2*ry^2 - 120*rx^2*rz^2+15*ry^4-120*ry^2*rz^2+40*rz^4)/rexp11 - 11*(15*rx^4*rz+30*rx^2*ry^2*rz -40*rx^2*rz^3 + 15*ry^4*rz-40*ry^2*rz^3+8*rz^5)*rz/rexp13)];
% 
% 
% 
% 
% 
% AJ3x_check = [-5/2*J3*R^3*mu*((9*rx^2*rz + 3*ry^2*rz - 4*rz^3)/rexp9 -9*(3*rx^3*rz + 3*ry^2*rx*rz - 4*rx*rz^3)*rx/rexp11), -5/2*J3*R^3*mu*((6*ry*rx*rz)/rexp9 -9*(3*rx^3*rz + 3*ry^2*rx*rz - 4*rx*rz^3)*ry/rexp11), -5/2*J3*R^3*mu*((3*rx^3 + 3*ry^2*rx - 12*rx*rz^2)/rexp9 -9*(3*rx^3*rz + 3*ry^2*rx*rz - 4*rx*rz^3)*rz/rexp11)];
% AJ3y_check = [-5/2*J3*R^3*mu*(6*rx*ry*rz/rexp9 -9*(3*rx^2*ry*rz + 3*ry^3*rz - 4*ry*rz^3)*rx/rexp11), -5/2*J3*R^3*mu*((3*rx^2*rz + 9*ry^2*rz - 4*rz^3)/rexp9 -9*(3*rx^2*ry*rz + 3*ry^3*rz - 4*ry*rz^3)*ry/rexp11), -5/2*J3*R^3*mu*((3*rx^2*ry + 3*ry^3 - 12*ry*rz^2)/rexp9 -9*(3*rx^2*ry*rz + 3*ry^3*rz - 4*ry*rz^3)*rz/rexp11)];
% AJ3z_check = [1/2*J3*R^3*mu*((12*rx^3 + 12*rx*ry^2 - 48*rx*rz^2)/rexp9 - 9*(3*rx^4 + 6*rx^2*ry^2 - 24*rx^2*rz^2 + 3*ry^4 - 24*ry^2*rz^2 + 8*rz^4)*rx/rexp11), 1/2*J3*R^3*mu*((12*rx^2*ry + 12*ry^3 - 48*ry*rz^2)/rexp9 - 9*(3*rx^4 + 6*rx^2*ry^2 - 24*rx^2*rz^2 + 3*ry^4 - 24*ry^2*rz^2 + 8*rz^4)*ry/rexp11), 1/2*J3*R^3*mu*((-48*rx^2*rz - 48*ry^2*rz + 32*rz^3)/rexp9 - 9*(3*rx^4 + 6*rx^2*ry^2 - 24*rx^2*rz^2 + 3*ry^4 - 24*ry^2*rz^2 + 8*rz^4)*rz/rexp11)];
% 
% 
% AJ2x_check = 3/2 * J2 * mu * R^2*[5*rz^2/rexp7-35*rz^2*rx^2/rexp9-1/rexp5+5*rx^2/rexp7, -35*rz^2*rx*ry/rexp9 + 5*rx*ry/rexp7,                 5*rx*rz/rexp7 + 10*rx*rz/rexp7 - 35*rz^3*rx/rexp9];
% AJ2y_check = 3/2 * J2 * mu * R^2*[5*ry*rx/rexp7 - 35*rz^2*ry*rx/rexp9,                  5*rz^2/rexp7-35*rz^2*ry^2/rexp9-1/rexp5+5*ry^2/rexp7, 5*ry*rz/rexp7 + 10*ry*rz/rexp7 - 35*rz^3*ry/rexp9];
% AJ2z_check = 3/2 * J2 * mu * R^2*[15*rz*rx/rexp7 - 35*rz^3*rx/rexp9,                    15*rz*ry/rexp7 - 35*rz^3*ry/rexp9,                    15*rz^2/rexp7 - 35*rz^4/rexp9-3/rexp5+15*rz^2/rexp7];
% 
% 
% 
% Atbp_x_check = [-mu/rexp3 + 3*mu/rexp5*rx^2, 3*mu/rexp5*rx*ry, 3*mu/rexp5*rx*rz];
% Atbp_y_check = [3*mu/rexp5*rx*ry, -mu/rexp3+3*mu/rexp5*ry^2,3*mu/rexp5*ry*rz];
% Atbp_z_check = [3*mu/rexp5*rx*rz, 3*mu/rexp5*ry*rz, -mu/rexp3+3*mu/rexp5*rz^2];
% 
% 
% for i = 1:3
% fprintf('Difference J4x:  %0.3e, %0.3e, %0.3e\n', AJ4x(:, :, i) - AJ4x_check)
% fprintf('Difference J4y:  %0.3e, %0.3e, %0.3e\n', AJ4y(:, :, i) - AJ4y_check)
% fprintf('Difference J4z:  %0.3e, %0.3e, %0.3e\n', AJ4z(:, :, i) - AJ4z_check)
% 
% fprintf('Difference J3x:  %0.3e, %0.3e, %0.3e\n', AJ3x(:, :, i) - AJ3x_check)
% fprintf('Difference J3y:  %0.3e, %0.3e, %0.3e\n', AJ3y(:, :, i) - AJ3y_check)
% fprintf('Difference J3z:  %0.3e, %0.3e, %0.3e\n', AJ3z(:, :, i) - AJ3z_check)
% 
% fprintf('Difference J2x:  %0.3e, %0.3e, %0.3e\n', AJ2x(:, :, i) - AJ2x_check)
% fprintf('Difference J2y:  %0.3e, %0.3e, %0.3e\n', AJ2y(:, :, i) - AJ2y_check)
% fprintf('Difference J2z:  %0.3e, %0.3e, %0.3e\n', AJ2z(:, :, i) - AJ2z_check)
% 
% fprintf('Difference TBPx: %0.3e, %0.3e, %0.3e\n', Atbp_x(:, :, i) - Atbp_x_check)
% fprintf('Difference TBPy: %0.3e, %0.3e, %0.3e\n', Atbp_y(:, :, i) - Atbp_y_check)
% fprintf('Difference TBPz: %0.3e, %0.3e, %0.3e\n\n', Atbp_z(:, :, i) - Atbp_z_check)
% end



AJ4 = [AJ4x; AJ4y; AJ4z];
AJ3 = [AJ3x; AJ3y; AJ3z];
AJ2 = [AJ2x; AJ2y; AJ2z];
ATBP = [Atbp_x; Atbp_y; Atbp_z];

A = AJ4 + AJ3 + AJ2 + ATBP;

end