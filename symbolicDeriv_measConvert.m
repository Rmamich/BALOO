%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symbolic deriv for transforming meas uncertainty
% 
%
% Created: 08/16/23
% Last Modified: 08/16/23
% Author: Rachel Mamich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all



% sigmarho = sym('sigmarho');
% sigmaA = sym('sigmaA');
% sigmaB = sym('sigmaB');

sigmarho = 1000;
sigmaA = 32*pi/180;
sigmaB = 12*pi/180;

sigmaz = sqrt(sigmarho^2*cos(sigmaA)^2 / (1 + cos(sigmaA)^2*tan(sigmaB)^2));
sigmax = tan(sigmaA)*sigmaz;
sigmay = tan(sigmaB)*sigmaz;

denomsigmaZ = (1 + cos(sigmaA)^2*tan(sigmaB)^2);
rho2 = sqrt(sigmax^2 + sigmay^2 + sigmaz^2);
A2 = atan2(sigmax, sigmaz)*180/pi;
B2 = atan2(sigmay, sigmaz)*180/pi;

dzdrho = diff(sigmaz, sigmarho); % sigmarho*cos(sigmaA)^2/denomsigmaZ/sigmaz
dzdA = diff(sigmaz, sigmaA); % sigmarho^2/sigmaz*cos(sigmaA)*sin(sigmaA)*(cos(sigmaA)^2/denomsigmaZ^2*tan(sigmaB)^2 - 1/denomsigmaZ)
dzdB = diff(sigmaz, sigmaB); % -sigmarho^2/denomsigmaZ^2*tan(sigmaB)*cos(sigmaA)^4*sec(sigmaB)^2/sigmaz


dxdrho = simplify(diff(sigmax, sigmarho)); % tan(sigmaA)*dzdrho
dxdA = simplify(diff(sigmax, sigmaA)); % dzdA*tan(sigmaA) + sigmaz*sec(sigmaA)^2
dxdB = simplify(diff(sigmax, sigmaB)); % dzdB*tan(sigmaA)

dydrho = diff(sigmay, sigmarho); % tan(sigmaB)*dzdrho
dydA = diff(sigmay, sigmaA); % tan(sigmaB)*dzdA
dydB = diff(sigmay, sigmaB); % tan(sigmaB)*dzdB + sec(sigmaB)^2*sigmaz











