syms rx ry rz vx vy vz mu mu1 mu2 L n t
r1=sqrt((rx+mu)^2+ry^2+rz^2);
r2=sqrt((rx-1+mu)^2+ry^2+rz^2);
c1=(1-mu)/r1^3;
c2=mu/r2^3;
xvec=[rx; ry; rz; vx; vy; vz];
f=[vx; vy; vz; 2*vy+rx-c1*(rx+mu)-c2*(rx-1+mu); -2*vx+ry-c1*ry-c2*ry; -c1*rz-c2*rz];
fx=jacobian(f,xvec);
fxx=diff(fx,rx);
fxx(:,:,2)=diff(fx,ry);
fxx(:,:,3)=diff(fx,rz);
fxx(:,:,4)=diff(fx,vx);
fxx(:,:,5)=diff(fx,vy);
fxx(:,:,6)=diff(fx,vz);



% r1expneg3 = sqrt((rx+mu)^2+ry^2+rz^2)^-3;
% r1expneg5 = sqrt((rx+mu)^2+ry^2+rz^2)^-5;


% 
% r2expneg3 = sqrt((rx-1+mu)^2+ry^2+rz^2)^(-3);
% r2expneg5 = sqrt((rx-1+mu)^2+ry^2+rz^2)^(-5);
% 
% xmu = rx+mu;
% x1mu = rx-1+mu;
% 
% 
% A(1, 1) = 1 - (1-mu)*r1expneg3 + 3*(1-mu)*(xmu)^2*r1expneg5 - mu*r2expneg3 + 3*mu*(x1mu)^2 *r2expneg5;
% A(1, 2) = 3*(1-mu)*(xmu)*ry*r1expneg5 + 3*ry*mu*(x1mu)*r2expneg5;
% A(1, 3) = 3*(1-mu)*(xmu)*rz*r1expneg5 + 3*rz*mu*(x1mu)*r2expneg5;
% A(2, 1) = A(1, 2);
% A(2, 2) = 1 - (1-mu)*r1expneg3 + 3*(1-mu)*ry^2*r1expneg5 - mu*r2expneg3 + 3*mu*ry^2*r2expneg5;
% A(2, 3) = 3*(1-mu)*ry*rz*r1expneg5 + 3*mu*ry*rz*r2expneg5;
% A(3, 1) = A(1, 3);
% A(3, 2) = A(2, 3);
% A(3, 3) = -(1-mu)*r1expneg3 + 3*(1-mu)*rz^2*r1expneg5 - mu*r2expneg3 + 3*mu*rz^2*r2expneg5;



% [dadr] = gravJac_CRTBP([rx; ry; rz], mu);
% Ac(4:6, 1:3) = dadr;
% Ac(1:3, 4:6) = eye(3);
% Ac(4:6, 4:6) = [0, 2, 0; -2, 0, 0; 0, 0, 0];



%% Checking 3BP with inetial info
r = [rx; ry; rz];
R1 = [-mu*L*cos(n*t); -mu*L*sin(n*t); 0];
R2 = [(1-mu)*cos(n*t); (1-mu)*sin(n*t); 0];

r1 = r - R1;
r1Mag = sqrt(r1(1)^2+r1(2)^2+r1(3)^2);
r2 = r - R2;
r2Mag = sqrt(r2(1)^2+r2(2)^2+r2(3)^2);

r2 = r - R2;

F1 = -mu1/(r1Mag)^3*r1;
F2 = -mu2/(r2Mag)^3*r2;


f_a = F1 + F2;
ax = -mu1*r1(1)/r1Mag^3 - mu2*r2(1)/r2Mag^3;
A = jacobian(f_a, r);
daxdrx = -mu1/ (r1Mag^3) + mu1*(3/(r1Mag^5) * r1(1)^2) - mu2/ (r2Mag^3) + mu2*(3*r2Mag^-5 * r2(1)^2);
daxdry = mu1*(3/(r1Mag^5) * r1(1)*r1(2)) + mu2*(3*r2Mag^-5 * r2(1)*r2(2));
daxdrz = mu1*(3/(r1Mag^5) * r1(1)*r1(3)) + mu2*(3*r2Mag^-5 * r2(1)*r2(3));
daydrx = mu1*(3/(r1Mag^5) * r1(1)*r1(2)) + mu2*(3*r2Mag^-5 * r2(1)*r2(2));
daydry = -mu1/ (r1Mag^3) + mu1*(3/(r1Mag^5) * r1(2)^2) - mu2/ (r2Mag^3) + mu2*(3*r2Mag^-5 * r2(2)^2);
daydrz = mu1*(3/(r1Mag^5) * r1(2)*r1(3)) + mu2*(3*r2Mag^-5 * r2(2)*r2(3));
dazdrz = -mu1/ (r1Mag^3) + mu1*(3/(r1Mag^5) * r1(3)^2) - mu2/ (r2Mag^3) + mu2*(3*r2Mag^-5 * r2(3)^2);


% simplify(daxdrx - A(1, 1))
% simplify(daxdry - A(1, 2))
% simplify(daxdrz - A(1, 3))
% simplify(daydrx - A(2, 1))
% simplify(daydry - A(2, 2))
% simplify(daydrz - A(2, 3))
% simplify(dazdrz - A(3, 3))
