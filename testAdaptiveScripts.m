%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking my algorithms versus Rahul's 
%
% Author: Rachel Mamich
% Date Created: 04/05/23
% Date Modified: 04/06/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% House Keeping
clear all
close all

tau = 1e4;
Qk = [3, 1; 1, 2];
Rk = 2;
x0 = [0; 0];
n = length(x0);
nq = size(Qk,1);
nr = size(Rk,1);
mq = nq*(nq+1)/2;
mr = nr*(nr+1)/2;
p = size(Rk,1);
P0 = 10*eye(n);

% Number of measurements stacked
m = 3;

l = 1; % Number of control inputs


% Loading Data
tf = .02;
dt = .0001;
N = tf/dt;
n_Monte	= 1;
rng(32)

% Initialization
xhatk = x0 + chol(P0)'*randn(n,1);
xhatkb = xhatk;
xhatkm = xhatk;
xk = x0;
Pk = P0;
Pkb = P0;
Pkm = P0;
PQRhatk = 1*eye(mq+mr);
PQRhatkm = 1*eye(mq+mr);
% PQhatkm = 1*eye(nq); PRhatkm = 1*eye(nr);
RW = 1;
Rm = RW;
QW = .1;
Qm = QW*eye(nq+nr);
theta_true = [vech(Qk); vech(Rk)];
lambda_min = 0.1;		% lower bound on the eigenvalue of Qk and Rk
Qhatk = Qk; Qhatk(1) = 1;
Rhatk = Rk; Rhatk = 1;
Q0 = Qhatk; R0 = Rhatk;
Qhatkm = Qhatk;
Rhatkm = Rhatk;
Qhatkr = Qhatk;
Rhatkr = Rhatk;
Fkwind = zeros(n, n, m);
Hkwind = zeros(p, n, m+1);
Gkwind = zeros(n, l, m);
Ziwind = zeros(n, m+1);
Awind = zeros(n, p, m+1, m+1);		% Coefficients hist
Bwind = zeros(n, n, m, m);			% Coefficients hist
ukwind = zeros(1,m+1);
nmeas = 0;
mu = 0;


reverseStr = '';

% Filter
for k =1:N
% Calculate system matrices
[Fk, Gk, Hk] = systemMatrices(k*dt);
uk = randn(size(Gk, 2), 1);
% uk = 0;
ukwind = [uk, ukwind(:,1:end-1)];

% Find true state and generate measurements
xkp1 = Fk*xk + Gk*uk + chol(Qk)*randn(n,1);
ykp1 = Hk*xkp1 + chol(Rk)*randn(p,1);
xhist(:,k) = xkp1;
yhist(:,k) = ykp1;

% Kalman Update using estimated
[xhatk, Pk, Pbark] = kalmanPropUpdate(xhatk, Pk, Fk, Hk, Gk, uk, ykp1, Qhatk, Rhatk);

% Baseline Kalman Filter
[xhatkb, Pkb, Pbarb] = kalmanPropUpdate(xhatkb, Pkb, Fk, Hk, Gk, uk, ykp1, Qk, Rk);

% Saving time history
xhathist(:,k) = xhatk;
Pbarhist(:,:,k) = Pbark;

xhathistb(:,k) = xhatkb;
Pbarhistb(:,:,k) = Pbarb;

xhathistm(:,k) = xhatkm;
% Pbarhistm(:,:,k) = Pbarkm;

Fkwind = cat(3, Fk, Fkwind(:,:,1:m-1));
Hkwind = cat(3, Hk, Hkwind(:,:,1:m));
Gkwind = cat(3, Gk, Gkwind(:,:,1:m-1));

% Adaptive Filter
if k > m
	[A, B] = calcTimeSeriesCoeffs(Fkwind, Gkwind, Hkwind, m);
    [ARachel, BRachel] = constructAB(Fkwind, Hkwind, m);
    
    for qq = 1:m+1
        fprintf('ARahul(:, :, %d) - ARachel(:, :, %d) = [%0.6e, %0.6e] \n', qq, qq, ARachel(:,:,qq) - A(:, :, qq))
    end
    
    for qq = 1:m
        fprintf('BRahul(:, :, %d) - BRachel(:, :, %d) = [%0.6e, %0.6e; %0.6e, %0.6e] \n', qq, qq, BRachel(:,:,qq) - B(:, :, qq))
    end
    
	Awind(:,:,:,2:end) = Awind(:,:,:,1:end-1);
	Bwind(:,:,:,2:end) = Bwind(:,:,:,1:end-1);
	Awind(:,:,:,1) = A;
	Bwind(:,:,:,1) = B;
	% Yk 		= [yhist(:,k-1);yhist(:,k-2)];
	% Ykp1 	= [yhist(:,k);yhist(:,k-1)];
end

if k > m + 3
	% nmeas = nmeas + 1;
	% cov_Zi 	= (k-m-1)*cov_Zi/(k-m) + Zi*Zi'/(k-m);
	% if nmeas > 1
	%     cov_Zi 	= (nmeas-2)*cov_Zi/(nmeas-1) + (Zi*Zi' - Zi*mu' - mu*Zi' + mu*mu')/nmeas;
	% else
	%     cov_Zi = zeros(l);
	% end
	% mu		= (nmeas-1)*mu/nmeas + Zi/nmeas;
	Zik = A(:,:,1)*yhist(:,k)+A(:,:,2)*yhist(:,k-1)+A(:,:,3)*yhist(:,k-2)+A(:,:,4)*yhist(:,k-3)...
		- B(:,:,1)*Gkwind(:,:,1)*ukwind(:,1) - B(:,:,2)*Gkwind(:,:,2)*ukwind(:,2)...
		- B(:,:,3)*Gkwind(:,:,3)*ukwind(:,3);
%     ZikParts = [A(:,:,1)*yhist(:,k), A(:,:,2)*yhist(:,k-1), A(:,:,3)*yhist(:,k-2), A(:,:,4)*yhist(:,k-3), - B(:,:,1)*Gkwind(:,:,1)*ukwind(:,1), - B(:,:,2)*Gkwind(:,:,2)*ukwind(:,2), - B(:,:,3)*Gkwind(:,:,3)*ukwind(:,3)]
    [ZRachel] = meask_adaptivePSD(A, B, flip(yhist(:,k-3:k)), ukwind, Gkwind);
    
	Zihist(:,k-m) = Zik;
	Ziwind = [Zik, Ziwind(:,1:end-1)];
	HH = formMeasMatrix(Awind, Bwind);
    [HRachel, ABlock] = formMeasMatrix_variableBufferSize(Awind, Bwind, m);
    
    fprintf('H matrix difference\n')
    HH - [HRachel, ABlock]
    
	bb = [kronh(Ziwind(:,1)); kron(Ziwind(:,2), Ziwind(:,1));...
		kron(Ziwind(:,3), Ziwind(:,1)); kron(Ziwind(:,4), Ziwind(:,1))];
    [bRachel] = Rachel_b(m, n, Ziwind);
    
    fprintf('b difference\n')
    bb - bRachel
	% rank(HH), cond(HH)
	% pause
	[thetahat, PQRhatk] = recursiveLS(HH, bb, Qhatk, Rhatk, PQRhatk, QW, RW);
	% [thetahatm] = riemUpdate(HH, bb, Qhatkm, Rhatkm, PQhatkm, PRhatkm, Qm, Rm, k);
% 	[Qhatkm, Rhatkm, PQRhatkm] = riemOptUpdate(HH, bb, p1, Qhatkm, Rhatkm, PQRhatkm, Q0, R0, QW, RW, lambda_min, alpha, beta, k);
	Qhatk = [thetahat(1), thetahat(2); thetahat(2), thetahat(3)];
	Rhatk = thetahat(4);
	% Qhatkm = [thetahatm(1), thetahatm(2); thetahatm(2), thetahatm(3)];
	% Rhatkm = thetahatm(4);
	% eigQ = eig(Qhatkm)
	% pause
end
Qhist(:,:,k) = Qhatk;
Rhist(:,:,k) = Rhatk;
% Qhistm(:,:,k) = Qhatkm;
Rhistm(:,:,k) = Rhatkm;
% Qhistr(:,:,k) = Qhatkr;
% Rhistr(:,:,k) = Rhatkr;

% Update states
xk = xkp1;

% Progress bar
% percentDone = 100 * k / N;
% msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
% fprintf([reverseStr, msg]);
% reverseStr = repmat(sprintf('\b'), 1, length(msg));

end

% Calc eigenvalues of Qhat, Rhat
for i=1:size(Qhist,3)
	lambdaQhist(:,i) = sort(eig(Qhist(:,:,i)));
	lambdaQhistm(:,i) = sort(eig(Qhistm(:,:,i)));
	lambdaPhist(:,i) = sort(eig(Pbarhist(:,:,i)));
	lambdaPhistb(:,i) = sort(eig(Pbarhistb(:,:,i)));
	lambdaPhistm(:,i) = sort(eig(Pbarhistm(:,:,i)));
	Pnormhist(i) = norm(Pbarhist(:,:,i) - Pbarhistb(:,:,i), 'fro');
	Pnormhistm(i) = norm(Pbarhistm(:,:,i) - Pbarhistb(:,:,i), 'fro');
	Qnormhist(i) = norm(Qhist(:,:,i) - Qk, 'fro');
	Qnormhistm(i) = norm(Qhistm(:,:,i) - Qk, 'fro');
	Rnormhist(i) = norm(Rhist(:,:,i) - Rk, 'fro');
	Rnormhistm(i) = norm(Rhistm(:,:,i) - Rk, 'fro');
end
lambdaQtrue = sort(eig(Qk));

fprintf('\n');
thist = dt*(1:N);
plotting_ltv




function [Fk, Gk, Hk] = systemMatrices(t)
	ck = -0.7 + 0.2*cos(2*pi*t);
	ak = ck + i*(.4+.2*sin(2*pi*t));
	bk = ck - i*(.4+.2*sin(2*pi*t));
	dk = 2*sin(10*pi*t);
	Fk = [0, 1; -ak*bk, -(ak + bk)];
	Gk = [1; 1];
	Hk = [1, dk];

	% % LTI system
	% Fk = diag([1, 2]); Fk(1,2) = 2;
	% Hk = [1, 0];
end




% kalmanPropUpdate: One step propagation and update of the linear Kalman filter
function [xhatkp1, Pkp1, Pbarkp1] = kalmanPropUpdate(xhatk, Pk, Fk, Hk, Gk, uk, ykp1, Qk, Rk)

n = length(xhatk);
p = length(ykp1);
xbark = Fk*xhatk;
Pbarkp1 = Fk*Pk*Fk' + Qk;
Wk = Pbarkp1*Hk'*((Hk*Pbarkp1*Hk' + Rk)\eye(p));
Pkp1 = (eye(n)-Wk*Hk)*Pbarkp1*(eye(n)-Wk*Hk)' + Wk*Rk*Wk';
nu = ykp1 - Hk*xbark;
xhatkp1 = xbark + Wk*nu;

end



% calcTimeSeriesCoeffs: Calculates the time series coefficients at time k using
% 						Fk, Gk, Hk and m (number of measurements stacked)
function [A, B] = calcTimeSeriesCoeffs(Fkwind, Gkwind, Hkwind, m)

% Fkwind is n x n x m
% Hkwind is p x n x m+1
n = size(Fkwind, 1);
p = size(Hkwind, 1);
A = zeros(n, p, m+1);
B = zeros(n, n, m);

Mk_km1 = zeros(n);
Mkm1_km = zeros(n);

% TODO: do it for a general m and not just m = 3
Mk_km1 = Hkwind(:,:,3)'*Hkwind(:,:,3) + Fkwind(:,:,2)'*Hkwind(:,:,2)'*Hkwind(:,:,2)*Fkwind(:,:,2) +...
		Fkwind(:,:,2)'*Fkwind(:,:,1)'*Hkwind(:,:,1)'*Hkwind(:,:,1)*Fkwind(:,:,1)*Fkwind(:,:,2);
Mkm1_km = Hkwind(:,:,4)'*Hkwind(:,:,4) + Fkwind(:,:,3)'*Hkwind(:,:,3)'*Hkwind(:,:,3)*Fkwind(:,:,3) +...
		Fkwind(:,:,3)'*Fkwind(:,:,2)'*Hkwind(:,:,2)'*Hkwind(:,:,2)*Fkwind(:,:,2)*Fkwind(:,:,3);
% Mkm1_km = Hkwind(:,:,4)'*Hkwind(:,:,4) + Fkwind(:,:,3)'*Hkwind(:,:,3)'*Hkwind(:,:,3)*Fkwind(:,:,3) + Fkwind(:,:,2)'*Hkwind(:,:,2)'*Hkwind(:,:,2)*Fkwind(:,:,2);
phi_k_km1 = Fkwind(:,:,1)*Fkwind(:,:,2);
phi_km1_km = Fkwind(:,:,2)*Fkwind(:,:,3);
A(:,:,1) = Mk_km1\(phi_k_km1'*Hkwind(:,:,1)');
A(:,:,2) = (Mk_km1\eye(n) - Fkwind(:,:,m)*(Mkm1_km\(Fkwind(:,:,m)')))*(Fkwind(:,:,2)'*Hkwind(:,:,2)');
A(:,:,3) = (Mk_km1\eye(n) - Fkwind(:,:,m)*(Mkm1_km\(Fkwind(:,:,m)')))*Hkwind(:,:,3)';
A(:,:,4) = -Fkwind(:,:,m)*(Mkm1_km\Fkwind(:,:,m)')*Hkwind(:,:,m+1)';

B(:,:,1) = Mk_km1\(phi_k_km1'*Hkwind(:,:,1)'*Hkwind(:,:,1));
Mk_m1 = Hkwind(:,:,2)'*Hkwind(:,:,2) + Fkwind(:,:,1)'*Hkwind(:,:,1)'*Hkwind(:,:,1)*Fkwind(:,:,1);
Mkm1_m2 = Hkwind(:,:,3)'*Hkwind(:,:,3) + Fkwind(:,:,2)'*Hkwind(:,:,2)'*Hkwind(:,:,2)*Fkwind(:,:,2);
B(:,:,2) = Mk_km1\(Fkwind(:,:,2)'*Mk_m1) - Fkwind(:,:,m)*(Mkm1_km\(Fkwind(:,:,3)'*Fkwind(:,:,2)'))*Hkwind(:,:,2)'*Hkwind(:,:,2);
B(:,:,3) = eye(n) - Fkwind(:,:,m)*(Mkm1_km\(Fkwind(:,:,m)'))*Mkm1_m2;

end



function HH = formMeasMatrix(Awind, Bwind)

AAkk = 0; BBkk = 0; AAkkm1 = 0; BBkkm1 = 0; BBkkm2 = 0; AAkkm2 = 0; AAkkm3 = 0;

for i=1:size(Bwind,3)
	BBkk = BBkk + kron(Bwind(:,:,i,1), Bwind(:,:,i,1));
end
for i=1:size(Awind,3)
	AAkk = AAkk + kron(Awind(:,:,i,1), Awind(:,:,i,1));
end
% Remove columns corr. to symmetric elements 
BBkk = [BBkk(:,1), BBkk(:,2)+BBkk(:,3), BBkk(:,4)];
% Remove rows corr. to symmetric elements as well
BBkk = [BBkk(1,:); BBkk(2,:); BBkk(4,:)];
AAkk = [AAkk(1,:); AAkk(2,:); AAkk(4,:)];

for i=1:size(Bwind,3)-1
	BBkkm1 = BBkkm1 + kron(Bwind(:,:,i+1,1), Bwind(:,:,i,2));
end
for i=1:size(Awind,3)-1
	AAkkm1 = AAkkm1 + kron(Awind(:,:,i+1,1), Awind(:,:,i,2));
end
% Remove columns corr. to symmetric elements 
BBkkm1 = [BBkkm1(:,1), BBkkm1(:,2)+BBkkm1(:,3), BBkkm1(:,4)];

for i=1:size(Bwind,3)-2
	BBkkm2 = BBkkm2 + kron(Bwind(:,:,i+2,1), Bwind(:,:,i,3));
end
for i=1:size(Awind,3)-2
	AAkkm2 = AAkkm2 + kron(Awind(:,:,i+2,1), Awind(:,:,i,3));
end
BBkkm2 = [BBkkm2(:,1), BBkkm2(:,2)+BBkkm2(:,3), BBkkm2(:,4)];

for i=1:size(Awind,3)-3
	AAkkm3 = AAkkm3 + kron(Awind(:,:,i+3,1), Awind(:,:,i,4));
end

% BBkk, AAkk, BBkkm1, AAkkm1, AAkkm2
n = size(Bwind, 1);
m = n*(n+1)/2;
HH = [BBkk, AAkk;...
      BBkkm1, AAkkm1;...
      BBkkm2, AAkkm2;...
      zeros(n*n,m), AAkkm3];
% rank(HH), size(HH), cond(HH)
% pause

end

function [b] = Rachel_b(bufferSize, nA, ZWindow)
endFirst = nA*(nA + 1)/2;
b = zeros((bufferSize)*nA*nA+ nA*(nA+1)/2, 1);
b(1:endFirst) = kronh(ZWindow(:, 1));
for i = 2:(bufferSize + 1)
    b((nA*nA*(i-2)+1+ endFirst):(nA*nA*(i-1) + endFirst)) = kron(ZWindow(:, i), ZWindow(:, 1));
end
end

% recursiveLS: Returns the updated estimate using the recursive least squares
function [thetahatkp1, Phatkp1] = recursiveLS(HH, bb, Qhatk, Rhatk, Phatk, QW, RW)
% RW is the variance in each of the measurements

thetahat = [uvec(Qhatk); uvec(Rhatk)];
n = length(thetahat);

% Phatk = Phatk + QW*eye(n);
p = length(bb);

Lk = Phatk*HH'*((RW*eye(p) + HH*Phatk*HH')\eye(p));

thetahatkp1 = thetahat + Lk*(bb - HH*thetahat);
Phatkp1 = (eye(n) - Lk*HH)*Phatk*(eye(n) - Lk*HH)' + Lk*Lk';

% (HH*thetahatkp1 - bb)'*(HH*thetahatkp1 - bb)/RW + (thetahatkp1 - thetahat)'*(Phatk\(thetahatkp1 - thetahat))
% Lk*(bb - HH*thetahat)

end
