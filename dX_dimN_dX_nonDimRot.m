function [dX_dimN_dX_nonDimR] = dX_dimN_dX_nonDimRot(X_N, t, gravParam)
dX_dimN_dX_nonDimR = zeros(9);
L = gravParam(2);
n = gravParam(3);
tCharacteristic = gravParam(4);

NR = nonDimRotational_dimInertial_conversion(X_N, t, L, n, tCharacteristic, 'DCM_R2N');

drNdrR = NR*L;

dvNdrR = [-n*NR(2, :); n*NR(1, :); zeros(1, 3)]*L;

dvNdvR = NR*L*n;

daNdaR = NR*(L)/tCharacteristic/tCharacteristic;

dX_dimN_dX_nonDimR(1:3, 1:3) = drNdrR;
dX_dimN_dX_nonDimR(4:6, 1:3) = dvNdrR;
dX_dimN_dX_nonDimR(4:6, 4:6) = dvNdvR;
dX_dimN_dX_nonDimR(7:9, 7:9) = daNdaR;


end