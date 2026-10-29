function [dX_nonDimR_dX_dimN] = dX_nonDimRot_dX_dimN(X_N, t, gravParam)
dX_nonDimR_dX_dimN = zeros(9);
L = gravParam(2);
n = gravParam(3);
tCharacteristic = gravParam(4);

NR = nonDimRotational_dimInertial_conversion(X_N, t, L, n, tCharacteristic, 'DCM_R2N');
RN = NR';

drRdrN = RN/L; % good 

dvRdrN = RN*[0, 1, 0; -1, 0, 0; 0, 0, 0]/L; % good 

dvRdvN = RN/L/n; 

daRdaN = RN/(L)*tCharacteristic^2;

dX_nonDimR_dX_dimN(1:3, 1:3) = drRdrN;
dX_nonDimR_dX_dimN(4:6, 1:3) = dvRdrN;
dX_nonDimR_dX_dimN(4:6, 4:6) = dvRdvN;
dX_nonDimR_dX_dimN(7:9, 7:9) = daRdaN;


end