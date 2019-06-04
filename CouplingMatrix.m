function C = CouplingMatrix(txSpacing,fc, numAnt)
%txSpacing = 0.1;
%rxSpacing = 0.1;
Nt=numAnt;
lambda = physconst('lightspeed')/fc;
antElement = dipole( ...
    'Length', lambda/2, ...
    'Width',  lambda/100);
txArray = linearArray( ...
    'Element',        antElement,...
    'NumElements',    Nt,...
    'ElementSpacing', txSpacing*lambda);      
arrayModel=txArray;
% Calculate impedance matrix
S = sparameters(arrayModel, fc);
Ztx = s2z(squeeze(S.Parameters));

% Enforce symmetry on the impedance matrix
Zu = triu(Ztx);
Zl = Zu.';
Zl([1,4]) = 0;
Ztx = Zl + Zu;

% Form coupling matrix as per Eq. (6) in paper
Zload = Ztx(1,1)';
Zlm = Zload .* eye(prod(numAnt));
C = (Zload + Ztx(1,1)) .* inv((Ztx + Zlm));

% [EOF]