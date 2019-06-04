%couplingmatrix
clear;
fc =2.4e9;   
nt=[4];
%sizeofrect=[2 2]
for jj=1:length(nt)
Nt = nt(jj);
if Nt==4
    sizeofrect=[Nt/2 Nt/2];
else
    sizeofrect=[Nt/4 Nt/2];
end
stepsize=0.01;
txSpacing1=0.1:stepsize:1.0;
%    txSpacing=0.5;
lambda = physconst('lightspeed')/fc;
%radiu=txSpacing/(2*(sin (pi/Nt)));
C=zeros(Nt);

for ii=1:length(txSpacing1)
    txSpacing=txSpacing1(ii);
%radiu=txSpacing/(2*(sin (pi/Nt))); %n must be >3
    antElement = dipole('Length', lambda/2,'Width',  lambda/100);
    txArray = rectangularArray('Element',antElement,'size',sizeofrect,'RowSpacing',...
    txSpacing*lambda,'ColumnSpacing',txSpacing*lambda);

numAnt=[1 Nt];

% Calculate impedance matrix
S = sparameters(txArray, fc);
Ztx = s2z(squeeze(S.Parameters));

% Enforce symmetry on the impedance matrix
Zu = triu(Ztx);
Zl = Zu.';
Zl([1,4]) = 0;
Ztx = Zl + Zu;

% Form coupling matrix as per Eq. (6) in paper
Zload = Ztx(1,1)';
Zlm = Zload .* eye(prod(numAnt));
%C = (Zload + Ztx(1,1)) .* inv((Ztx + Zlm));
C(:,:,ii) = (Zload + Ztx(1,1)) .* inv((Ztx + Zlm));
end

normffcc=zeros(length(txSpacing1));
%figure;
for kk=1:length(txSpacing1)
    normffcc1(kk)=norm(C(:,:,kk));
   % plot(txSpacing1(jj),normffcc)
   % hold on
end
plot(txSpacing1,normffcc1,'LineWidth',2)
hold on
end

title('dipole circular array 2.4GHz AE l=\lambda/2& w=\lambda/100')
grid on
xlabel('Tx separation (in \lambda)')
ylabel('Norm of Coupling Matrix')
hold on
 
% figure;
% show(txArray)
% figure;
% layout(txArray)
%reset(C)