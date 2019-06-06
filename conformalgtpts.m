fc=2.4e9;%Nt=8;
lambda = physconst('lightspeed')/fc;
%fig=figure;
objectwidthinm=0.15;
objectheightinm=0.075;
xx=[0   objectheightinm    objectheightinm      0   0]*lambda;
yy=[0   0  objectwidthinm objectwidthinm 0]*lambda;
fig=figure;
%imshow('E:\Picture1.png')

plot(xx,yy,'r','linewidth',2) 
xlim([-0.01 0.02])
ylim([0 0.02])
hold on
[Xpos, Ypos]=getpts(fig);
    Nt=length(Xpos);

tota=zeros(Nt,3);
for ii=1:Nt
    Zpos(ii)=0;
arrpos=[Xpos(ii) Ypos(ii) Zpos(ii)];
tota(ii,:)=arrpos;
end
%totalpos=[xpos;ypos;zpos;tpos;rpos;spos;qpos]


antElement = dipole('Length', lambda/2,'Width',  lambda/100);
    txArray = conformalArray('Element',antElement,...
   'ElementPosition',tota);
numAnt=[1 Nt];
show(txArray) 
figure
layout(txArray)
title('conformal array')
grid on
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
C = (Zload + Ztx(1,1)) .* inv((Ztx + Zlm));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%capacity_analysis_%%%%%%%%%
txcoupmat=C;
txCorrMtx=eye(Nt);
txMCCorrMtx = txcoupmat * txCorrMtx * txcoupmat';
% number of channel realization
It = 10000;

% SNR range in dB
SNRdBvalues = -15:30;

% initialize variables: ergodic capacity for each value of SNRdBvalues
Cmimo  = zeros(1, length(SNRdBvalues));
CmimoMC  = zeros(1, length(SNRdBvalues));
Mt=Nt;
Mr=Mt;

for kk=1:It
        SNRidx = 0;  % SNR index
        % generate channel realization
    Hmimo = ( randn(Mr,Mt) + 1i*randn(Mr,Mt) )/sqrt(2);  % mimo
       for SNRdB = SNRdBvalues
                SNR = 10.^(SNRdB./10);  % linear scale
                SNRidx = SNRidx + 1;        
        % MIMO
        Cmimo(SNRidx) = Cmimo(SNRidx) + log2(real(det( eye(Mr) + SNR/Mt*Hmimo*Hmimo' )));
          CmimoMC(SNRidx) = CmimoMC(SNRidx) + log2(real(det( eye(Mr) + SNR/Mt*Hmimo*txMCCorrMtx*Hmimo' )));             
        end
end

% Compute average over all channel realizations

Cmimo = Cmimo/It;
CmimoMC= CmimoMC/It;

% plot
figure(1)
plot(SNRdBvalues, Cmimo,'b','linewidth',2)
hold on
plot(SNRdBvalues, CmimoMC,'r--','linewidth',2)
hold on
grid on
xlabel('Average SNR \rho (dB)')
ylabel('Ergodic Capacity (bits/transmission)')
title('Ergodic Capacity for Rayleigh fast fading channel')
legend('Without Coupling Effect','With Coupling Effect')