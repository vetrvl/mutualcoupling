%clear 
fc=2.4e9;
% number of channel realization
It = 10000;

% SNR range in dB
SNRdBvalues = -15:3:30;

% range of values for each SNR to estimate CDF
%xi_range = [1 2 4 8 10 12 14 18 20 25 30 33 35 40 42 45];
xi_range = [25 30 33 35 40 42 45 75 90 123 145 245 260 285 300 310];
% initialize variables
Cmimo = zeros(1,It);
CmimoMC=zeros(1,It);
Numt=[2 4 8 16];
init=0;
for jk=1:length(Numt)
    Nt=Numt(jk);
    init=init+1;
SNRidx = 0;
numAnt=Nt;
txCorrMtx = eye(Nt);
txcoupmat=CouplingMatrix(0.5,fc, numAnt);
txMCCorrMtx = txcoupmat * txCorrMtx * txcoupmat';
for SNRdB = SNRdBvalues
    
    SNRdB
    SNR = 10.^(SNRdB./10);  % linear scale
    SNRidx = SNRidx + 1;
    
    % collect realizations of the maximal achievable rate
    for kk=1:It
        
        Hmimo = ( randn(Nt) + 1i*randn(Nt) )/sqrt(2);
        Cmimo(kk) =  log2(real(det( eye(Nt) + SNR/Nt*(Hmimo)*(Hmimo)' )));
        CmimoMC(kk) = log2(real(det( eye(Nt) + SNR/Nt*Hmimo*txMCCorrMtx*Hmimo' )));  
       
        
    end
    
    % CDF estimates
    xi = linspace(0,xi_range(SNRidx),10001);
    CDFmimo = ksdensity(Cmimo,xi,'function','cdf');
   CDFMimoMC = ksdensity(CmimoMC,xi,'function','cdf');
    
    % Capacoty with outage 10%
    [no, idx] = max (sign(CDFmimo - 0.1));
    Cout_mimo(SNRidx) = xi(idx-1);
    [nq, idy] = max (sign(CDFMimoMC - 0.1));
    Cout_mimoMC(SNRidx) = xi(idy-1);
   
    
end
plot(SNRdBvalues,Cout_mimo,'linewidth',2);
hold on
plot(SNRdBvalues,Cout_mimoMC,'linewidth',2,'linestyle','--');
hold on
end
%plot
% figure(1)
% plot(SNRdBvalues,Cout_mimo,'r');
% hold on

xlabel('Average SNR \rho (dB)')
ylabel('Capacity with outage (bits/transmission)')
title('Capacity with 10% outage for i.i.d. Rayleigh slow fading channel')
legend('2X2','2X2 with MC', '4x4','4x4 with MC', '8X8','8X8 with MC', '16X16','16X16 with MC')
grid
hold off
