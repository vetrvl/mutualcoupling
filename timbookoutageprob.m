clear 

% number of channel realization. This number needs to be increased when
% estimating a small outage probability
It = 50000;

fprintf('Total number of channel realizations: %d\n', It);

% SNR range in dB
SNRdBvalues = -15:1:30;

% Selected transmission rate: 10bits/s
R=20;
Mtloop=[2 4 8 16];
for jj=1:length(Mtloop)


% initialize variables = number of channel realizations for which the
% selected rate exceeds the maximal achievable rate
N_mimo = zeros(1,length(SNRdBvalues));
N_mimocoup= zeros(1,length(SNRdBvalues));
Mt=Mtloop(jj);
%Mt=8;%Mrloop(jj);
fc=2.4e9;
numAnt=Mt;
txCorrMtx = eye(Mt);
txcoupmat=CouplingMatrix(0.5,fc, numAnt);
txMCCorrMtx = txcoupmat * txCorrMtx * txcoupmat';


for kk=1:It
    
    if mod(kk,10000) == 0
        fprintf('Number of channel realizations: %d\n', kk);
    end
    
    % generate channel realization
    Hmimo4x4 = ( randn(Mt) + 1i*randn(Mt) )/sqrt(2);
      
    SNRidx = 0;
    for SNRdB = SNRdBvalues
        
        SNR = 10.^(SNRdB./10);  % linear scale
        SNRidx = SNRidx + 1;
                HHcap=(Hmimo4x4)*(Hmimo4x4)';
                HHcapcoup=(Hmimo4x4)*(txMCCorrMtx)*(Hmimo4x4)';
        % number of channel realizations for which the selected rate exceeds the maximal achievable rate
        % MIMO 4x4
        Cmimo4x4 =  log2(real(det( eye(Mt) + SNR/Mt*HHcap )));
        N_mimo(SNRidx) = N_mimo(SNRidx) + (1-sign(Cmimo4x4 - R))/2;
        Cmimocoup =  log2(real(det( eye(Mt) + SNR/Mt*HHcapcoup )));
        N_mimocoup(SNRidx) = N_mimocoup(SNRidx) + (1-sign(Cmimocoup - R))/2;
                
    end
    
end

Pout_mimo4x4 = N_mimo/It;
Pout_mimocoup = N_mimocoup/It;
%plot
%figure
semilogy(SNRdBvalues,Pout_mimo4x4,'linewidth',2);
hold on
semilogy(SNRdBvalues,Pout_mimocoup,'linestyle','-.','linewidth',2);

end
xlabel('Average SNR \rho (dB)')
ylabel('Outage Probability')
title('Outage Probability vs SNR for slow Rayleigh fading channel')
grid

hold off
