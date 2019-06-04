
%clear 

%Mtloop = [2 4 8];  % number of transmit antennas
Mtloop=4;
Mrloop = Mtloop;  % number of receive antennas
fc=2.4e9;
for jj=1:length(Mtloop)
Mt=Mtloop(jj);
Mr=Mrloop(jj);
numAnt=Mt;
txCorrMtx = eye(Mt);
txcoupmat=CouplingMatrix(0.5,fc, numAnt);
txMCCorrMtx = txcoupmat * txCorrMtx * txcoupmat';

% number of channel realization
It = 10000;

% SNR range in dB
SNRdBvalues = -15:30;

% initialize variables: ergodic capacity for each value of SNRdBvalues
Cmimo  = zeros(1, length(SNRdBvalues));
CmimoMC  = zeros(1, length(SNRdBvalues));


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
%figure(1)
plot(SNRdBvalues, Cmimo)
hold on
plot(SNRdBvalues, CmimoMC)
hold on
end
xlabel('Average SNR \rho (dB)')
ylabel('Ergodic Capacity (bits/transmission)')
title('Ergodic Capacity for Rayleigh fast fading channel')
%legend('2x2 MIMO without coupling', '2x2 with coupling')
grid
hold on
 