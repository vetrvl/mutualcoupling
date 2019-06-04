clear 

% Maximum number of transmit and receive antennas
Mmax = 10;

% number of channel realization
It = 10000;


% initialize variables: ergodic capacity for each number of antennas
Cmimo  = zeros(1, Mmax);
Csimo  = zeros(1, Mmax);
Cmiso  = zeros(1, Mmax);
Csiso  = zeros(1, Mmax);

SNRdB  = 0;  % in dB
SNR = 10.^(SNRdB./10);  % linear scale

for kk=1:It
    
    for M = 1:Mmax
        
        % MIMO
        Hmimo = ( randn(M)   + j*randn(M)   )/sqrt(2);
        Cmimo(M) = Cmimo(M) + log2(real(det( eye(M) + SNR/M*Hmimo*Hmimo' )));
        
        %% SIMO
        hsimo = ( randn(M,1) + j*randn(M,1) )/sqrt(2);
        Csimo(M) = Csimo(M) + log2( 1 + SNR*norm(hsimo)^2);
        
        %% MISO
        hmiso = ( randn(M,1)  + j*randn(M,1)  )/sqrt(2);  
        Cmiso(M) = Cmiso(M) + log2( 1 + SNR/M*norm(hmiso)^2);
        
        %% SISO
        hsiso = ( randn + j*randn )/sqrt(2);
        Csiso(M) = Csiso(M) + log2( 1 + SNR*abs(hsiso)^2);
        
    end
    
end

% Compute average over all channel realizations
Csiso = Csiso/It
Csimo = Csimo/It
Cmiso = Cmiso/It
Cmimo = Cmimo/It

% plot
figure(1)
plot(1:Mmax, Cmimo,'r')
hold on
plot(1:Mmax, Csimo,'b')
plot(1:Mmax, Cmiso,'k')
plot(1:Mmax, Csiso,'m')
xlabel('Number of antennas M')
ylabel('Ergodic Capacity (bits/transmission)')
title('Ergodic Capacity for i.i.d. Rayleigh fast fading channel - SNR=0dB')
legend('4x4 MIMO', '1x4 SIMO', '4x1 MISO', '1x1 SISO')
grid
hold off