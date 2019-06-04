clear 

% number of channel realization
It = 10000;

% SNR range in dB
SNRdBvalues = -10:5:15;

% range of values for each SNR to estimate CDF
xi_range = [1 2 4 8 12 18 25 35 40];

% initialize variables
Cmimo4x4 = zeros(1,It);
Csimo1x4 = zeros(1,It);
Cmiso4x1 = zeros(1,It);
Csiso = zeros(1,It);


SNRidx = 0;
for SNRdB = SNRdBvalues
    
    SNRdB
    SNR = 10.^(SNRdB./10);  % linear scale
    SNRidx = SNRidx + 1;
    
    % collect realizations of the maximal achievable rate
    for kk=1:It
        
        Hmimo4x4 = ( randn(4) + 1i*randn(4) )/sqrt(2);
        Cmimo4x4(kk) =  log2(real(det( eye(4) + SNR/4*Hmimo4x4*Hmimo4x4' )));
        
        Hsimo1x4 = ( randn(4,1) + j*randn(4,1) )/sqrt(2);
        Csimo1x4(kk) =  log2(real(det( 1 + SNR*norm(Hsimo1x4)^2 )));
        
        Hmiso4x1 = ( randn(1,4) + j*randn(1,4) )/sqrt(2);
        Cmiso4x1(kk) =  log2(real(det( 1 + SNR/4*norm(Hmiso4x1)^2 )));
        
        Hsiso = ( randn + j*randn )/sqrt(2);
        Csiso(kk) =  log2(real(det( 1 + SNR*norm(Hsiso)^2 )));
        
    end
    
    % CDF estimates
    xi = linspace(0,xi_range(SNRidx),10001);
    CDFmimo4x4 = ksdensity(Cmimo4x4,xi,'function','cdf');
    CDFsimo1x4 = ksdensity(Csimo1x4,xi,'function','cdf');
    CDFmiso4x1 = ksdensity(Cmiso4x1,xi,'function','cdf');
    CDFsiso    = ksdensity(Csiso,xi,'function','cdf');
    
    % Capacoty with outage 10%
    [no, idx] = max (sign(CDFmimo4x4 - 0.1));
    Cout_mimo4x4(SNRidx) = xi(idx-1);
    
    [no, idx] = max (sign(CDFsimo1x4 - 0.1));
    Cout_simo1x4(SNRidx) = xi(idx-1);
    
    [no, idx] = max (sign(CDFmiso4x1 - 0.1));
    Cout_miso4x1(SNRidx) = xi(idx-1);
    
    [no, idx] = max (sign(CDFsiso - 0.1));
    Cout_siso1x1(SNRidx) = xi(idx-1);
    
end


%plot
figure(1)
plot(SNRdBvalues,Cout_mimo4x4,'r');
hold on
plot(SNRdBvalues,Cout_simo1x4,'b');
plot(SNRdBvalues,Cout_miso4x1,'k');
plot(SNRdBvalues,Cout_siso1x1,'m');
xlabel('Average SNR \rho (dB)')
ylabel('Capacity with outage (bits/transmission)')
title('Capacity with 10% outage for i.i.d. Rayleigh slow fading channel')
legend('4x4 MIMO', '1x4 SIMO', '4x1 MISO', '1x1 SIMO')
grid
hold off
