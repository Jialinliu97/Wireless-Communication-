%% 1×2  MRC Scheme With Rayleigh fading
clc
close all
clear all
%SISO vs. Receive Diversity
%% Theoretical BER result for Rayleigh with 8-psk :equation(50)
SNR = 0:1:20;                %SNR/3=Eb/No from 5 to 20
snr=(1/3)*10.^(SNR/10);		%Trans it to linear scale (SNR) for 8PSK
a=3*snr*((cos(pi/8))^2);
b=tan(pi/8);
c=3*snr*((cos(3*pi/8))^2);
d=tan(3*pi/8);
e=sqrt(1-1./(a+1));
f=sqrt(1-1./(c+1));
a1=1/3*e.*(0.5+1/pi*atan(e*b));
b1=1/3*f.*(0.5+1/pi*atan(f*d));
Ber_8psk_Ray=0.5-a1-b1;
%% Theoretical BER result for AWGN with 8-psk :equation(51)
snr2=(1/3)*10.^(SNR/10);		%linear format
    a2=(1/3)*(erfc(sqrt(3* snr2)*cos(pi/8))); %Equation(51)
    b2=(1/3)*(erfc(sqrt(3* snr2)*cos(pi*3/8)));
    Ber_8psk_Awgn=a2+b2.*(1-1.5*a2);        % BER for 8psk_AWGN 
%% Use Built-in function for  square root raised cosine filter
%b = rcosdesign(beta,span,sps,shape);12000 symbols
%The number of sampling points per symbol is 3 %the factor of roll-off is 0.3
rolloff = 0.3;       % Rolloff factor
span = 13;           % Filter span in symbols-13    from -6Ts to 6Ts 
sps = 10;             % Samples per symbol-10
h = rcosdesign(0.3,12,10,'sqrt');%the total samples 12*10=120
%% Square root raised cosine filter/Matched filter (not built-in function)
Ts=1; %Symbol interval
alpha=0.3; %The roll-off factor
t=-6*Ts:Ts/1000:6*Ts;%Filter delay is 12 symbols time interval width
%g=(sin((pi*t)./Ts)./(pi*t)).*cos((pi*alpha*t)./Ts)./(1-power(4*alpha*t/(2*Ts),2)); 
a1=cos((1+alpha)*pi*t/Ts);
b1=sin((1-alpha)*pi*t/Ts).*((4*alpha*t/Ts).^-1);
c1=pi*sqrt(Ts)*(1-16*(alpha^2)*(t.^2)/(Ts^2));
g=4*alpha*(a1+b1)./c1 ; %Root raised cosine Filter
%%  common simulation parameters
frmLen = 100;       % frame length
numPackets = 1200;  % number of packets /total bit is 120000
snr3 = 0:1:20;      % SNR from 0 to 20 dB
N = 2;              % maximum number of Tx antennas
M = 2;              % maximum number of Rx antennas
%% Set up simulation
% Create comm.8PSKModulator and comm.8PSKDemodulator System objects(TM)
P = 8;				% modulation order 8PSK
psk8Mod=comm.PSKModulator(8,'PhaseOffset',0);   %MOdulation without PhaseOffset
psk8Demod=comm.PSKDemodulator(8,'PhaseOffset',0);%Demodulation without PhaseOffset

% Create comm.OSTBCEncoder and comm.OSTBCCombiner System objects
ostbcEnc = comm.OSTBCEncoder;
ostbcComb = comm.OSTBCCombiner;
awgn1Rx = comm.AWGNChannel(...'NoiseMethod', 'Signal to noise ratio (Eb/No)', ...
'SignalPower', 1/3);
awgn2Rx = clone(awgn1Rx);

% Create comm.ErrorRate calculator System objects to evaluate BER.
errorCalc1 = comm.ErrorRate;
errorCalc2 = comm.ErrorRate;
errorCalc3 = comm.ErrorRate;

% Pre-allocate variables for speed
H = zeros(frmLen, N, M);
ber_noDiver  = zeros(3,length(snr3));
ber_Alamouti = zeros(3,length(snr3));
ber_MaxRatio = zeros(3,length(snr3));

% Set up a figure for visualizing BER results
fig = figure;
grid on;
ax = fig.CurrentAxes;
hold(ax,'on');
ax.YScale = 'log';
xlim(ax,[snr3(1), snr3(end)]);
ylim(ax,[1e-4 1]);
xlabel(ax,'SNR (dB)');
ylabel(ax,'BER');
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
fig.Name = 'MRC receive diversity Simulation';
title(ax,'8PSK-Receive diversity Simulation');
set(fig, 'DefaultLegendAutoUpdate', 'off');
fig.Position = figposition([15 50 25 30]);

% Loop over several SNR points
for idx = 1:length(snr3)
    reset(errorCalc1);
    reset(errorCalc2);
    reset(errorCalc3);
     awgn1Rx.EbNo= snr3(idx);% Set the EbNo property of the AWGNChannel System objects
     awgn2Rx.EbNo = snr3(idx);
    % Loop over the number of packets
    for packetIdx = 1:numPackets
        % Generate data vector per frame
        data = randi([0 P-1], frmLen, 1);% Creat symbols
        % Modulate data(8PSK)
        modData = psk8Mod(data);
        % Alamouti Space-Time Block Encoder
        encData = ostbcEnc(modData);
        % Create the Rayleigh distributed channel response matrix
        %   for  transmit and receive antennas
        H(1:N:end, :, :) = (randn(frmLen/2, N, M) + ...
                           1i*randn(frmLen/2, N, M))/sqrt(2);
        H(2:N:end, :, :) = H(1:N:end, :, :);

        % Extract part of H to represent the 1x1, 2x1 and 1x2 channels
        H11 = H(:,1,1);
        H21 = H(:,:,1)/sqrt(2);
        H12 = squeeze(H(:,1,:));

        % Pass through Rayleigh fading channels
        chanOut11 = H11 .* modData;
        chanOut21 = sum(H21.* encData, 2);
        chanOut12 = H12 .* repmat(modData, 1, 2);

        % Add AWGN
        rxSig11 = awgn1Rx(chanOut11);
        rxSig21 = awgn1Rx(chanOut21);
        rxSig12 = awgn2Rx(chanOut12);

        % Alamouti Space-Time Block Combiner
        decData = ostbcComb(rxSig21, H21);

        % ML Detector (minimum Euclidean distance)
        demod11 = psk8Demod(rxSig11.*conj(H11));
        demod21 = psk8Demod(decData);
        demod12 = psk8Demod(sum(rxSig12.*conj(H12), 2));%line sum

        % Calculate and update BER for current EbNo value
        %   for uncoded 1x1 system
        ber_noDiver(:,idx)  = errorCalc1(data, demod11);
        %   for Alamouti coded 2x1 system
        ber_Alamouti(:,idx) = errorCalc2(data, demod21);
        %   for Maximal-ratio combined 1x2 system
        ber_MaxRatio(:,idx) = errorCalc3(data, demod12);

    end % end of FOR loop for numPackets

    % Plot results
    %Simulation result
    semilogy(ax,snr3(1:idx), ber_noDiver(1,1:idx), 'r*', ...
         snr3(1:idx), ber_MaxRatio(1,1:idx), 'bs');
    hold on
    %Theoretical Result
    semilogy(SNR,Ber_8psk_Awgn,'k-o');
    hold on
    semilogy (SNR,Ber_8psk_Ray,'r-o');
          legend(ax,'No Diversity Simulation(SISO)', 'Receive Diversity (1Tx, 2Rx)', ...
           'AWGN Theoretical(SISO)','Rayleigh Theoretical(SISO)');

    drawnow;
end  % end of for loop for EbNo

% Perform curve fitting and replot the results
fitBER11 = berfit(snr3, ber_noDiver(1,:));
fitBER21 = berfit(snr3, ber_Alamouti(1,:));
fitBER12 = berfit(snr3, ber_MaxRatio(1,:));
semilogy(ax,snr3, fitBER11, 'r', snr3, fitBER12, 'b');
hold(ax,'off');

% Restore default stream
% rng(s);