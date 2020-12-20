%% Project 3

%Demodulate the re-sampled symbols and calculate BER with signal-to-noise-ratio (SNR) being
% 5 to 20 dB. Note: you may need more than 10,000 symbols for each SNR run.
% Ts=1;%symbol interval
% Tsp=0.1Ts;%the sampling period

clc;
clear all;
close all;
%% Square root raised cosine filter ( not built-in function)
Ts=1; %Symbol interval
alpha=0.3; %The roll-off factor
t=-6*Ts:Ts/1000:6*Ts;%Filter delay is 12 symbols time interval width
%g=(sin((pi*t)./Ts)./(pi*t)).*cos((pi*alpha*t)./Ts)./(1-power(4*alpha*t/(2*Ts),2)); 
a1=cos((1+alpha)*pi*t/Ts);
b1=sin((1-alpha)*pi*t/Ts).*((4*alpha*t/Ts).^-1);
c1=pi*sqrt(Ts)*(1-16*(alpha^2)*(t.^2)/(Ts^2));
g=4*alpha*(a1+b1)./c1 ; %Squared Root raised cosine Filter
figure(1);
plot(t,g);
xlabel('t'),ylabel('g(t)'); 
title('Impulse Response of Raised cosine roll filter'); 
axis([-6 6 -0.4 1.2]);
legend("\beta=0.3");

load Hw20.dat  % the matrix for the SNR=20dB ,Fd=80Hz
%% Use Built-in function for  square root raised cosine filter
%b = rcosdesign(beta,span,sps,shape);12000 symbols
%The number of sampling points per symbol is 3 %the factor of roll-off is 0.3
rolloff = 0.3;       % Rolloff factor
span = 13;           % Filter span in symbols-13    from -6Ts to 6Ts 
sps = 10;             % Samples per symbol-10
h = rcosdesign(0.3,12,10,'sqrt');%the total samples 12*10=120
figure(2);
plot(h);
%axis
grid on;
xlabel('(n)Samples')
ylabel('Amplitude');
title(' Raised cosine rolloff filter');
legend("\beta=0.3");
%% Theoretical BER result for Rayleigh with 8-psk :equation(50)
SNR = 0:1:20;           %SNR/3=Eb/No from 5 to 20
snr=(1/3)*10.^(SNR/10);		%Trans it to linear format(SNR) for 8PSK
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
%% Simulation Parameters 
N=360000;                    %the number of binary numbers
M=8;                        % 8PSK
Fd = 1;    %Symbol sampling frequency
Fs = 10;   %Raised cosine shaping filter sampling frequency
a = 0.3;   % different roll-off factor
delay = 6; %Filter delay
x=randi([0,1],1,N);         %Creat binary numbers
x1=fix(2*rand(1,30));       %30 random binary data
dataInMatrix=reshape(x,length(x)/3,3);
dataSymbolsIn=bi2de(dataInMatrix); %Creat symbols
R=raylrnd(0.5,N/3,1);         %Creat Rayleigh signal
h=pskmod(dataSymbolsIn,M); % 8PSK modulation
hR=h.*R;                   % add Rayleigh channel
%% Impulse Response of the filter
[num,den] = rcosine(Fd,Fs,'default',a,delay); %Design filter
srcFilter = rcosine(Fd,Fs,'sqrt',a,delay);
figure;
subplot(211)
impz(num,den);
xlabel('(n)Samples')
ylabel('Amplitude');
%The impulse response of the filter, num is the numerator coefficient of the transfer function, 
%  and 1 is the denominator coefficient of the transfer function
title('Impulse response of rms raised cosine filter with \beta=0.3');
[y,ty]=rcosflt(x,Fd,Fs,'filter',num,den); %Pulse shaping of binary data
t=delay:length(x)+delay-1;
subplot(212);
stem(t,x,'-b');
hold on %Draw binary data
plot(ty,y,'r') %draw data after Pulse shaping
legend('Original Binary data','\beta=0.3 Shaped message waveform')
axis([-1 40 -0.5 2])

%Carry out raised cosine shaping of the mapped 8PSK baseband signal
%% AWGN And Rayleigh
for i=1:length(SNR)
    %AWGN
    s1 = rcosflt(h,Fd,Fs,'filter',srcFilter);   %Pulse shaping after  modulation
    yAn=awgn(s1,SNR(i),'measured') ;
    rxAddNosieSrcFilter1 = rcosflt(yAn,Fd,Fs,'Fs/filter',srcFilter);%Matched filter
    rxDownsample1= downsample(rxAddNosieSrcFilter1,Fs);
    %Down-sampling, because the sender has been up-sampled
    % and the receiver performs down-sampling and recovery
    rxEnd1 = rxDownsample1(2*delay+1:end-2*delay); % Remove filter delay
    yA=pskdemod(rxEnd1,M);   %8PSK demodulation
    [bit_A,l]=biterr(dataSymbolsIn,yA);
     PSK_s_AWGN(i)=bit_A/N;
    
    %Rayleigh
    s2 = rcosflt(hR,Fd,Fs,'filter',srcFilter);%Pulse shaping after  modulation
    yRn=awgn(s2,SNR(i),'measured'); % Add AWGN          
    rxAddNosieSrcFilter2 = rcosflt(yRn,Fd,Fs,'Fs/filter',srcFilter);%Matched filter
    rxDownsample2 = downsample(rxAddNosieSrcFilter2,Fs)  ;% Remove filter delay
    rxEnd2 = rxDownsample2(2*delay+1:end-2*delay);
    yR=pskdemod(rxEnd2,M) ;   %8PSK demodulation
    [bit_R,ll]=biterr(dataSymbolsIn,yR);
    PSK_s_Ray(i)=bit_R/N; 
end


%% Compared BER with AWGN and Rayleigh fading
figure(4)
semilogy (SNR,Ber_8psk_Ray,'r-o');hold on
semilogy(SNR,PSK_s_Ray,'r-*');hold on
semilogy(SNR,Ber_8psk_Awgn,'b-o');hold on;
semilogy(SNR,PSK_s_AWGN,'b-*');hold on;
% axis([0 20 10^-5 1]);
axis([0 20 1e-3 1e0]);
xlabel('SNR[dB]');% Bit signal to noise ratio: Eb/N0
ylabel('BER');
%axis([-1 20 10^-4 1]);
legend('Rayleigh Theoretical','Rayleigh Simulation','AWGN Theoretical','AWGN Simulation');
title('8PSK BER with AWGN/Rayleigh');


