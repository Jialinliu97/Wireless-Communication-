%% Rayleigh fading channel
%   R = Rayleigh(Ns, Fd, Ts)
%   - simulates Rayleigh fading with a Doppler frequency
%     of Fd Hz, returning an N-(row)vector of samples taken 
%     Ts seconds apart
%   - returns the complex channel response, R, with unit
%     variance

function R = Rayleigh(Ns, Fd, Ts)
% function R =Rayleigh(5000, 1, 0.01)
  clear all
  rand('seed',0);
  Fd=1  ;
  Ts=1/100;
  Ns=360/Ts; % the number of samples is 36000

%% fading channel parameters 

      M=8;
      N=4*M;
      dop_gain = sqrt(2/M) * ones(1,M);
      theta=(2*rand-1)*pi;
      doppler = Fd * cos(2*pi/N*(1:M)+theta/N-pi/N);
      phi=(2*rand(M,1)-1)*pi;
      varphi=(2*rand(M,1)-1)*pi;
      state = zeros(M,1);
      dop_update = (2*pi*doppler * Ts).';	
      
%% generate fading channel samples 

    R = zeros(1,Ns);
    for (k=1:Ns)
        R(k) = dop_gain * [cos(state+phi)+1i*sin(state+varphi)];
        state = dop_update + state;
    end  
return


     


