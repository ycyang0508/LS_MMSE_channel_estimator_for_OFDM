%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Author:      Vinay Uday Prabhu 
%   E-mail:      vinay_u_prabhu@yahoo.co.uk 
%   Function:    Comparison of the performances of the LS and the MMSE channel estimators 
%                for a 64 sub carrier OFDM system based on the parameter of Mean square error 
%  Assumptions: The channel is assumed to be g(t)=delta(t-0.5 Ts)+delta(t-3.5 Ts) 
%               {Fractionally spaced} 
%For more information on the theory and formulae used , please do refer to the paper On 
%"Channel Estimation In OFDM systems" By Jan-Jaap van de Beek, Ove Edfors, Magnus Sandell 
% Sarah Kate wilson and Petr Ola Borjesson In proceedings Of VTC'95 Vol 2 pg.815-819 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clc; 
clear all; 
%Generation of a naive training sequence.. 
%Assuming BPSK modulation ...symbols:+1/-1 
X=zeros(64,64); 
d=rand(64,1); 
      for i=1:64 
       if(d(i)>=0.5) 
           d(i)=+1; 
       else 
           d(i)=-1; 
       end 
    end 
 for i=1:64 
     X(i,i)=d(i); 
 end 
%Calculation of G[The channel Matrix] 
 %The channnel is...  
  tau=[0.5 3.5];%The fractionally spaced taps.. 
%Generation of the G matrix... 
for k=1:64 
      s=0; 
      for m=1:2 
         s=s+(exp(-j*pi*(1/64)*(k+63*tau(m))) * (( sin(pi*tau(m)) / sin(pi*(1/64)*(tau(m)-k))))); 
         %Go through the above cited paper for the theory behind the formula 
      end 
g(k)=s/sqrt(64); 
end 
G=g';%Thus, the channel vector is evaluated.. 
H=fft(G);% In the freq domain.. 
u=rand(64,64); 
F=fft(u)*inv(u);% 'F' is the twiddle factor matrix.. 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Evaluation of the autocovariance matrix of G-Rgg 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
gg=zeros(64,64); 
for i=1:64 
    gg(i,i)=G(i); 
end 
gg_myu = sum(gg, 1)/64;                     
gg_mid = gg - gg_myu(ones(64,1),:);         
sum_gg_mid= sum(gg_mid, 1); 
Rgg = (gg_mid' * gg_mid- (sum_gg_mid'  * sum_gg_mid) / 64) / (64 - 1); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Running for a dozen trials to try and average out the results.. 
for m=1:12 
     
for n=1:5 
 
SNR_send=5*n; 
XFG=X*H; 
n1=ones(64,1); 
n1=n1*0.000000000000000001i;%Just to ensure that the function awgn adds 'complex gaussian noise'.. 
noise=awgn(n1,SNR_send); 
variance=var(noise); 
N=fft(noise); 
Y=XFG+N; 
%Evaluating the mean squared error for the LS estimator.. 
mean_squared_error_ls=LS_MSE_calc(X,H,Y); 
%Evaluating the mean squared error for the MMSE estimator.. 
mean_squared_error_mmse=MMSE_MSE_calc(X,H,Y,Rgg,variance); 
SNR(n)=SNR_send; 
mmse_mse(m,n)=mean_squared_error_mmse; 
ls_mse(m,n)=mean_squared_error_ls; 
end; 
 
end; 
ls_mse 
mmse_mse 
mmse_mse_ave=mean(mmse_mse); 
ls_mse_ave=mean(ls_mse); 
%Now just the display part..... 
semilogy(SNR,mmse_mse_ave,'k-'); 
grid on; 
xlabel('SNR in DB'); 
ylabel('mean squared error'); 
title('PLOT OF SNR V/S MSE FOR AN OFDM SYSTEM WITH MMSE/LS ESTIMATOR BASED RECEIVERS'); 
 
hold on; 
semilogy(SNR,ls_mse_ave,'b*'); 
semilogy(SNR,ls_mse_ave,'b-'); 
semilogy(SNR,mmse_mse_ave,'kv'); 
grid on; 
xlabel('SNR in DB'); 
ylabel('mean squared error'); 
title('PLOT OF SNR V/S MSE FOR AN OFDM SYSTEM WITH MMSE/LS ESTIMATOR BASED RECEIVERS'); 
