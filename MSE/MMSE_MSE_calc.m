%Function Declaration: 
function ms_error=MMSE_MSE_calc(X,H,Y,Rgg,variance); 
%This function generates mean squared error for the the MMSE estimator.. 
%EVALUATION OF Hmmse 
%Hmmse=F*Rgg*inv(Rgy)*Y; 
u=rand(64,64); 
F=fft(u)*inv(u);%The 64 X 64 twiddle factor matrix.. 
I=eye(64,64); 
Rgy=Rgg * F'* X'; 
Ryy=X * F * Rgg * F' *X' + variance * I; 
for i=1:64 
    yy(i,i)=Y(i); 
end 
Gmmse=Rgy * inv(Ryy)* Y; 
Hmmse=fft(Gmmse); 
 
ms_error_mat=mean(((abs(H)-abs(Hmmse))/abs(H)).^2); 
for i=1:64 
    if(ms_error_mat(i)~=0) 
        ms_error=ms_error_mat(i); 
    end 
end
