%Function Declaration: 
function ms_error=LS_MSE_calc(X,H,Y); 
%This function generates mean squared error for the the LS estimator.. 
%EVALUATION OF Hls 
Hls =(inv(X)) * Y; 
%The simplest of 'em all indeed.. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
ms_error_mat=mean(((abs(H-Hls))/abs(H)).^2); 
for i=1:64 
    if(ms_error_mat(i)~=0) 
        ms_error=ms_error_mat(i); 
    end 
end 
 
