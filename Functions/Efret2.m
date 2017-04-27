function [ E ] = Efret2(r,fa)
%EFRET Summary of this function goes here
%   Detailed explanation goes here
gamma=0.61;
lambda=1.1; 
t= 2.13e14 ; % this is estimated from BC4.0 #18 (0,0,16,20)
sigma_a=t*fa; 
nad=pi*r^2*sigma_a; 

% FRET equation from Koppel et. al eq (39)
E= (nad/gamma)^lambda/(1+(nad/gamma)^lambda); 

% and conversely 
%sigma_a= (gamma/(pi*r^2)) * (E/(1-E))^(1/lambda)



end

