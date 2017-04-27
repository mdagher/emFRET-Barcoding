function [ D , E ] = MFM_SBR( f )

%CRUNCHMODEL takes in the input stoichiometry of fluorphores (LOs) (f) and computes the expected
%detector intensities in SBR (D) and also returns the emFRET efficiencies
%(E)

% Output D = [D1, D2, D3, D4]
% Output E = [E1e,E12,E13,E14,E2e,E23,E24,E34]

% This is strictly for the four dyes (FAM,Cy3,Cy5,Cy5.5) on FACS CANTO II 

[m,n]=size(f);

% The effective forster radii for donors with multiple acceptors 
R1e=zeros(m,1);
R2e=zeros(m,1);

% The inter-dye FRET efficiencies that must be computed. 
E1e=zeros(m,1);
E12=zeros(m,1);
E13=zeros(m,1);
E14=zeros(m,1);
E2e=zeros(m,1);
E23=zeros(m,1);
E24=zeros(m,1);
E34=zeros(m,1);

% The detector SBR  
D=zeros(m,n);

%% Parameters of model extracted from fitting using Prism and empirical fits to the equations(Check BC4.0)

% The mu matrix -ie. direct excitation
mu=[0.8845 0.11 4.296 1.961]; %

% The beta matrix-ie. bleedthrough-
beta=   [1 0.1607 0 0;
      0.0728 1 0 0 ; 
      0 0 1  0.3527; 
      0 0 0.0164 1];

% The Forster radius matrix - calculated from spectra 

R= [NaN 5.5e-9 4.5e-9 4.3e-9; 
    NaN    NaN 5.3e-9 4.9e-9;
    NaN    NaN    NaN 6.7e-9;
    NaN    NaN    NaN   NaN];

% The FRET proportionality matrix
alpha= [NaN  0.27   NaN  NaN; 
       NaN    NaN   NaN   NaN;
       NaN    NaN   NaN   0.6; 
       NaN    NaN   NaN   NaN]; 
   

% Constants
%gamma=0.61; 
%lambda=1.1;
%t=1e7; 


%% Calculation of parameters and estimated output response
for i = 1:m

% Calculate effective Forster radii and emFRET efficiencies   
R1e(i)=sqrt( (f(i,2)./(f(i,2)+f(i,3)+f(i,4))).*R(1,2)^2 + (f(i,3)/(f(i,2)+f(i,3)+f(i,4))).*R(1,3)^2 + (f(i,4)/(f(i,2)+f(i,3)+f(i,4))).*R(1,4)^2); 
E1e(i)=Efret2(R1e(i), f(i,2)+f(i,3)+f(i,4));
E12(i)=E1e(i).*(f(i,2).*R(1,2)^2)./((f(i,2)+f(i,3)+f(i,4)).*R1e(i)^2); 
E13(i)=E1e(i).*(f(i,3).*R(1,3)^2)./((f(i,2)+f(i,3)+f(i,4)).*R1e(i)^2); 
E14(i)=E1e(i).*(f(i,4).*R(1,4)^2)./((f(i,2)+f(i,3)+f(i,4)).*R1e(i)^2); 

R2e(i)=sqrt((f(i,3)./(f(i,3)+f(i,4))).*R(2,3)^2 + (f(i,4)/(f(i,3)+f(i,4))).*R(2,4)^2 );
E2e(i)=Efret2(R2e(i), f(i,3)+f(i,4)); 
E23(i)=E2e(i).*(f(i,3).*R(1,3)^2)./((f(i,3)+f(i,4)).*R1e(i)^2); 
E24(i)=E2e(i).*(f(i,4).*R(1,4)^2)./((f(i,3)+f(i,4)).*R1e(i)^2); 

E34(i)=Efret2(R(3,4), f(i,4));

if f(i,2)+f(i,3)+f(i,4)==0
    E1e(i)=0;
    E12(i)=0;
end

if f(i,3)+f(i,4)==0
    E2e(i)=0; 
end


% Calculate SBR for all four detectors
D(i,1)= 1 + mu(1)*f(i,1)*(1-E1e(i)) + beta(2,1) *( mu(2).*f(i,2) + alpha(1,2) .* E12(i).* mu(1) .* f(i,1)) .* (1- E2e(i));

D(i,2)= 1 + beta(1,2).*mu(1).*f(i,1).*(1-E1e(i)) + ( mu(2) * f(i,2) + alpha(1,2).* E12(i).*mu(1).*f(i,1)).*(1-E2e(i)); 

D(i,3)= 1 + mu(3).*f(i,3)*(1-E34(i)) + beta(4,3) .* ( mu(4).*f(i,4) + alpha(3,4).*E34(i).*mu(3).*f(i,3) );

D(i,4)= 1 + beta(3,4).*mu(3).*f(i,3)*(1-E34(i)) + mu(4).*f(i,4) + alpha(3,4).*E34(i).*mu(3).*f(i,3) ; 

end

% Concatenate the emFRET efficiencies
E=cat(2,E1e,E12,E13,E14,E2e,E23,E24,E34);

end

