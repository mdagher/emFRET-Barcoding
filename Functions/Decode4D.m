function [ Unique_Red, I_norm_model_R, Dat_Final_Result, Clusters, posterior_red, Posterior_Blue_Final, num_clust_blue_output ] = Decode4D( BC, BeadbyBead, BCnum_col, Var_Col_nb, Beta, BB_MFI, Model )
%This function performs 4D classification and, if model is used, barcode
%identification simultaneously 
%   Detailed explanation goes here



%% FIND THE UNIQUE RED CLASSIFIERS

%TO DO MICHAEL: This works for 4p4 (2 red classifiers) but not for many red
%classifiers - 
Unique_Red (1,:) = BC(1,3:4); %changed from 3 to 1
for i=2:size(BC,1)
    count = 0;
    for j=1:size(Unique_Red,1)
        if BC(i,3:4) ~= Unique_Red(j,:)
            count = count + 1;
        end
    end
    if count == size(Unique_Red,1)
        Unique_Red (size(Unique_Red,1)+1,:) = BC(i,3:4);
    end
end
%Unique_Red = BC(:,3:4);

%% Crunch the model for red classifiers to get expected values, and laser normalized values
concentration = zeros(size(Unique_Red,1),4);
concentration(:,3:4) = Unique_Red;
[I_model_SBR_R, E_model_R]= MFM_SBR(concentration); %the expected value is stored in I_model_SBR_4p4_LR


m=size(Unique_Red,1);
n=4;
for i=1:m
    for j=1:n
        I_model_R(i,j) = I_model_SBR_R(i,j) .* BB_MFI(j) ;
        I_norm_model_R(i,j)= Beta(j).*(I_model_R(i,j)-BB_MFI(j)) + BB_MFI(j) ;
    end 
end 

%% Cluster groups for the red and append the cluster index to Dat_Final_Result
% Cluster red with model even when not using model when testing for
% clustering errors


[clusterX_red, posterior_red, gmfit_red] = BC_ClusterFunc (BeadbyBead(:,Var_Col_nb(7:8)), size(Unique_Red,1), I_norm_model_R, 3, 4, 1);  %it shouldnt be 2, should be number of red classifiers  
Dat_Final_Result = zeros(size(BeadbyBead,1),6);
Dat_Final_Result (:,1:4) = BeadbyBead(:,Var_Col_nb(5:8));
Dat_Final_Result (:,5) = clusterX_red;

%% for each of these groups/clusters for the reds (that we have), see the number of unique blues

%if Model==0
%Unique_Red=sortrows(Unique_Red,-1);
%end 

cumsum = 0;
Clusters=[];
num_clust_blue_output = [];
for z=1: size(Unique_Red,1)
    num_clust_blue = size(find(BC(:,3) == Unique_Red(z,1) & BC(:,4) == Unique_Red(z,2)),1);
    num_clust_blue_output = [num_clust_blue_output; num_clust_blue];
    blue_BC_ind = find(BC(:,3) == Unique_Red(z,1) & BC(:,4) == Unique_Red(z,2));
    
    %get the expected values for each blue (using MFM)
    [I_model_SBR_B, E_model_B]= MFM_SBR(BC); %not just blue '_B'
    m=num_clust_blue;
    %blue_BC_ind = find(BC(:,3) == Unique_Red(z,1) & BC(:,4) == Unique_Red(z,2));
    blue_BC = BC(blue_BC_ind,1:2);
    n=4; % number of dyes 
    
    %get normalized model params
    I_norm_model_B=[];
    for i=1:m
        for j=1:n
            I_model_B(i,j) = I_model_SBR_B(i+cumsum,j) .* BB_MFI(j) ;
            I_norm_model_B(i,j)= Beta(j).*(I_model_B(i,j)-BB_MFI(j)) + BB_MFI(j) ;
        end 
    end 
    
    %run the GMM with the ICs specified by the model
    Dat_Blue_Ind = find(clusterX_red == z);
    [Cluster_Blue , posterior_blue, gmfit_blue] = BC_ClusterFunc (BeadbyBead(Dat_Blue_Ind,Var_Col_nb(5:6)), num_clust_blue, I_norm_model_B, 1, 2, Model);
    Posterior_Blue_Final(Dat_Blue_Ind,1+cumsum:num_clust_blue+cumsum) = posterior_blue; %MK
    Dat_Final_Result (Dat_Blue_Ind,6) = Cluster_Blue;
    %add red barcode concentrations to table
    Dat_Final_Result (Dat_Blue_Ind,9) = Unique_Red(z,1);
    Dat_Final_Result (Dat_Blue_Ind,10) = Unique_Red(z,2);
    cumsum = cumsum + num_clust_blue;
    
    
    %add blue barcode concentrations to table
    for p=1:num_clust_blue
        clust_blue_ind = find(Dat_Final_Result(:,6) == p & Dat_Final_Result(:,5) == z);
        Dat_Final_Result (clust_blue_ind,7) = blue_BC(p,1);
        Dat_Final_Result (clust_blue_ind,8) = blue_BC(p,2);
        Cluster_add = [p,z,gmfit_blue.mu(p,1), gmfit_blue.mu(p,2),gmfit_red.mu(z,1), gmfit_red.mu(z,2)];
        Clusters = [Clusters ; Cluster_add];
    end
    
    %remove vars that are used in loop 
    clear I_model_4p4_LB;
    clear I_norm_model_4p4_LB
end
%%


Dat_Final_Result (:,11) = BeadbyBead(:,BCnum_col);
% put expected barcode number in table; %%TODO: FIX THIS
for i=1:size(Dat_Final_Result,1)
    for j=1:size(BC,1)
        if Dat_Final_Result(i,7) == BC(j,1) && Dat_Final_Result(i,8) == BC(j,2) && ...
                Dat_Final_Result(i,9) == BC(j,3) && Dat_Final_Result(i,10) == BC(j,4)
            Dat_Final_Result(i,12) = j;
        end 
    end
end


end

