function [ BeadbyBead_Singlets, BeadbyBead_beads ] = SingletGate( BB_MFI, BeadbyBead, Var_Col_Nb  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

BbB_Fluorescence= [BeadbyBead(:,Var_Col_Nb(5)),BeadbyBead(:,Var_Col_Nb(6)),...
    BeadbyBead(:,Var_Col_Nb(7)),BeadbyBead(:,Var_Col_Nb(8))];

%% Step 1: Removes 'beads' (read: crap) with intensity lower than bare beads. 

% This removes all beads with intensity <  0.5*bare beads
ind_to_clean = [];
for i=1:size(BbB_Fluorescence,1)
            count = 0;
    for j=1:4
        if BbB_Fluorescence(i,j) < 0.5*BB_MFI(j) 
            count = count + 1;
        end
    end
        if count >= 1
            ind_to_clean = [ind_to_clean; i];
        end
    
end

% remove background indices
BeadbyBead (ind_to_clean, :) = [];
BeadbyBead_beads=BeadbyBead; % Report back on number of beads (what is not crap) 
BbB_Scatter = [ BeadbyBead(:,Var_Col_Nb(1)), BeadbyBead(:,Var_Col_Nb(2)), ...
    BeadbyBead(:,Var_Col_Nb(3)), BeadbyBead(:,Var_Col_Nb(4))];


%% Step 2: CLUSTER TO FIND SINGLETS

% Clustering to find singlets - Cluster based on FSC-A and SSC-A
X = [BbB_Scatter(:,1), BbB_Scatter(:,3)];
numberClusters = 2;
gmfit = fitgmdist(X,numberClusters);
P = posterior(gmfit,X);
clusterX = cluster(gmfit,X); 


% Plot that shows the two clusters were appropriately found 
figure;
h1 = gscatter(X(:,1),X(:,2),clusterX); 
hold on;
plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',20)
hold on

%% Step 3: Doublets have higher values, KILL DOUBLETS!

cluster_to_kill = 2;
if gmfit.mu(2,1) < gmfit.mu(1,1)
    cluster_to_kill = 1;
end

% Get rid of doublets 
ind_doublets = find(clusterX == cluster_to_kill);
BeadbyBead_Singlets = BeadbyBead;
BeadbyBead_Singlets (ind_doublets, :) = []; 

% Verification that you're getting a single cluster here. 
figure;
scatter (BeadbyBead_Singlets(:,1), BeadbyBead_Singlets(:,3))

end
