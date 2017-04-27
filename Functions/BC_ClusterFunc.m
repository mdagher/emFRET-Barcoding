function [ ClusterX, P, GMfit ] = BC_ClusterFunc( X, numberClusters, Exp, IndBot, IndTop, Model)
%GAUSSIAN_CLUSTERFUNC Summary of this function goes here
%   X: data to cluster
%   numbmerClusters: number of clusters to assign datapoints to
%   Exp: expected means for each cluster (from the barcoding model)
%   IndBot/IndTop specify which intensities of the model to use for the
%   initial cluster mean
%   For clustering red data: IntBot=3, IndTop=4
%   For clustering blue data: IntBot=1, IndTop=2
%   Model: 1 for with model, 0 for without model

CV = 0.1;
options = statset('MaxIter',5000, 'TolFun',1e-7); % Increase number of EM iterations
if Model == 1
    mu = Exp(:,IndBot:IndTop);
else
    Indices = round(rand(numberClusters,1)*size(X,1));
    for i=1:numberClusters
        mu(i,:) = X(Indices(i),:);
    end
end


Comp_Prop = 1/numberClusters*ones(numberClusters,1); % define
Sigma = zeros(2,2,numberClusters);
for i=1:numberClusters
    SD1 = mu(i,1)*CV;
    SD2 = mu(i,2)*CV;
    var1 = SD1^2;
    var2 = SD2^2;
    Sigma(1,1,i) = var1;
    Sigma(2,2,i) = var2;   
end
        
S = struct('mu', mu, 'Sigma', Sigma, 'ComponentProportion', Comp_Prop); %initial conditions
GMfit = fitgmdist(X,numberClusters,'Options', options, 'Start', S);
P = posterior(GMfit,X);
ClusterX = cluster(GMfit,X); %NOTE TO MILAD: function to investigate for cluster numbering


%plot
figure;
h1 = gscatter(X(:,1),X(:,2),ClusterX); %this just plots X, not X_Good
%h1 = gscatter(X(:,1),X(:,2));
hold on;
loglog(GMfit.mu(:,1),GMfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',20)
ax=gca;
set(gca, 'xscale', 'log', 'yscale','log', 'XLim', [0, 3e5], 'YLim', [0, 1e5]);
ax.XTick=[ 1e2 1e3 1e4 1e5 ]; 
ax.YTick=[ 1e1 1e2 1e3 1e4]; 
scatter(mu(:,1),mu(:,2));
%axis([0, 4e4, 0 1.5e4])
hXLabel = xlabel('Fluorescence 1');
hYLabel = ylabel('Fluorescence 2');
hold off

end

