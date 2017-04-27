function [ rejection_rate_vector, error_rate_vector, percentage_final ] = BC_PosteriorProbabilityFigures( Dat_Final_Result,num_blue_clusters, posterior_blue, posterior_red )
%BC_POSTERIORPROBABILITYFIGUES Summary of this function goes here
%   Detailed explanation goes here
%Modified April 9, 2017

%% Data Initialization
%Make Posterior probability table for the 4D decode, and store it in P
%The global posterior probability table is created by multiplying the posterior
%probability for the red and blue classifiers
index = 0;
for j=1:size(num_blue_clusters,1) 
    for k=1:num_blue_clusters(j)
        index = index + 1;
        for i=1:size(posterior_red,1)
            P(i,index) = posterior_blue(i,index) * posterior_red(i,j);
        end
    end
end

assigned_real_BC = Dat_Final_Result(:,11:12);
assigned_real_BC = [assigned_real_BC, P]; %concatenate with the posterior probability table, so that the table is sorted in the same way
sorted_clusters = sortrows(assigned_real_BC,2);
P = sorted_clusters(:,3:2+index);
clusters = sorted_clusters(:,1);
orig_clusters = sorted_clusters(:,2);
for i=1:index
    numberPerPop(i) = length(find(Dat_Final_Result(:,12) == i));
end
n = sum(num_blue_clusters); %number of total clusters found (40)
num_beads = size(clusters,1);

%% Varying Posterior Probability Threshold
%iterate over posterior probabilities from 0.5-0.99, and reject beads that
%lie below the posterior probability threshold (the likelihood of belonging
%to a particular cluster
for z=50:99
%z=50;
    %filter by posterior probability
    clusterX_Good = zeros(num_beads,1);
    accepted_count =0;
    for i= 1:n
       for j=1:(num_beads)
           if P(j,i) > z/100 %confidence value
               clusterX_Good(j) = clusters(j);
               accepted_count = accepted_count + 1;
           end
       end
    end
    total_rejected = num_beads - accepted_count;


    %get the count and percentage matrix (where each cluster was assigned)
    percentage = zeros(n,n); %heat map
    count = zeros(n,n);
    index = 0;
    for i=1:n
        for j=1:numberPerPop(i)
            index = index + 1;
            for k=1:n
                if clusterX_Good(index) == k
                    count(i,k) = count(i,k) + 1;
                end
            end
        end
    end

    %calculate numberPerPop remaining
    numberPerPop_remaining = zeros(n,1);
    index = 0;
    for i=1:n
        count_zero =0;
        for j=1:numberPerPop(i)
            index = index + 1;
            if clusterX_Good(index) == 0;
                count_zero = count_zero + 1;
            end
        end
        numberPerPop_remaining(i) = numberPerPop(i) - count_zero;
    end

    %calculate percentages for bead identification
    for i=1:n
        for j=1:n
            percentage(i,j) = 100*(count(i,j)/numberPerPop_remaining(i));
        end
    end

   
    %This part is quite confusing. Basically, each cluster is identified
    %from the initial cluster (from the model) for which there is the maximum overlap
    %(it could also be measured by percentage overlap, as opposed to
    %absolute overlap, or 'count'). Once one cluster was used, it cannot be
    %reassigned, the rows, and columns are removed
    col_inds = zeros(n,1);
    row_inds = zeros(n,1);
    rows = transpose(1:n);
    cols = transpose(1:n);

    copy_percentage = count;
    for k=1:n
        [Y,I] = max(copy_percentage);
        [J,M] = max(Y);
        row_inds(k) = rows(I(M));
        col_inds(k) = cols(M);
        copy_percentage(I(M),:) = []; 
        copy_percentage(:,M) = [];
        rows(I(M)) = []; %delete row
        cols(M) = []; %delete column
    end

    %identify cluster according to algorithm described above
    %note that clusters = 0 corresponds to rejection
    count_error = 0;
    count_rejection = 0;
    properCluster_inds = zeros(n,1);
    for i=1:n
        for j=1:n
            if row_inds(j) == i
                properCluster_inds(i) = col_inds(j);
            end
        end
    end

    %assign each bead to its "proper cluster"
    properCluster = zeros(num_beads,1);
    properCluster(1:numberPerPop(1)) = properCluster_inds(1);
    for i=2 : n
        properCluster(sum(numberPerPop(1:i-1)) + 1 : sum(numberPerPop(1:i))) = properCluster_inds(i);
    end

    %See what differred between beads (after thresholding with posterior
    %prob and propercluster)
    for i= 1:num_beads
        if clusterX_Good(i) ~= properCluster(i) && clusterX_Good(i) ~= 0
            count_error = count_error + 1;
        elseif clusterX_Good(i) == 0
            count_rejection = count_rejection + 1;
        end
    end
    
    %get rates
    total_rejected = num_beads - accepted_count;
    rejection_rate = 100*(total_rejected/(num_beads));
    error_rate = 100*(count_error/((num_beads) - total_rejected));
    error_rate_vector(z) = error_rate;
    rejection_rate_vector(z) = rejection_rate;
    percentage_final(:,:,z-49) = percentage;

end %z
%total_error = sum(errors);
%error_rate = 100*(total_error/num_beads);

%% Plots
xAxis = linspace(0.5,0.99,50);
figure;
plot(xAxis,error_rate_vector(50:99));
xlabel('Posterior Probability, P'); ylabel('Error Rate (%)');
figure;
plot(xAxis,rejection_rate_vector(50:99))
xlabel('Posterior Probability, P'); ylabel('Rejection Rate (%)');


end

