% Milad Dagher, McGill University, April 1st 2017. 

% This file validates the automated decoding (model vs no model) whereby we
% use digitally concatenated barcodes (i.e. we have knowledge of barcode
% for every bead), scramble them, and characterize the automatic decoding's
% ability in (1) finding clusters and (2) identifying clusters 

% (1) Automatically decoding (w/ or w/o model) beads using GMM beads
% (2) Calculating the % of correctly clustered beads (i.e. have they been clustered with the right cluster)
% wrt to the scanning posterior probability - and plotting that. 
% (3) Decoding % = clustering + identification. Plotting as a heat map.
% PS. Identification without model is done by sorting the data

% The data used for this validation is BC4.4 (2 red clusters)

% PS. For fca_readfcs to work properly, MATLAB must be pointed to the working folder

clear all 

%% USE MODEL? 
model = 0;


%% Input 1: POINT TO BARCODES USED (XLS)

%load barcodes we want to Matlab (Barcode Repertoire)
bc=xlsread('mfiles/BarcodingManuscriptFigures/DecValidation/DecodingValidation_BCrep',2,'C3:F28');



%% Input 2: POINT TO FCS OF BARE BEADS (BBs), AND NORMALIZATION BEADS (NBs)  

% Bare Beads
[bb_fcsdat,bb_fcshdr,bb_scaled]=fca_readfcs('FCS_4p3/20160415_1.fcs');

% Normalization Beads 
% To do, check if one normalization bead, properly chosen, leads to similar
% normilzation values. 


%% Input 3: POINT TO FCS OF BARCODES - HEY WANNA SCRAMBLE?

filepath = 'FCS_4p3/';
filename_start = '20160415_';
filename_end = '.fcs';
populations = [2:8 10:22 24:27 29:31]; 

% Scramble the beads? 
scramble=1; 


%% Step 1: FCS FILES COLUMN INFORMATION
% Finds the columns of the parameters of interest. 

% How many parameters exported overall? 
fcscolsize= size(bb_fcsdat,2); 

% What are thy names? 
colnames= extractfield(bb_fcshdr.par, 'name');

% Where are the ones we are interested in? Specific for FCS CANTO. Must be edited for another cytometer.  
var_col_nb= FCScol_info_Canto(colnames); 

% Les voila, each of these variables is the column number of the parameter.
% e.g. I1col is the column number of the detector 1 intensity. mkay? 
FSCcol_A=var_col_nb(1); FSCcol_H=var_col_nb(2); SSCcol_A=var_col_nb(3); SSCcol_H=var_col_nb(4);
I1col=var_col_nb(5); I2col=var_col_nb(6); I3col=var_col_nb(7); I4col=var_col_nb(8);




%% Step 2: CONCATENATE THEM FCS FILES (validation specific) 

beadbybead = [];
numberPerPop = zeros(size(populations,2),1);  % MD: seems useless?
%beadbybead = concatenated (and scrambled) table
%beadperFCS = how many beads for every FCS file (barcode)
%bcnum_col = column number in beadbybead in which the original barcode
%number is specified
[beadbybead, beadperFCS, bcnum_col] = ConcatFCSfile( filepath, filename_start, filename_end, populations, scramble);


%% Step 3: GATE YOUR SINGLETS! 

% Fluorescent background of Bare beads. 
bb_median = median(bb_fcsdat,1);

% ATTENTION: the MFI of bare beads is being calculated without being
% singled out first. should gate it then calculate MFI
%bb_MFI= [bb_median(:,I1col),bb_median(:,I2col),bb_median(:,I3col),bb_median(:,I4col)]  ;   
bb_MFI=[816 1307 723 131];

% Call singlet function, specifying bacgkround, and the column number of
% relevant variables in the tables. 
[beadbybead_singlets, beadbybead_beads] = SingletGate( bb_MFI, beadbybead, var_col_nb);

%% Step 4: Find normalizaiton constants

% To-do: 
beta= [1 1 1.68 1.68];  %  laser Normalization constant found from normalizing the S-B in exp BC4p3 (on excel)


%% Step 5: DECODING

%Output is in dat_Final_Result (col1:4 = I1-I4, col5: red cluster nb,
% col6: blue cluster nb (Depends on red)
% col7-10: n1-n4 decoded
% col11: original barcode nb
% col12: decoded and assigned barcode

%added posterior_red & posterior_blue (MK)
[unique_red, Inorm_model_R, dat_Final_Result, clusters, posterior_red, posterior_blue, num_clusters_blue  ] = Decode4D( bc, beadbybead_singlets, bcnum_col, var_col_nb, beta, bb_MFI, model );

% For sorting purposes. 
clusters(:,7)=sqrt(clusters(:,3).^2 + clusters(:,4).^2);
clusters(:,8)=sqrt(clusters(:,5).^2 + clusters(:,6).^2);


%% Step6: Sort if not using model

if model==0
    clusters_sorted=sortrows(clusters, [8 7]);
   for i=1:size(populations,2)
    bcs_assigned(i,:) = [i , bc(i,:)];
   end
   bcs_assigned(:,6)= sqrt(bcs_assigned(:,2).^2+bcs_assigned(:,3).^2);
   bcs_assigned(:,7)= sqrt(bcs_assigned(:,4).^2+bcs_assigned(:,5).^2);
   bcs_assigned=sortrows(bcs_assigned, [7 6]);
   clusters_assign= [bcs_assigned(:,1), clusters_sorted(:,1:2)];
   for j= 1:size(dat_Final_Result,1)
   sluters_indeces=  find(clusters_assign(:,3) == dat_Final_Result(j,5) & clusters_assign(:,2) == dat_Final_Result(j,6));
   dat_Final_Result(j,12)=clusters_assign(sluters_indeces);
   end  
end



%% Step7: Plot HeatMap: 


Dat_Final_Result=dat_Final_Result;

n=size(populations,2);
number = zeros (n,n); 
for i=1:size(Dat_Final_Result,1)
    number(Dat_Final_Result(i,11),Dat_Final_Result(i,12))= number(Dat_Final_Result(i,11),Dat_Final_Result(i,12))+1;
end

    %calculate percentages
    for i=1:n
        beadnum=0;
        for j=1:n 
            beadnum=beadnum+number(i,j);        
        end
        for j=1:n
            percentage(i,j)=100*number(i,j)/beadnum; 
        end
    end
    
    
figure
colormap('parula');
colorbar; 
imagesc(percentage)

%% Vary the posterior proability threshold and get corresponding figures
%  (added 4/9/2017)

%[rejection_rate_vector, error_rate_vector, percentage_posterior] = BC_PosteriorProbabilityFigures(Dat_Final_Result, num_clusters_blue, posterior_blue, posterior_red);

autoArrangeFigures();

    

