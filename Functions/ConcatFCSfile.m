function [ BeadbyBead, BeadperFCS, BCnum_col, Minimum_bead_number] = ConcatFCSfile( Filepath, Filename_start, Filename_end, Populations, Scramble )
%ConcatFCS reads and concatenates FCS files, saving their original file as
%an extra column, with the option of scrambling the data too. 
% Returns:
% BeadbyBead is the table of all the bead values (all detectors) + the
% assigned barcode number as per "Populations"
% BeadperFCS lists how many beads per FCS files
% BCnum_col is the column number at which the BC number is saved at
% numberPerPop is amount of beads in each digitally concatenated file,
% necessary for posterior probability graph (MK)
    % within Bead_Table
    
    
BeadbyBead = [];
BeadperFCS = [];


% Minimum bead number for all populations
Minimum_bead_number=1e6;

for i=1:size(Populations,2)
        [fcsdat,fcshdr,scaled]=fca_readfcs(strcat(Filepath,Filename_start, num2str(Populations(i)), Filename_end));
        Beadnumber=size(fcsdat,1);
        if Beadnumber<Minimum_bead_number
        Minimum_bead_number=Beadnumber;
        end
end

% MD: the number of events we extract from every FCS file is equal to the number of events
% in the FCS with the LEAST events?
% may have been my mistake but maybe not too good a measure (what if one
% FCS file didn't have a lot of events..)

for i=1: size(Populations,2)
    
    % Read FCS file i  
    [fcsdat,fcshdr,scaled]=fca_readfcs(strcat(Filepath,Filename_start, num2str(Populations(i)), Filename_end));
    
    % Find number of beads in this FCS file and create an array with that
    % length, with content = BC number (that is, 'i')
    pop_assign = zeros(size(fcsdat(1:Minimum_bead_number,7),1),1);
    
    %pop_assign(:) = Populations(i);
    pop_assign(:) = i;
    
    % Extract the valuable columns from FCS files, and append the BC number
    %data_per_pop = [fcsdat(:,5), fcsdat(:,7), fcsdat(:,13), fcsdat(:, 15), pop_assign]; 
    data_per_pop = [fcsdat(1:Minimum_bead_number,:) , pop_assign]; 
   
    % BCnum_col is the column number at which the BC number is saved at
    % within Bead_Table
    BCnum_col = size(fcsdat,2) + 1 ; 

    % How many beads per FCS file? 
    BeadperFCS(i) = size(data_per_pop,1);  
    
    % Concatenate! 
    BeadbyBead = [BeadbyBead ; data_per_pop];
    
    % 
    %Gating_Table = [fcsdat(:,1), fcsdat(:,2), fcsdat(:,3)];
    %Gating_Table_Final = [Gating_Table_Final ; Gating_Table];
end

    % Scramble? 
    if Scramble==1 
        BeadbyBead= BeadbyBead(randperm(size(BeadbyBead,1)),:);
    end

end

