function [ FSC_cols ] = FCScol_info_Canto(ColNames)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here




Lookfor= {'FSC-A', 'FSC-H', 'SSC-A', 'SSC-H', 'Alexa Fluor 488-A', 'PE-A', 'APC-A', 'APC-Cy7-A'};

%ColNames= extractfield(BB_fcshdr_par, 'name');
colsize=size(ColNames,2); 

for i=1:colsize
    if isequal (ColNames(i), Lookfor(1))
    FSCcol_A=i; 
    end
    
    if isequal (ColNames(i), Lookfor(2))
    FSCcol_H=i; 
    end
    
    if isequal (ColNames(i), Lookfor(3))
    SSCcol_A=i; 
    end 
    
    if isequal (ColNames(i), Lookfor(4))
    SSCcol_H=i;
    end
    
    if isequal (ColNames(i), Lookfor(5))
    I1col=i;
    end
    
    if isequal (ColNames(i), Lookfor(6))
    I2col=i;
    end
    
    if isequal (ColNames(i), Lookfor(7))
    I3col=i;
    end   

    if isequal (ColNames(i), Lookfor(8))
    I4col=i;
    end 
end


FSC_cols = [ FSCcol_A, FSCcol_H, SSCcol_A, SSCcol_H, I1col, I2col, I3col, I4col];

end

