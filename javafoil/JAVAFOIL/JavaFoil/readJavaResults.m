function [alfa, cl, cd, cm, cp] = readJavaResults(process_id)
    % Inputs:
    % process_id : [1 x 1] (str) ID of current task (process) - suffix for W/R files
            
    if ~ischar(process_id)
        disp('readJavaResults: argument process_id is not valid.\n');
        return;
    end

    
    % Open results file from JavaFoil
    filename = [pwd,'/JavaFoil/output_',process_id,'.txt'];
    fid      = fopen(filename,'r');
    Intro    = textscan(fid,'%s',5,'Delimiter','\n');
    
    data = fscanf(fid,'%f',[11 Inf])'; %  Read the data
    fclose(fid);
    
    % Extract coefficients
    alfa = deg2rad(data(:,1)); 
    cl   = data(:,2);
    cd   = data(:,3);
    cm   = data(:,4);
    cp   = data(:,11);

end