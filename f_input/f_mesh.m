%__________________________________________________________________________
%
%       f_mesh.m    
%
%       INPUT: 
%       path          = path parameters
%       calc          = string defining which material to be used
%
%       OUTPUT: 
%       - none
%
%       DESCRIPTION: 
%       - Create FE Mesh within Tunnel with ANSYS (4 Node Elements)
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_mesh)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [] = f_mesh(path_aoad,calc_cont,calc)

% Initialize: Input
if strcmp(calc_cont.mesh,'new')
    % Call Ansys in Batch Mode for mesh creation
    % Save current path
    % ATTENTION: Maximum 32 characters (= max Ansys APDL string length)
    directory = cd;
    
    % Identify current installed Ansys Version
    Ansys_ordner = dir('C:\Program Files\ANSYS Inc');
    Anzahl_namen = length(Ansys_ordner);
    for j=1:1:Anzahl_namen
        name = double(Ansys_ordner(j).name);
        if name(1)==118  % 118 im double-Format entspricht 'v' im char-Format
            version = char(name(2:4));
        end
    end
    % Set Version of Ansys manually due to compatibility reasons
    % -> because Ansys 19.2 export different result files not working with this
    % code in the current verstion
    % version='181';
    
    % Get number of logical cores (Kerne)
    numcores = num2str(feature('numcores'));
    
    %% Write Batch File
    
    fid = fopen([directory '/ansys_mesh.bat'],'wt');    % Open file in textmode with permisson 'write'(wt)
    fprintf(fid,'%s',[...                               % Write a line in Format 'string' (%s)
        '"C:\Program Files\Ansys Inc\v' version '\ANSYS\bin\winx64\ansys' version '.exe" -p ' calc_cont.ansys_lic ' -dir '...
        '"' path_aoad '"' ' -j "Netz" -s noread -l en-us -np ' numcores ' -b -i '...
        '"' directory '\' calc.mesh_ansys '.mac" -o '...% ANSYS APDL Skript that controls the meshing!!!
        '"' path_aoad '\output_ansys.out"'
        ]);
    fclose(fid);
    
    % Commands for batch running
    % -p    Product aa_r = research, aa_t_a = teaching
    % -dir  working directory
    % -j    Jobname
    % -np   number of used cores
    % -b    causes input listing to be included in the output
    % -i    Specifies the name of the file to read input into Mechanical APDL for batch processing
    % -o    Specifies the name of the file to store the output from a batch execution of Mechanical APDL
    
    
    %% Execute Batch file
    % evalc(['!ansys_netz']);
    % Run Batch file
    system('SET KMP_STACKSIZE=2048k & ansys_mesh.bat');
    
    %% Check outputfile for errors
    
    fid = fopen([path_aoad '\output_ansys.out']);
    l = fgetl(fid);
    errornumber=0;
    while ~feof(fid)
        l = fgetl(fid);
        if ~isempty(strfind(l, '*** ERROR ***'))
            errornumber=errornumber+1;
            errordata=textscan(fid,'%s',3,'delimiter',',,');
            disp(['Fehler' num2str(errornumber) ':'])
            disp(errordata{1})
        end
    end
    fclose(fid);
    if errornumber~=0
        return;
    end
end
end