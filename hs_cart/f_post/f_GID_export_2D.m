%__________________________________________________________________________
%
%       f_GID_export_2D.m    
%
%       INPUT: 
%       path      = pathes
%       calc      = calculation setup
%       plot_cont = plot control parameters
%       loading   = load amplitude, width, position
%       dis_itm   = discretization 
%       geo       = geometry

%       OUTPUT: 
%
%       DESCRIPTION: 
%       The function writes two Outputfiles:
%
%       1. GiD_Output.post_2D.msh contains the nodes and elements with their
%       coordinates
%
%       2. GiD_Output.post_2D.res contains the displacements and stresses. This
%       skript has to be put into GiD.
%      -> Displacements at nodes: real for (x,y,z,omega) and (x,y,z,t)
%      -> Stresses at GP: real for (xx,yy,zz,omega) and (xx,yy,zz,t)
% Caution: Both Files have to have the same name and the same location.
%
%
%       REMARK: 
%       Original Author: Julian Freisinger  (f_GID_export_2D)
%       Modified: Tom Hicks               
%       Date:    23-05-2023 
%       Changed: 23-05-2023 - Hicks   
%__________________________________________________________________________



function [] = f_GID_export_2D(path,loading,dis_itm,x,y,uz_xyt_z0,name,n_node)

%% Initialize: Input

% Just use ever n_node node for visualisation 
timesteps  = dis_itm.t_harm;


%% Outputfile for nodes 
texfile=fopen([path.GiD '\GiD_' name '.post_2D.msh'],'w','n','UTF-8');

% Dimension: 3 (z=0);
% Element:   3-Triangluar-Element (at z=0m)
fprintf(texfile,['MESH "Triangular_surf_elements" dimension 3 ElemType Triangle Nnode 3 \n' 'Coordinates \n']);



[X_OF,Y_OF]     = ndgrid(x(1:n_node:end),y(1:n_node:end)); % all x,y points at z=0

num_nodes_ITM   = length(x(1:n_node:end))*length(y(1:n_node:end));
numbering_ITM   = (1:1:num_nodes_ITM).';
% Transpose -> over rows -ymax...+ymax over columns: -xmax...+xmax
X_OF = X_OF.';
Y_OF = Y_OF.';

% Row vector
X_OF = X_OF(:);
Y_OF = Y_OF(:);

nodes_x_ITM = X_OF;
nodes_y_ITM = Y_OF;
nodes_z_ITM = zeros(size(X_OF,1),1);

nodes_ITM_z0_GiD   =[numbering_ITM nodes_x_ITM nodes_y_ITM nodes_z_ITM];

% [ Node-number | x | y | z ]
fprintf(texfile,'%d %f %f %f \n',nodes_ITM_z0_GiD.');

fprintf(texfile,['End Coordinates \n' 'Elements \n']);

%% Outputfile for elements

DT          = delaunayTriangulation(X_OF,Y_OF);
tri_connect = DT.ConnectivityList;

% Numbering of ITM-elements continues with the numbering of FEM-elements
numbering_elem_ITM  = (1:1:size(tri_connect,1)).';

elements_ITM = [numbering_elem_ITM tri_connect];

% [ Element-number | Kn1 | Kn2 | Kn3 ]
fprintf(texfile,'%d %d %d %d\n',elements_ITM.');

fprintf(texfile,['End Elements \n']);

fclose(texfile);

%% Outputfile for displacements

texfile = fopen([path.GiD '\GiD_' name '.post_2D.res'],'w','n','UTF-8');

fprintf(texfile,['GiD Post Results File 1.0 \n' ]);

NNt = length(timesteps);

% Results in (x,y,z=0,t)
for ii = 1:NNt
    
    fprintf(texfile,['Result "Displacement (real)" "Time" ' num2str(timesteps(ii)) ' Vector OnNodes \n' 'Values \n']);
    
    u_x_t_real = real(squeeze(zeros(size(uz_xyt_z0(1:n_node:end,1:n_node:end,ii)))));
    u_y_t_real = real(squeeze(zeros(size(uz_xyt_z0(1:n_node:end,1:n_node:end,ii)))));
    u_z_t_real = real(squeeze(uz_xyt_z0(1:n_node:end,1:n_node:end,ii)));
    
    u_x_t_real = u_x_t_real.';
    u_x_t_real = u_x_t_real(:);
    
    u_y_t_real = u_y_t_real.';
    u_y_t_real = u_y_t_real(:);
    
    u_z_t_real = u_z_t_real.';
    u_z_t_real = u_z_t_real(:);
    
    displ_t_itm_real = [numbering_ITM u_x_t_real u_y_t_real u_z_t_real];
    
    fprintf(texfile,'%u %e %e %e \n',displ_t_itm_real.');

    fprintf(texfile,['End Values \n']);
end
clear ii

fclose(texfile);

end


