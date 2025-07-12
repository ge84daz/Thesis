%__________________________________________________________________________
%
%        f_GiD_export_3D.m    
%
%       INPUT: 
%
%       OUTPUT:
%
%       DESCRIPTION: 
%       The function writes two Outputfiles:
%
%       1. filename.msh contains the nodes and elements with their coordinates
%       2. filename.res contains the displacements and stresses. This
%          skript has to be put into GiD.
%          -> Displacements at nodes: real for (x,y,z,f) and (x,y,z,t)
%          -> Stresses at GP: real for (xx,yy,zz,omega) and (xx,yy,zz,t)
%
%       Caution: Both Files have to have the same name and the same location.
%
%       - Imaginary values of stresses and displacements not exported as they are
%         zero (or maybe e-16 due to numerical errors)
%       - Absolute values of stresses and displacements not exported as the
%         absolute value u_abs = (ux^2 + u_y^2 + u_z^2)^(1/2) is automatically
%         calculated by GiD as |u|
%       - Results are calculated 2.5D for each kx independently but visualized
%         three dimensional after retaining all information with the IFFT)
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_GiD_export_3D)
%       Modified:   Tom Hicks                 
%       Date:       01-12-2021 - Freisinger 
%       Changed:    04-12-2024 - Hicks   
%__________________________________________________________________________

function [] = f_GiD_export_3D(load_dir,filename,dis_itm,dis_fem_T1,geo_T1,loading,displ,path,calc,calc_cont)

%% Initialize: Input

% Load Variables
% Tunnel 1
Knoten_T1    = dis_fem_T1.Knoten;
Elemente_T1  = dis_fem_T1.Elemente;
num_elem_T1  = dis_fem_T1.num_elem;
num_nod_T1   = dis_fem_T1.num_nod;
r0_T1        = geo_T1.radius;
y_Tc_T1      = geo_T1.y_Tc;
y_OF_T1      = geo_T1.y_OF;
H_T1         = geo_T1.H;

% General discretizazion
x         = geo_T1.x;
y         = geo_T1.y;
f         = dis_itm.f;
Nf        = dis_itm.Nf;
Nx_hs     = dis_itm.Nx_hs;

% Time discretization
switch loading.type
    case loading.harm
        timesteps = dis_itm.t_harm;
        Nt        = dis_itm.Nt_harm;
    case loading.trans
        timesteps = dis_itm.t_trans;
        Nt        = dis_itm.Nt_trans;
end

% Load and select results
% Displ hscyl ITM u(x,y,z,f)
ux_xyf_hscyl_itm_z0 = getfield(displ.hscyl,calc.system,'itm',   ['ux' load_dir '_xyf_hscyl_itm_z0']);
uy_xyf_hscyl_itm_z0 = getfield(displ.hscyl,calc.system,'itm',   ['uy' load_dir '_xyf_hscyl_itm_z0']);
uz_xyf_hscyl_itm_z0 = getfield(displ.hscyl,calc.system,'itm',   ['uz' load_dir '_xyf_hscyl_itm_z0']);
ux_xyt_hscyl_itm_z0 = getfield(displ.hscyl,calc.system,'itm_Pt',['ux' load_dir '_xyt_hscyl_itm_z0']);
uy_xyt_hscyl_itm_z0 = getfield(displ.hscyl,calc.system,'itm_Pt',['uy' load_dir '_xyt_hscyl_itm_z0']);
uz_xyt_hscyl_itm_z0 = getfield(displ.hscyl,calc.system,'itm_Pt',['uz' load_dir '_xyt_hscyl_itm_z0']);
 
% Displ hscyl FEM u(x,y,z,f)
uxPi_xyf_fem_all_T1   = getfield(displ.hscyl,calc.system,'fem','T1', ['ux' load_dir '_xyf_hscyl_fem_all']);
uyPi_xyf_fem_all_T1   = getfield(displ.hscyl,calc.system,'fem','T1', ['uy' load_dir '_xyf_hscyl_fem_all']);
uzPi_xyf_fem_all_T1   = getfield(displ.hscyl,calc.system,'fem','T1', ['uz' load_dir '_xyf_hscyl_fem_all']);

uxPi_xyt_fem_all_T1   = getfield(displ.hscyl,calc.system,'fem_Pt','T1',['ux' load_dir '_xyt_hscyl_fem_all']);
uyPi_xyt_fem_all_T1   = getfield(displ.hscyl,calc.system,'fem_Pt','T1',['uy' load_dir '_xyt_hscyl_fem_all']);
uzPi_xyt_fem_all_T1   = getfield(displ.hscyl,calc.system,'fem_Pt','T1',['uz' load_dir '_xyt_hscyl_fem_all']);

%% 1.)   Outputfile for nodes and elements

texfile = fopen([path.GiD '\' filename '.msh'],'w','n','UTF-8');

%% 1.1.) Header FEM T1: Nodes, Elements

% Dimension: 3; Element: 8-Node-Element (8node_Hexahedra)
fprintf(texfile,['MESH "FEM_T1" dimension 3 ElemType Hexahedra Nnode 8 \n' 'Coordinates \n']);

% --------------------
% Nodes FEM
% --------------------
% Coordinates of nodes out of vector 'Knoten'
% [ Node-number | x | y | z ]
num_node_FEM_T1 = num_nod_T1;

% x-, y- and z-Coordinates FEM
% Copy coordinates from Ansys for each layer from -xmax...xmax. The (y,z)
% coordinates are equal in each layer.
for ii = 1:Nx_hs
    switch calc.system
        case 'trench';   nodes_FEM_GiD_T1((ii-1)*num_node_FEM_T1+1:ii*num_node_FEM_T1,:) = [Knoten_T1(:,1) + (ii-1)*num_node_FEM_T1, x(ii)*ones(num_node_FEM_T1,1), Knoten_T1(:,2)+y_Tc_T1, Knoten_T1(:,3)];
        case 'tunnel';   nodes_FEM_GiD_T1((ii-1)*num_node_FEM_T1+1:ii*num_node_FEM_T1,:) = [Knoten_T1(:,1) + (ii-1)*num_node_FEM_T1, x(ii)*ones(num_node_FEM_T1,1), Knoten_T1(:,2)+y_Tc_T1, Knoten_T1(:,3) + H_T1];
    end
end
fprintf(texfile,'%d %f %f %f\n',nodes_FEM_GiD_T1.');
fprintf(texfile,['End Coordinates \n' 'Elements \n']);

% --------------------
% Elements FEM
% --------------------
num_elem_FEM_T1 = num_elem_T1;

% Connect each element in layer i (for one x coord) with the corresponding
% element in layer j (for x coord from before +dx) following the GiD
% element definition for 8 node Hexahedral elements (1. bottom layer
% clockwise, 2. top layer clockwise).
%
% Therefore use Element connectivities for each layer from Ansys and
% add maximal number of FEM nodes for each layer.
% Apply for all x from -xmax...+xmax.
%
% First num_elem_FEM rows for elements between [-xmax; -xmax+dx].
% Then for all x-layers until [+xmax-dx; +xmax].

for ii = 1:Nx_hs-1
    elements_FEM_GiD_T1((ii-1)*num_elem_FEM_T1+1:ii*num_elem_FEM_T1,:) = [Elemente_T1(:,1)+(ii-1)*num_elem_FEM_T1, Elemente_T1(:,5)+ii*max(max(Elemente_T1(:,2:5))), Elemente_T1(:,5)+(ii-1)*max(max(Elemente_T1(:,2:5))), Elemente_T1(:,4)+(ii-1)*max(max(Elemente_T1(:,2:5))), Elemente_T1(:,4)+ii*max(max(Elemente_T1(:,2:5))), Elemente_T1(:,2)+ii*max(max(Elemente_T1(:,2:5))), Elemente_T1(:,2)+(ii-1)*max(max(Elemente_T1(:,2:5))), Elemente_T1(:,3)+(ii-1)*max(max(Elemente_T1(:,2:5))), Elemente_T1(:,3)+ii*max(max(Elemente_T1(:,2:5)))];
end

% Elements out fo vector 'Elemente_T1'
% [ Element-number | Kn1 | Kn2 | Kn3 | Kn4 | Kn5 | Kn6 | Kn7 | Kn8]
fprintf(texfile,'%d %d %d %d %d %d %d %d %d\n',elements_FEM_GiD_T1.');
fprintf(texfile,'End Elements \n');

%% 1.2.) Header ITM: Nodes, Elements
%**************************************************************************
% Dimension: 3 (z=0);
% Element:   3-Triangluar-Element (at z=0m) 
fprintf(texfile,['MESH "ITM" dimension 3 ElemType Triangle Nnode 3 \n' 'Coordinates \n']);

% --------------------
% Nodes ITM
% --------------------
% Creating triangular elements at z=0m in ITM-area [Node-number x y z]
% y-nodes at z=0 
% - Trench: only outside trenches T1 (& T2)
% - Tunnel: all y nodes
switch calc_cont.mesh_type
    case 'half_cyl'
            num_nod_FEM_tot = nodes_FEM_GiD_T1(end,1);
            y_OF_temp = y_OF_T1;
    case 'full_cyl'
            num_nod_FEM_tot = nodes_FEM_GiD_T1(end,1);
            y_OF_temp = y;
end

% Sort y-coord of nodes [-ymax...0...+ymax]
temp1 = sortrows(y_OF_temp(y_OF_temp>=0));
temp2 = sortrows(y_OF_temp(y_OF_temp<0));
y_OF  = [temp2 temp1];
        
% Numbering of ITM-nodes continues from the total number of FEM-nodes
% -> num_nod_FEM_tot = maximum FEM node number considering all x-layers: 
num_nodes_ITM = size(y_OF,2)*Nx_hs;
numbering_ITM = (num_nod_FEM_tot+1:1:num_nod_FEM_tot+num_nodes_ITM).';
[X_OF,Y_OF]   = ndgrid(x,y_OF); % all y,z points at z=0

% Transpose -> over rows -ymax...+ymax over columns: -xmax...+xmax
X_OF = X_OF.';
Y_OF = Y_OF.';

% Row vector
X_OF = X_OF(:);
Y_OF = Y_OF(:);

nodes_x_ITM = X_OF;
nodes_y_ITM = Y_OF;
nodes_z_ITM = zeros(size(X_OF,1),1);

nodes_ITM_z0_GiD  = [numbering_ITM nodes_x_ITM nodes_y_ITM nodes_z_ITM];

% [ Node-number | x | y | z ]
fprintf(texfile,'%d %f %f %f \n',nodes_ITM_z0_GiD.');

fprintf(texfile,['End Coordinates \n' 'Elements \n']);

% --------------------
% Elements ITM
% --------------------

% Delaunay triangulation based on points at z=0
DT    = delaunayTriangulation(X_OF,Y_OF);
switch calc.system
    case {'trench'}
        % Calculate incenter points of triangles
        C    = incenter(DT);
        % Select all incenter points within halfcylinder
        % C_tr = abs(C(:,2)) < radius_T1;
        C_tr = C(:,2) > -r0_T1+y_Tc_T1 & C(:,2)< +r0_T1+y_Tc_T1;
       
        % Row numbers to delete in node connectivity
        row_num        = (1:1:size(C_tr)).';
        row_num_delete = row_num(C_tr);
        
        % Delete triangles inside cyl
        triangle_connectivity = DT.ConnectivityList;
        triangle_connectivity(row_num_delete,:) = [];
        
    case {'tunnel'}
        % no elements deleted for tunnel
        triangle_connectivity = DT.ConnectivityList;
end

% Add number of nodes of FEM mesh to adapt node connectivities to ITM nodes
% defined above (Numbering of ITM nodes continiues with FEM node numbers
triangle_connectivity = triangle_connectivity + num_nod_FEM_tot;

% Numbering of ITM-elements continues with the numbering of FEM-elements
num_elem_FE_tot = num_elem_T1;


numbering_elem_ITM  = (num_elem_FE_tot+1:1:num_elem_FE_tot+size(triangle_connectivity,1)).';
elements_ITM        = [numbering_elem_ITM triangle_connectivity];

% [ Element-number | Kn1 | Kn2 | Kn3 ]
fprintf(texfile,'%d %d %d %d\n',elements_ITM.');
fprintf(texfile,'End Elements \n');
fclose(texfile);

%% 2.)  Outputfile for displ & stresses
texfile = fopen([path.GiD '\' filename '.res'],'w','n','UTF-8');

%% 2.1.) Header Displacements
fprintf(texfile,'GiD Post Results File 1.0 \n' );

% ------------------------------------------
% Displ. ITM/FEM uz(x,y,z,f) - FREQUENCY REAL
% ------------------------------------------
for nnf = 1:Nf
    fprintf(texfile,['Result "Displacement (real)" "Frequency" ' num2str(f(nnf)) ' Vector OnNodes \n' 'Values \n']);
    
    %% FEM Tunnel 1
    % Real part of displacements in FE-area
    ux_xyf_fem_real_T1 = squeeze(real(uxPi_xyf_fem_all_T1(:,:,nnf)));
    uy_xyf_fem_real_T1 = squeeze(real(uyPi_xyf_fem_all_T1(:,:,nnf)));
    uz_xyf_fem_real_T1 = squeeze(real(uzPi_xyf_fem_all_T1(:,:,nnf)));
    ux_xyf_fem_real_T1 = ux_xyf_fem_real_T1.';
    uy_xyf_fem_real_T1 = uy_xyf_fem_real_T1.';
    uz_xyf_fem_real_T1 = uz_xyf_fem_real_T1.';
    ux_xyf_fem_real_T1 = ux_xyf_fem_real_T1(:);
    uy_xyf_fem_real_T1 = uy_xyf_fem_real_T1(:);
    uz_xyf_fem_real_T1 = uz_xyf_fem_real_T1(:);
    
    displ_fem_f_real_T1 = [nodes_FEM_GiD_T1(:,1) ux_xyf_fem_real_T1 uy_xyf_fem_real_T1 uz_xyf_fem_real_T1];
    
    fprintf(texfile,'%u %e %e %e \n',displ_fem_f_real_T1.');

   
    %% ITM Real part of displ at z=0 over (x,y)
    switch calc.system
        case {'trench'}
            ux_xyf_real = real(squeeze(ux_xyf_hscyl_itm_z0(:,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnf)));
            uy_xyf_real = real(squeeze(uy_xyf_hscyl_itm_z0(:,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnf)));
            uz_xyf_real = real(squeeze(uz_xyf_hscyl_itm_z0(:,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnf)));
        case {'tunnel'}
            ux_xyf_real = real(squeeze(uy_xyf_hscyl_itm_z0(:,:,nnf)));
            uy_xyf_real = real(squeeze(uy_xyf_hscyl_itm_z0(:,:,nnf)));
            uz_xyf_real = real(squeeze(uz_xyf_hscyl_itm_z0(:,:,nnf)));
    end
    % transpose
    ux_xyf_real = ux_xyf_real.';
    ux_xyf_real = ux_xyf_real(:);
    
    uy_xyf_real = uy_xyf_real.';
    uy_xyf_real = uy_xyf_real(:);
    
    uz_xyf_real = uz_xyf_real.';
    uz_xyf_real = uz_xyf_real(:);
    
    displ_omega_itm_real = [numbering_ITM ux_xyf_real uy_xyf_real uz_xyf_real];
    
    fprintf(texfile,'%u %e %e %e \n',displ_omega_itm_real.');
    fprintf(texfile,'End Values \n');
end

% ------------------------------------------
% Displ. ITM/FEM uz(x,y,z,t) - TIME REAL
% ------------------------------------------

for nnt = 1:Nt
    fprintf(texfile,['Result "Displacement (real)" "Time" ' num2str(timesteps(nnt)) ' Vector OnNodes \n' 'Values \n']);
    
    %% FEM Tunnel 1
    % Real part of displacements in FE-area
    ux_xyt_fem_real_T1 = squeeze(real(uxPi_xyt_fem_all_T1(:,:,nnt)));
    uy_xyt_fem_real_T1 = squeeze(real(uyPi_xyt_fem_all_T1(:,:,nnt)));
    uz_xyt_fem_real_T1 = squeeze(real(uzPi_xyt_fem_all_T1(:,:,nnt)));
    ux_xyt_fem_real_T1 = ux_xyt_fem_real_T1.';
    uy_xyt_fem_real_T1 = uy_xyt_fem_real_T1.';
    uz_xyt_fem_real_T1 = uz_xyt_fem_real_T1.';
    ux_xyt_fem_real_T1 = ux_xyt_fem_real_T1(:);
    uy_xyt_fem_real_T1 = uy_xyt_fem_real_T1(:);
    uz_xyt_fem_real_T1 = uz_xyt_fem_real_T1(:);
    
    displ_fem_f_real_T1 = [nodes_FEM_GiD_T1(:,1)  ux_xyt_fem_real_T1  uy_xyt_fem_real_T1  uz_xyt_fem_real_T1];
    
    fprintf(texfile,'%u %e %e %e \n',displ_fem_f_real_T1.');
    
    
    %% ITM Real part of displ at z=0 over (x,y)
    switch calc.system
        case {'trench'}
            ux_xyt_real = real(squeeze(ux_xyt_hscyl_itm_z0(:,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnt)));
            uy_xyt_real = real(squeeze(uy_xyt_hscyl_itm_z0(:,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnt)));
            uz_xyt_real = real(squeeze(uz_xyt_hscyl_itm_z0(:,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnt)));
        case {'tunnel'}
            ux_xyt_real = real(squeeze(ux_xyt_hscyl_itm_z0(:,:,nnt)));
            uy_xyt_real = real(squeeze(uy_xyt_hscyl_itm_z0(:,:,nnt)));
            uz_xyt_real = real(squeeze(uz_xyt_hscyl_itm_z0(:,:,nnt)));
    end
    % transpose
    ux_xyt_real = ux_xyt_real.';
    ux_xyt_real = ux_xyt_real(:);
    uy_xyt_real = uy_xyt_real.';
    uy_xyt_real = uy_xyt_real(:);
    uz_xyt_real = uz_xyt_real.';
    uz_xyt_real = uz_xyt_real(:);
    
    displ_omega_itm_real = [numbering_ITM ux_xyt_real uy_xyt_real uz_xyt_real];
    
    fprintf(texfile,'%u %e %e %e \n',displ_omega_itm_real.');
    
    fprintf(texfile,'End Values \n');
    
end
fclose(texfile);

end
