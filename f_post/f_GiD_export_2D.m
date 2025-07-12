%__________________________________________________________________________
%
%        f_GiD_export_2D.m    
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
%       - Export of results and mesh in cutting plane at x=0 over (y,z)
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_GiD_export_2D)
%       Modified:   Tom Hicks                 
%       Date:       01-12-2021 - Freisinger 
%       Changed:    04-12-2024 - Hicks   
%__________________________________________________________________________

function [] = f_GiD_export_2D(load_dir,filename,dis_itm,dis_fem_T1,geo_T1,loading,displ,path,calc,calc_cont)

%% 0.)   Initialize: Input

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
H_cyl_tot_T1 = geo_T1.H_cyl_tot;

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

%% 1.1.) Header FEM T1: Nodes, Elements Tunnel 
% Dimension: 3 (Ux=0); Element: 4-Node-Element "ShellElement3d4n"
fprintf(texfile,['MESH "FEM_T1" dimension 3 ElemType Quadrilateral Nnode 4 \n' 'Coordinates \n']);

% --------------------
% Nodes FEM
% --------------------
% Coordinates of nodes out of vector 'Knoten'
% [ Node-number | x | y | z ]
Knoten_x_T1 = zeros(num_nod_T1,1);

% y- and z-Coordinates FEM
switch calc.system
    case {'trench'} 
        nodes_FEM_GiD_T1 = [Knoten_T1(:,1), Knoten_x_T1, Knoten_T1(:,2)+y_Tc_T1, Knoten_T1(:,3)];
    case 'tunnel'
        nodes_FEM_GiD_T1 = [Knoten_T1(:,1), Knoten_x_T1, Knoten_T1(:,2)+y_Tc_T1, Knoten_T1(:,3) + H_T1];   
end 
    
fprintf(texfile,'%d %f %f %f\n',nodes_FEM_GiD_T1.');

fprintf(texfile,['End Coordinates \n' 'Elements \n']);

% --------------------
% Elements FEM
% --------------------
% Elements out of vector 'Elemente'
% [ Element-number | Kn1 | Kn2 | Kn3 | Kn4]

fprintf(texfile,'%d %d %d %d %d\n',Elemente_T1.');

fprintf(texfile,'End Elements \n');

%% 1.2.) Header ITM: Nodes, Elements

% Dimension: 3 (Ux=0);
% Element: 2-Node-Element (line elements on ground at z=0m) LineElement3d2n_Mesh_2

fprintf(texfile,['MESH "ITM" dimension 3 ElemType Line Nnode 2 \n' 'Coordinates \n']);

% --------------------
% Nodes ITM
% --------------------
% y-discrtization at x=z=0
% Trench: only outside trench area of T1 (& T2)
% Tunnel: all surface nodes
switch calc_cont.mesh_type
    case 'half_cyl'
          y_OF = y_OF_T1;
    case 'full_cyl'          
          y_OF = y;
end

% Number of ITM discretization nodes at z=0
num_nodes_ITM = size(y_OF,2);

% Numbering of ITM-nodes continues from total number of FEM-nodes
numbering_ITM = (num_nod_T1+1:1:num_nod_T1+num_nodes_ITM).';

% x coord = 0 for all nodes
% z coord = 0 for all nodes on ground surface
nodes_x_ITM   = zeros(num_nodes_ITM,1);
nodes_z_ITM   = zeros(num_nodes_ITM,1);

% Sort y-coord of nodes [-ymax...0...+ymax]
temp1 = sortrows(y_OF(y_OF>=0));
temp2 = sortrows(y_OF(y_OF<0));
nodes_y_ITM = [temp2 temp1];

% ITM nodes at z=0m in ITM-area [Node-number x y z]
nodes_ITM_z0_GiD = [numbering_ITM nodes_x_ITM nodes_y_ITM.' nodes_z_ITM];

% [ Node-number | x | y | z ]
fprintf(texfile,'%d %f %f %f \n',nodes_ITM_z0_GiD.');

fprintf(texfile,['End Coordinates \n' 'Elements \n']);

% --------------------
% Elements ITM
% --------------------
% Creating line-elements at z=0m in ITM-area [Element-number Kn1 Kn2]
% Numbering of ITM-elements continues from total number of FEM-elements
num_elem_FEM_tot = num_elem_T1;

switch calc_cont.mesh_type
    case 'half_cyl'
            numbering_elem_ITM  = (num_elem_FEM_tot+1:1:num_elem_FEM_tot+num_nodes_ITM-2).'; % no element within T1  
    case 'full_cyl'                
            numbering_elem_ITM  = (num_elem_FEM_tot+1:1:num_elem_FEM_tot+num_nodes_ITM-1)';
end

% Sorting ITM-nodenumbers, that belong to the ITM-element numbers
elem_y_z_ITM = [numbering_ITM(1):1:numbering_ITM(1)+num_nodes_ITM-2; numbering_ITM(1)+1:1:numbering_ITM(1)+num_nodes_ITM-1].';

% find id of elements connecting nodes left and right of T1 (and T2) and delet them
if strcmp(calc_cont.mesh_type,'half_cyl')
    id = find(nodes_ITM_z0_GiD(:,3)==y_Tc_T1-r0_T1);
    elem_y_z_ITM(id,:) = [];
end

elements_ITM = [numbering_elem_ITM elem_y_z_ITM];

% [ Element-number | Kn1 | Kn2 ]
fprintf(texfile,'%d %d %d \n',elements_ITM.');
fprintf(texfile,'End Elements \n');
fclose(texfile);

%% 2.)   Outputfile for displ & stresses     

texfile = fopen([path.GiD '\' filename '.res'],'w','n','UTF-8');

%% 2.1.) Header Displacements
fprintf(texfile,'GiD Post Results File 1.0 \n' );

% ------------------------------------------
% Displ. ITM/FEM uz(x,y,z,f) - FREQUENCY REAL
% ------------------------------------------
for nnf = 1:Nf
    fprintf(texfile,['Result "Displacement (real)" "Frequency" ' num2str(f(nnf)) ' Vector OnNodes \n' 'Values \n']);
    
    %% Tunnel T1
    % Real part of displacements in FE-area at x=0
    ux_x0yzf_fem_re_T1  = squeeze(real(uxPi_xyf_fem_all_T1(x==0,:,nnf)));
    uy_x0yzf_fem_re_T1  = squeeze(real(uyPi_xyf_fem_all_T1(x==0,:,nnf)));
    uz_x0yzf_fem_re_T1  = squeeze(real(uzPi_xyf_fem_all_T1(x==0,:,nnf)));
    
    displ_fem_f_re_T1 = [nodes_FEM_GiD_T1(:,1) ux_x0yzf_fem_re_T1.' uy_x0yzf_fem_re_T1.' uz_x0yzf_fem_re_T1.'];
    fprintf(texfile,'%u %e %e %e \n',displ_fem_f_re_T1.');

    
    %% ITM
    % Real part of displacements in ITM-area at x=z=0 over y
    switch calc.system
        case {'trench'}
            ux_x0yz0_itm = squeeze(real(ux_xyf_hscyl_itm_z0(x==0,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnf)));
            uy_x0yz0_itm = squeeze(real(uy_xyf_hscyl_itm_z0(x==0,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnf))); 
            uz_x0yz0_itm = squeeze(real(uz_xyf_hscyl_itm_z0(x==0,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnf)));
        case {'tunnel'}
            ux_x0yz0_itm = squeeze(real(ux_xyf_hscyl_itm_z0(x==0,:,nnf)));
            uy_x0yz0_itm = squeeze(real(uy_xyf_hscyl_itm_z0(x==0,:,nnf)));  
            uz_x0yz0_itm = squeeze(real(uz_xyf_hscyl_itm_z0(x==0,:,nnf)));
    end
    
    displ_itm_f_re  = [numbering_ITM ux_x0yz0_itm.' uy_x0yz0_itm.' uz_x0yz0_itm.'];
    
    fprintf(texfile,'%u %e %e %e \n',displ_itm_f_re.');
    
    fprintf(texfile,'End Values \n');
    
end
clear ii

% ------------------------------------------
% Displ. ITM/FEM uz(x,y,z,t) - TIME REAL
% ------------------------------------------
for nnt = 1:Nt
    fprintf(texfile,['Result "Displacement (real)" "Time" ' num2str(timesteps(nnt)) ' Vector OnNodes \n' 'Values \n']);
    
    %% FEM Tunnel T1
    % Real part of displacements in FE-area at x=0
    ux_x0yzt_fem_re_T1 = squeeze(real(uxPi_xyt_fem_all_T1(x==0,:,nnt)));
    uy_x0yzt_fem_re_T1 = squeeze(real(uyPi_xyt_fem_all_T1(x==0,:,nnt)));
    uz_x0yzt_fem_re_T1 = squeeze(real(uzPi_xyt_fem_all_T1(x==0,:,nnt)));
    
    displ_fem_t_re_T1  = [nodes_FEM_GiD_T1(:,1) ux_x0yzt_fem_re_T1.' uy_x0yzt_fem_re_T1.' uz_x0yzt_fem_re_T1.'];
    fprintf(texfile,'%u %e %e %e \n',displ_fem_t_re_T1.');
    
    %% ITM
    % Real part of displacements in ITM-area at x=z=0 over y
    switch calc.system
        case {'trench'}
            ux_x0yz0_itm = squeeze(real(uy_xyt_hscyl_itm_z0(x==0,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnt)));
            uy_x0yz0_itm = squeeze(real(uy_xyt_hscyl_itm_z0(x==0,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnt)));              
            uz_x0yz0_itm = squeeze(real(uz_xyt_hscyl_itm_z0(x==0,y<=+y_Tc_T1-r0_T1 | y>=y_Tc_T1+r0_T1,nnt)));  
        case {'tunnel'}
            ux_x0yz0_itm = squeeze(real(ux_xyt_hscyl_itm_z0(x==0,:,nnt))); 
            uy_x0yz0_itm = squeeze(real(uy_xyt_hscyl_itm_z0(x==0,:,nnt)));  
            uz_x0yz0_itm = squeeze(real(uz_xyt_hscyl_itm_z0(x==0,:,nnt)));
    end
    
    displ_itm_t_re  = [numbering_ITM ux_x0yz0_itm.' uy_x0yz0_itm.' uz_x0yz0_itm.'];
    
    fprintf(texfile,'%u %e %e %e \n',displ_itm_t_re.');
    fprintf(texfile,'End Values \n');
    
end
clear ii
fclose(texfile);
        
end
