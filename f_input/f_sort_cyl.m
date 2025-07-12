%__________________________________________________________________________
%
%       f_sort_cyl.m    
%
%       INPUT: 
%       path.aoad     = Path of Ansys folder
%       geo.radius    = Tunnel radius
%
%       OUTPUT: 
%       NDOF          = Number of DOFs
%       Knoten        = Nodes with Coordinates |Node ID| Ny(i)| Nz(i)|
%       num_elem      = Number of elements
%       num_nod       = Number of nodes
%       elem_coord    = Elements with coordinates
%       elem_found    = Elements of foundation
%       Edof          = Mapping Matrix
%       phi           = circumferential angle
%       num_nod_bd    = Number of nodes on boundary
%       num_nod_in    = Number of nodes inside Tunnel
%       nod_midline   = Nodes on Midline in y direction sorted from -radius+(radius/2/div) ... radius-(radius/2/div)
%
%       DESCRIPTION: 
%       - Import of Nodes and elements from ansys mesh
%       - Sorting of Nodes and Elements on boundary and inside cylinder
%       - Node and element connectivities
%       - element list for different fem structures
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_sort_cyl)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [dis_fem] = f_sort_cyl(geo,path_aoad,dis_fem,calc_cont,calc)

radius = geo.radius;

%% Import Ansys export file into Matlab
% ALL NODES and ELEMENTS:

% Import Node List
% |Node ID | Nx(i) | Ny(i) | Nz(i) |
nodes_import  = importdata([path_aoad '/Knoten_import.dat']);
num_nod = size(nodes_import,1);
% Cancel Nx(i) = x-coordinate due to plane strain problem all x=0
nodes_import(:,2) = [];

% Import Element List
% |Element ID | Nelem(i,1) | Nelem(i,2) | Nelem(i,3) | Nelem(i,4) |
% Node IDs of Element Nodes ANSYS
Elemente_import = importdata([path_aoad '/Elemente_import.dat']);


%% Sorting of nodes 
% 1. Nodes on the FEM boundary
% 2. Nodes inside FEM domain
t = 1; % set counter
r = 1; % set counter
epsilon = 0.001; % tolerance for radius

for s=1:num_nod
    if sqrt(nodes_import(s,2).^2+nodes_import(s,3).^2)>=radius-epsilon && sqrt(nodes_import(s,2).^2+nodes_import(s,3).^2)<=radius+epsilon
        nod_bd_us(t,:) = nodes_import(s,:);   % Boundary nodes unsorted
        t=t+1;
    else
        nod_in(r,:) = nodes_import(s,:);      % Nodes inside domain
        r=r+1;
    end
end
% check if all nodes included
if size(nodes_import,1) ~= size(nod_bd_us,1)+size(nod_in,1)
    error('Error when Sorting Nodes. Not all nodes included')
end
clear t r epsilon s


% Sort Nodes in dependency of angle phi
% no change of Coordinate System -> Cartesian
phi_us              = f_calc_phi(nod_bd_us, radius);   % Calculate angle phi of nodes on boundary
nod_bd_us_phi       = [nod_bd_us phi_us];              % Add column with phi
nod_bd_phi          = sortrows(nod_bd_us_phi,4);       % Sort phi 0... 2pi
nod_bd              = nod_bd_phi(:,1:3);               % Neglect column with phi
phi                 = nod_bd_phi(:,4);                 % Save angels phi

switch calc_cont.mesh_type
    case 'half_cyl'     
        % Increase the number of boundary nodes as if it was a full circle
        % Number of boundary nodes for full circle
        num_nod_bd = (size(nod_bd,1)-1)*2;
        
        % Number of nodes of the first quarter of the circle (beginning with phi=0)
        Nphi_quart = num_nod_bd/4+1;
        
    case 'full_cyl'
        % Number of boundary nodes
        num_nod_bd = size(nod_bd,1);  
end

% Number of nodes inside domain
num_nod_in          = size(nod_in,1);                

% Create vector with
% 1. Nodes on the FEM boundary
% 2. Nodes inside FEM domain
Knoten_sort = [nod_bd; nod_in];

% Rename Nodes
% Node ID of sorted Nodes renamed starting by one [1,2,...]
Knoten      = Knoten_sort;
Knoten(:,1) = 1:num_nod;

% Mapping Matrix
% Old Node ID ANSYS - NEW Node ID
mapping = [Knoten_sort(:,1) Knoten(:,1)];

% Delete nodes on the intersection point of midline and tunnel boundary
% Vector with nodes and coordinates on center line (z=0) sorted from -radius+(radius/2/div) ... radius-(radius/2/div)
% Konten |Node ID| Ny(i)| Nz(i)|
nod_midline = sortrows(Knoten((Knoten(:,3)==0),:),2);
    
% Number of nodes on center line
num_nod_midline = size(nod_midline,1);

% delete nodes on tunnel boundary - only nodes within tunnel remain
nod_midline(num_nod_midline,:) = [];
nod_midline(1,:) = [];
num_nod_midline = num_nod_midline-2;

%% Sorting of elements
% Number of Elements
num_elem = size(Elemente_import,1);

Elemente      = zeros(num_elem,5);            % prestore new Element Matrix
Elemente(:,1) = Elemente_import(:,1);         % Element IDs 1st column

% Replace old Node IDs in Element List by new Node IDs after sorting
% |Element ID | Nelem(i,1) | Nelem(i,2) | Nelem(i,3) | Nelem(i,4) |
% Node IDs of Element Nodes after sorting
for r=1:num_elem % loop over rows
    for s=2:5    % loop over columns
        for t=1:num_nod
            if Elemente_import(r,s) == mapping(t,1)
               Elemente(r,s) = mapping(t,2);
            end
        end
    end
end

% Elements Matrix with Coordinates
% |Element ID | Nelem(i,1) |Ny(i)|Nz(i)| Nelem(i,2) |Ny(i)|Nz(i)| Nelem(i,3) |Ny(i)|Nz(i)| Nelem(i,4) |Ny(i)|Nz(i)|

elem_coord = [Elemente(:,1) Elemente(:,2) Knoten(Elemente(:,2),2) Knoten(Elemente(:,2),3) ...
              Elemente(:,3) Knoten(Elemente(:,3),2) Knoten(Elemente(:,3),3) ...
              Elemente(:,4) Knoten(Elemente(:,4),2) Knoten(Elemente(:,4),3) ...
              Elemente(:,5) Knoten(Elemente(:,5),2) Knoten(Elemente(:,5),3)];

% DOFs at each Node
% u,v,w each in x,y,z direction
switch calc_cont.mesh_type
    case 'half_cyl'
        
        % Skipping the non-existing DOF of upper halfcircle for each direction
        DOF_hc_x = [1:3:3*Nphi_quart-2 (3*Nphi_quart-2)*3-2:3:3*(num_nod+(num_nod_bd/2)-1)-2].';
        DOF_hc_y = [2:3:3*Nphi_quart-1 (3*Nphi_quart-2)*3-1:3:3*(num_nod+(num_nod_bd/2)-1)-1].';
        DOF_hc_z = [3:3:3*Nphi_quart   (3*Nphi_quart-2)*3-0:3:3*(num_nod+(num_nod_bd/2)-1)-0].';
        DOFS     = [Knoten(:,1) DOF_hc_x DOF_hc_y DOF_hc_z];
        
        % Number of DOF for the halfcircle
        NDOF = 3*(num_nod+(num_nod_bd/2)-1);
        
    case 'full_cyl'
        DOFS = [Knoten(:,1) (1:3:3*num_nod).' (2:3:3*num_nod).' (3:3:3*num_nod).'];
        NDOF = (size(DOFS,2)-1)*size(DOFS,1);
        
end


% Edof Matrix
% |Element ID    | Nelem(i,1)  |  Nelem(i,2)  |  Nelem(i,3)  |  Nelem(i,4)|
%                |  u | v | w  |   u | v | w  |   u | v | w  |  u | v | w |

Edof = [Elemente(:,1) DOFS(Elemente(:,2),2) DOFS(Elemente(:,2),3) DOFS(Elemente(:,2),4)...
                      DOFS(Elemente(:,3),2) DOFS(Elemente(:,3),3) DOFS(Elemente(:,3),4)...
                      DOFS(Elemente(:,4),2) DOFS(Elemente(:,4),3) DOFS(Elemente(:,4),4)...
                      DOFS(Elemente(:,5),2) DOFS(Elemente(:,5),3) DOFS(Elemente(:,5),4)];

%%  Output
dis_fem.NDOF               = NDOF;
dis_fem.DOFS               = DOFS;

dis_fem.Knoten             = Knoten;
dis_fem.Elemente           = Elemente;
dis_fem.num_elem           = num_elem;
dis_fem.num_nod            = num_nod;
dis_fem.elem_coord         = elem_coord;

dis_fem.Edof               = Edof;
dis_fem.phi                = phi;
dis_fem.num_nod_bd         = num_nod_bd;
dis_fem.num_nod_in         = num_nod_in;
dis_fem.num_nod_midline    = num_nod_midline;
dis_fem.nod_midline        = nod_midline;

switch calc_cont.mesh_type
    case 'half_cyl'
        dis_fem.DOF_hc_x   = DOF_hc_x; % Output only if existing
        dis_fem.DOF_hc_y   = DOF_hc_y; % Output only if existing
        dis_fem.DOF_hc_z   = DOF_hc_z; % Output only if existing
    otherwise
end

end % function