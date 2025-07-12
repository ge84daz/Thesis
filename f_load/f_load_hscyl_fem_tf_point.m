%__________________________________________________________________________
%
%       f_load_hsycyl_itm_tf_point.m    
%
%       INPUT: 
%       - path          = pathes
%       - plot_cont     = plot control parameter
%       - calc_cont     = calculation control paramter
%       - loading       = load parameters
%       - dis_itm       = discretization ITM subsystem
%       - dis_fem       = discretization FEM subsystem
%       - geo           = geometry
%
%       Input arguments:
%
%       plot_load         = Parameter for Plots yes or no
%       P0_FE             = Load Amplitude
%       T                 = Time Period
%       num_nod           = Number of FE Nodes
%       Knoten            = Nodes (ID, Ny(i), Nz(i))
%       pos_load_x_fem    = Position of load midpoint
%       pos_load_y_fem    = Position of load midpoint
%       dx, dt            = Inkrements
%       Nx_hs             = Number of Points on halfspace surface
%       Bx                = Width of ground surface in x direction
%       NDOF              = Number of DOFs
%       radius            = Tunnel radius
%       div               = Number of elements per half of radius
%       x, t              = Geometry and time
%       omega             = frequency vector
%       
% 
%       OUTPUT: 
%       - P_FE          = FEM nodal load matrix in transformed domain for time harmonic Point laoding.
%
%       DESCRIPTION: 
%       The function determines the FEM nodal load vector for a time harmonic POINT load
%       within the Tunnel. A Point load in (x,y,z) direction on a horizontal line through
%       the tunnel midpoint for an arbitrary Point (x,y) within the tunnel can be applied.
%       The loadvector is transformed in the (kx,omega) domain.
%       The resultant load matrix contains:
%       Over the rows variation from -kmax... +kmax.
%       Over the columns variation from Node 1 (Px,Py,Pz)... Node n (Px,Py,Pz)
%       The 3rd dimension the frequency -fmax... 0
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_load_hscyl_itm_tf_point)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________


function [loading] = f_load_hscyl_fem_tf_point(path,plot_cont,calc_cont,loading,dis_itm,dis_fem,geo)

%% Initialize: Input

P0_FE           = loading.P0_FE;
pos_load_x_fem  = loading.pos_load_x_fem;
pos_load_y_fem  = loading.pos_load_y_fem;
load_dir        = loading.load_dir;

Knoten          = dis_fem.Knoten;
NDOF            = dis_fem.NDOF;
num_nod_midline = dis_fem.num_nod_midline;
nod_midline     = dis_fem.nod_midline;
dy_FE           = dis_fem.dy_FE;

dx              = dis_itm.dx;
kx              = dis_itm.kx;
Nx_hs           = dis_itm.Nx_hs;
Nphi_t          = dis_itm.Nphi_t;
f               = dis_itm.f;
Nf              = dis_itm.Nf;

Bx              = geo.Bx;
x               = geo.x;
H               = geo.H;
y_LR            = geo.y_LR;
z_LR            = geo.z_LR;

%% Creation of load matrix

if P0_FE > 0
    
    %% Predefinitions
    % For each x-coordinate -> in total Nx_hs points
    P_x_FE    = zeros(Nx_hs,1);
    P_xKn_FE  = zeros(Nx_hs,NDOF);
    Pz_xKnt_FE = zeros(Nx_hs,NDOF,Nf);
    
    
    % Nodal loads and load direction
    switch load_dir, case 'x', load_dir=2; case 'y', load_dir=1; case 'z',load_dir=0; end
    
    switch calc_cont.mesh_type
        case 'full_cyl'; DOF_skip = 0;
        case 'half_cyl'; DOF_skip = 3*(Nphi_t/2-1);
            % Number of boundary nodes on (not existing)
            % upper half of circle, that has to be skipped
    end
    
    %% Loaded Nodes
    % Nodes on horizontal line (y-direction) through midpoint of tunnel
    % with specific extension in x-direction (=longitudinal direction of tunnel)
    % here only at pos_load_x -> Single load
    
    % Loaded node in x direction
    sx = find(x==pos_load_x_fem);
    
    % node_id of node on midline corr. to pos_load_y_fem
    node_id = nod_midline(nod_midline(:,2) < pos_load_y_fem+0.01 & nod_midline(:,2) > pos_load_y_fem-0.01,1);
    if isempty(node_id); error('Error: pos_load_y_fem not included in nod_midline');end
    
    % ID| y | z of loaded nodes
    loaded_nodes            = Knoten(node_id,:);
    loaded_nodes            = sortrows(loaded_nodes,2);
    
    %% Creation of load matrix
    % load_ampl * load_collecting_width
    % good agreement with ITM point load for dx=dy=dy_FE!!
    P_x_FE(sx) = P0_FE*dy_FE;
    
    % Insert point load in total load matrix
    P_xKn_FE(:,3*node_id-load_dir+DOF_skip) = P_x_FE;
    
    % Structure:
    
    %           Node 1       Node2 ...       Node n
    %           u v w        u v w           u v w
    % -xmax   |                                     |
    %         |                                     |
    % x=0     |                          X          |
    %         |                                     |
    % +xmax   |                                     |
    %                             Pk = pos_load_y_fem
    
    %% Set constant f spectrum for load in 3rd dim
    % Predefine
    P_f = ones(Nx_hs,NDOF,Nf);
    
    % Load const. over f for TF
    P_xKnf_FE = P_xKn_FE.*P_f;
    
    % Fouriertransformation from x to kx
    % -> for each column (-xmax...xmax) over the row dimension
    P_kxKnf_FE = dx/Bx*fftshift(fft(fftshift(P_xKnf_FE,1),[],1),1);
    
    P_FE = P_kxKnf_FE;
    
else
    % Zero load vector
    % For full_cyl NDOF=3*num_nod
    % For half_cyl NDOF different -> same num_nod_bd but other num_nod_in
    P_FE = zeros(Nx_hs,NDOF,Nf);
end

%%  Output

loading.P_FE = P_FE;

%%  Plots

if strcmpi(plot_cont.load,'yes') && P0_FE > 0
    
    %% Plot P(x,y), P(kx,ky), P(f)
    %--------------------------------------------
    fig                     = figure;
    fig.Name                = 'Hs cyl FEM - POINT load';
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------
    % P_{FE}(y)
    subplot(2,2,1) 
    hold on
    plot(nod_midline(:,2),zeros(num_nod_midline,1),'.','Color',plot_cont.color_cyan)
    stem(pos_load_y_fem,P_x_FE(P_x_FE~=0),'Color','r','Marker','.')
    xlabel('$y_{\mathrm{tun}}$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{FE}(y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylim([0 1.2])

    
	% P_{FE}(x)
    subplot(2,2,2)
    hold on
    stem(x,zeros(Nx_hs),'Color',plot_cont.color_cyan,'Marker','.')
    stem(x(sx),P_x_FE(sx),'Color','r','Marker','.')
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{FE}(x)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylim([0 P0_FE+0.2*P0_FE])
    
    % P_{FE}(x,y)
    subplot(2,2,3)
    hold on
    
    [X,Y] = ndgrid(x,nod_midline(:,2));
    stem3(X,Y,zeros(Nx_hs,num_nod_midline),'Marker','.','Color',plot_cont.color_cyan)
    stem3(pos_load_x_fem,pos_load_y_fem,P_x_FE(P_x_FE~=0),'Marker','.','Color','r')
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$y_{\mathrm{tun}}$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    zlabel('$P_{FE}(x,y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    view(3)
    set(gca,'YDir','reverse')

    %P_{FE}(f)
    subplot(2,2,4)
    stem(f,squeeze(P_f(1,1,:)),'color',plot_cont.color_cyan,'Marker','.');
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{FE}(f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    xlim([min(f)*1.2 max(f)+0.2])
    ylim([0 1.2])
   

  
end

