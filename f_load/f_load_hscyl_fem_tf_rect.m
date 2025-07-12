%__________________________________________________________________________
%
%       f_load_hsycyl_fem_tf_rect.m    
%
%       INPUT: 
%       - path          = pathes
%       - plot_cont     = plot control parameter
%       - loading       = load parameters
%       - dis_itm       = discretization ITM
%       - dis_fem       = discretization FEM
%       - geo           = geometry
%       - calc_cont     = calculation control prameter
%       - calc          = calculation param
%
%       OUTPUT: 
%       - P_FE                              = FEM nodal load matrix in transformed domain for time harmonic Point laoding.
%       - P0_tot_nodes                      = Sum of nodal loads
%       - loaded_nodes_fem_block_x          = loaded nodes in x
%       - loaded_nodes_fem_block_y          = loaded nodes in y
%       - load_ampl_fem_line_by             = load amplitudes over y
%       - load_ampl_xy_fem_block            = load amplitudes of block
%
%       DESCRIPTION: 
%       The function determines the load vector P_FE for a harmonic rectangular
%       block load within the FEM domain. The resultant(=sum of the nodal forces)
%       of the load is equal to the volume of a rectangular distributed load.
%       Therefore nodal load distribution:
%
%         0.25  0.50  0.25
%
%         0.50  1.00  0.50
%
%         0.25  0.50  0.25
%
%       Structure Pz_xKn_FE:
%
%                Node 1       Node2 ...       Node n
%                   u v w        u v w           u v w
%       -xmax   |                                     |
%               |                         x x x       |
%       x=0     |                         x X x       | nodeID = pos_load_x_fem
%               |                         x x x       |
%       +xmax   |                                     |
%                                     nodeID = pos_load_y_fem
%       
%       REMARK: 
%       Original Author: Julian Freisinger (f_load_hscyl_fem_tf_rect)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________


function [loading] = f_load_hscyl_fem_tf_rect(path,plot_cont,loading,dis_itm,dis_fem,geo,calc_cont,calc)

%% Initialize: Input

dx       = dis_itm.dx;
dy       = dis_itm.dy;
kx       = dis_itm.kx;
Nx_hs    = dis_itm.Nx_hs;
Nphi_t   = dis_itm.Nphi_t;
Nf       = dis_itm.Nf;
f        = dis_itm.f;

Bx       = geo.Bx;
radius   = geo.radius;
x        = geo.x;
dy_FE    = dis_fem.dy_FE;                   % Element length in y-direction

P0_FE            = loading.P0_FE     ;      % Amplitude
pos_load_x_fem   = loading.pos_load_x_fem;  % Load position
pos_load_y_fem   = loading.pos_load_y_fem;
load_dir         = loading.load_dir;         % Load direction
bx_fem           = loading.bx_fem;           % Load width
by_fem           = loading.by_fem;

Knoten           = dis_fem.Knoten;           % FEM nodes + coord
NDOF             = dis_fem.NDOF;             % Number of DOFS
div              = dis_fem.div;              % Number of elements per half of radius
num_nod_midline  = dis_fem.num_nod_midline;  % Number of nodes at x=0
nod_midline      = dis_fem.nod_midline;      % Nodes at x=0 over y
% in cyl middle resp. top of concrete slab

%% Creation of load matrix

if P0_FE > 0
    
    %% Predefinition
    % Load vector and matrix
    Pz_x_FE    = zeros(Nx_hs,1);
    Pz_xKn_FE  = zeros(Nx_hs,NDOF);
    Pz_xKnt_FE = zeros(Nx_hs,NDOF,Nf);
    
    % Nodal loads and load direction
    % Direction of the load
    switch load_dir, case 'x', load_dir=2; case 'y', load_dir=1; case 'z',load_dir=0; end
    
    switch calc_cont.mesh_type
        case 'full_cyl'; DOF_skip=0;
        case'half_cyl';  DOF_skip=3*(Nphi_t/2-1);
            % Number of boundary nodes on (not existing)
            % upper half of circle, that has to be skipped
    end
    
    %% Loaded Nodes
    
    % Loaded node in x direction
    nbx = bx_fem/dx; % Number of dx of load block in x
    sx  = Nx_hs/2+1+(pos_load_x_fem/dx)-nbx:Nx_hs/2+1+(pos_load_x_fem/dx)+nbx; % symmetric to origin (pos_load_x,pos_load_y)
    if size(sx,2)>size(x,2), sx(size(sx,2))=[]; end
    
    % Loaded node in y direction -> node_ids of loaded nodes

    % Loadmidpoint
    loadmidpoint(1) = pos_load_y_fem;
    loadmidpoint(2) = 0;                % z=0
            
    [pv] = [loadmidpoint(1)-by_fem-dy_FE/2, loadmidpoint(2)+ 0.01; loadmidpoint(1)+by_fem+dy_FE/2, loadmidpoint(2)+ 0.01;...
             loadmidpoint(1)+by_fem+dy_FE/2, loadmidpoint(2)- 0.01; loadmidpoint(1)-by_fem-dy_FE/2, loadmidpoint(2)- 0.01];
            
    % Find all nodes within polygon around loadmindpoint
    [in]    = inpolygon(nod_midline(:,2),nod_midline(:,3),pv(:,1),pv(:,2));
    node_id = nod_midline(in);
    
    loaded_nodes     = Knoten(node_id,:);
    loaded_nodes     = sortrows(loaded_nodes,2);
    num_loaded_nodes = size(loaded_nodes,1);
    node_id_sort     = loaded_nodes(:,1);
    
     %% Creation of load matrix
    
    % Equivalent nodal loads
    % load_ampl * load_collecting_width
    P0_FE = P0_FE*dy_FE;
    
    % Load distribution along x
    Pz_x_FE(sx)      = P0_FE;
    if ~strcmpi(calc.structure,{'rw','rw_gab','rw_sub'}) % then load with full amplitude in x and y direction at top of the rails
    Pz_x_FE(sx(1))   = 0.5*P0_FE;
    Pz_x_FE(sx(end)) = 0.5*P0_FE;
    end
    % Insert Point Load in total load matrix
    for i=1:num_loaded_nodes
        Pz_xKn_FE(:,(3*node_id_sort(i)-load_dir+DOF_skip))=Pz_x_FE;
    end
    if ~strcmpi(calc.structure,{'rw','rw_gab','rw_sub'}) % then load with full amplitude in x and y direction at top of the rails
    Pz_xKn_FE(:,(3*node_id_sort(1)-load_dir+DOF_skip))                = Pz_x_FE*0.5;
    Pz_xKn_FE(:,(3*node_id_sort(num_loaded_nodes)-load_dir+DOF_skip)) = Pz_x_FE*0.5;
    end
    
    % Sum of nodal loads
    P0_tot_nodes = sum(sum(Pz_xKn_FE));
    
    % Amplitudes of loaded nodes over y
    load_ampl_fem_block_y = Pz_xKn_FE(Nx_hs/2+1,(3*node_id_sort'-load_dir+DOF_skip));
    
    % Amplitudes of loaded nodes over x and y
    load_ampl_fem_block_xy = zeros(Nx_hs,num_loaded_nodes); % Predefine
    
    % Get load amplitudes
    for i=1:num_loaded_nodes
        load_ampl_fem_block_xy(:,i) = Pz_xKn_FE(:,(3*node_id_sort(i)-load_dir+DOF_skip));
    end
    
    % Delete all zero rows
    load_ampl_fem_block_xy( all(~load_ampl_fem_block_xy,2), : ) = [];
    
    % Predefine constant f spectrum in 3rd dim
    P_f = ones(Nx_hs,NDOF,Nf);
    
    % Load const. over f for TF
    P_xKnf_FE = Pz_xKn_FE.*P_f;
    
    % Fouriertransformation
    P_kxKnf_FE     = dx/Bx*fftshift(fft(fftshift(P_xKnf_FE,1),[],1),1); % FFT x -> kx
    
    P_FE = P_kxKnf_FE;
    
else
    % For full_cyl NDOF=3*num_nod
    % For half_cyl NDOF different -> same num_nod_bd but other num_nod_in
    P_FE = zeros(Nx_hs,NDOF,Nf);
    P0_tot_nodes = 0;
    
end

%%  Output

loading.P_FE         = P_FE;
loading.P0_tot_nodes = P0_tot_nodes;

if P0_FE > 0
    loading.loaded_nodes_fem_block_x    = x(sx);
    loading.loaded_nodes_fem_block_y    = loaded_nodes(:,2);
    loading.loaded_nodes_fem_block_z    = loaded_nodes(:,3);
    loading.load_ampl_fem_block_y       = load_ampl_fem_block_y;
    loading.load_ampl_fem_block_xy      = load_ampl_fem_block_xy;
end

%%  Plots

if strcmpi(plot_cont.load,'yes') && P0_FE > 0
    %--------------------------------------------
    fig                     = figure;
    fig.Name                = 'HsCyl - Load FEM RECT';
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------
    
    % P_{FE}(y)
    subplot(2,2,1)
    hold on
    plot(nod_midline(:,2)',zeros(size(nod_midline,1)),'color',plot_cont.color_cyan)
    stem(loaded_nodes(:,2)',load_ampl_fem_block_y,'color','r','Marker','.')
    set(gca,'XDir','reverse');
    xlabel('$y$','Interpreter',plot_cont.interpreter)
    ylabel('$P_{FE}(y)$','Interpreter',plot_cont.interpreter)
    
    % P_{FE}(x)
    subplot(2,2,2)
    hold on
    Pplot = Pz_xKn_FE(:,3*node_id(2)+DOF_skip);
    stem(x,zeros(Nx_hs),'Color',plot_cont.color_cyan,'Marker','.')
    stem(x(Pplot~=0),Pplot(Pplot~=0),'color','r','Marker','.')
    set(gca,'XDir','reverse');
    xlabel('$x$','Interpreter',plot_cont.interpreter)
    ylabel('$P_{FE}(x)$','Interpreter',plot_cont.interpreter)
    
    % P_{FE}(x,y)
    subplot(2,2,3)
    hold on
    [X,Y] = ndgrid(x,nod_midline(:,2));
    [X_l,Y_l] = ndgrid(x(sx),loaded_nodes(:,2));
    plot3(X,Y,zeros(size(X)),'.','Color',plot_cont.color_cyan)
    stem3(X_l,Y_l,load_ampl_fem_block_xy,'Color','r','Marker','.')
    set(gca,'XDir','reverse');
    xlabel('$x$','Interpreter',plot_cont.interpreter)
    ylabel('$y_{tunnel}$','Interpreter',plot_cont.interpreter)
    view([-15.236170212766 22.7526760563381]);

    %P_{FE}(f)
    subplot(2,2,4)
    stem(f,squeeze(P_f(1,1,:)),'Marker','.')
    xlabel('$f$','Interpreter',plot_cont.interpreter)
    ylabel('$P_{\mathrm{FE}}(f)$','Interpreter',plot_cont.interpreter)
    xlim([min(f)*1.2 max(f)+0.2])
    ylim([0 1.2])
            
   
end
end

