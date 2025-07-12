%__________________________________________________________________________
%
%       f_load_hsycyl_itm_tf_line.m    
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
%       INPUT arguments:
%
%       plot_load       = Parameter for Plots yes or no
%       P0_FE           = Load Amplitude
%       num_nod         = Number of FE Nodes
%       Knoten          = Nodes (ID, Ny(i), Nz(i))
%       pos_load_x_fem  = Position of load midpoint in x
%       pos_load_y_fem  = Position of load midpoint in y
%       dx              = Inkrements in x
%       Nx_hs           = Number of Points on halfspace surface
%       Bx              = Width of ground surface in x direction
%       NDOF            = Number of DOFs
%       radius          = Tunnel radius
%       div             = Number of elements per half of radius
%       x               = x discretization
%       kx,ky           = Inkrements in transformed domain
%
%       OUTPUT: 
%       - P_FE = FEM nodal load matrix in transformed domain for time harmonic Line loading.
%
%       DESCRIPTION: 
%       The function determines the FEM nodal load vector for a time harmonic 
%       LINE load within the FEM subsystem. A LINE Load in (x,y,z), which is 
%       constant over all x can be applied to an arbitrary y coordinate at z=0. 
%       The load is distributed on the nodes at y=pos_load_y_fem an one left 
%       and right of this node. The loadvector is transformed into (kx,omega) 
%       domain. The resultant load matrix contains in every row the nodal loads 
%       for each DOF per node. In each column the wavnumbers kx.
%       Over the rows variation from -kmax... +kmax.
%       Over the columns variation from Node 1 (u,v,w,)... Node n (u,v,w)
%       The 3rd dimension varys form -omega_max ... +omega_max
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_load_hscyl_itm_tf_line)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________


function [loading] = f_load_hscyl_fem_tf_line(path,plot_cont,loading,dis_itm,dis_fem,geo,calc_cont,calc)

%% Initialize: Input

P0_FE               = loading.P0_FE;
pos_load_y_fem      = loading.pos_load_y_fem;
pos_load_z_fem      = loading.pos_load_z_fem;
by_fem              = loading.by_fem;
bz_fem              = loading.bz_fem;
load_dir            = loading.load_dir;

Knoten              = dis_fem.Knoten;
NDOF                = dis_fem.NDOF;
div                 = dis_fem.div;
num_nod_midline     = dis_fem.num_nod_midline;

nod_midline         = dis_fem.nod_midline;
dy_FE               = dis_fem.dy_FE;

dx                  = dis_itm.dx;
Nx_hs               = dis_itm.Nx_hs;
Nphi_t              = dis_itm.Nphi_t;
kx                  = dis_itm.kx;
f                   = dis_itm.f;
Nf                  = dis_itm.Nf;

Bx                  = geo.Bx;
radius              = geo.radius;
x                   = geo.x;


% Creation of load matrix

if P0_FE > 0
    
    %%  Predefinition
    % For each x-coordinate -> in total Nx_hs points
    Pz_x_FE    = zeros(Nx_hs,1);
    Pz_xKn_FE  = zeros(Nx_hs,NDOF);
    Pz_xKnf_FE = zeros(Nx_hs,NDOF,Nf);
    
    % if dy_FE > by_fem; error('Too less FEM nodes loaded -> increase width of FEM Line load'); end

    %% Nodal loads and load direction
    
    % Loaded nodes in x direction
    sx = 1:Nx_hs;   % Constant load over all points in x direction
    
    % Calculation of equivalent nodal loads for FEM that are equal to lineload of ITM system.
    % Nodal load = lineload*load collecting width (dy_FE=radius/(2*div))
    Pz_x_FE(sx) = P0_FE*dy_FE;
    
    % Direction of the load 
    switch load_dir; case 'x', load_dir=2; case 'y', load_dir=1; case 'z',load_dir=0; end
    
    switch calc_cont.mesh_type
        case 'full_cyl';  DOF_skip = 0;
        case 'half_cyl';   DOF_skip = 3*(Nphi_t/2-1);
            % Number of boundary nodes on (not existing)
            % upper half of circle, that has to be skipped   
    end
    
    %% Node IDs of loaded nodes
    % loaded nodes on horizontal line in y-dir through midpoint of tunnel
    % over the total tunnel length (x-dir)
    
    nby = floor(by_fem/(radius/2/div));
    
    switch pos_load_z_fem
        case 'z0'
            node_id = nod_midline((num_nod_midline-1)/2+1+pos_load_y_fem/(radius/2/div)-nby:(num_nod_midline-1)/2+1+pos_load_y_fem/(radius/2/div)+nby);
    end
    
    loaded_nodes           = Knoten(node_id,:);
    loaded_nodes           = sortrows(loaded_nodes,2);
    num_loaded_nodes       = size(loaded_nodes,1);
    loaded_nodes_fem_line  = loaded_nodes;
    node_id_sort           = loaded_nodes(:,1);
    
    for i=1:num_loaded_nodes
        Pz_xKn_FE(:,(3*node_id_sort(i)-load_dir+DOF_skip)) = Pz_x_FE;
    end
    
    % The outer nodes w.r.t. y-direction have only half of the contribution area dy_FE
    Pz_xKn_FE(:,(3*node_id_sort(1)-load_dir+DOF_skip))                = Pz_x_FE*0.5;
    Pz_xKn_FE(:,(3*node_id_sort(num_loaded_nodes)-load_dir+DOF_skip)) = Pz_x_FE*0.5;
    
    % Load Amplitudes
    load_ampl_fem_line    = Pz_xKn_FE(Nx_hs/2+1,(3*node_id_sort.'-load_dir+DOF_skip));
    load_ampl_fem_line_xy = ones(Nx_hs,1)*load_ampl_fem_line;
    
    
    %% Creation of load matrix
    
    % Predefine constant f spectrum in 3rd dim
    P_f = ones(Nx_hs,NDOF,Nf);
    
    % Load const. over f for TF
    Pz_xKnf_FE = Pz_xKn_FE.*P_f;
    
    % Fourier Transform from x to kx -> over 1st dim (rows) 
    % for all x from -xmax ... + xmax to -kmax ... kmax
    Pz_kxKnf_FE = dx/Bx*fftshift(fft(fftshift(Pz_xKnf_FE,1),[],1),1);   % over cloumn   x -> kx
    
    P_FE = Pz_kxKnf_FE;
    
else 
    P_FE = zeros(Nx_hs,NDOF,Nf);
end

%%  Output

loading.P_FE = P_FE;

if P0_FE > 0
    
    loading.loaded_nodes_fem_line    = loaded_nodes_fem_line;
    loading.load_ampl_fem_line       = load_ampl_fem_line;
    loading.load_ampl_fem_line_xy    = load_ampl_fem_line_xy;
end
%% Plots

if strcmpi(plot_cont.load,'yes') && P0_FE > 0
    
    switch calc.structure
        case {'tunnel'}
            
            %--------------------------------------------
            fig                     = figure;
            fig.Name                = 'FEM harmonic LINE Load';
            fig.Position            = plot_cont.fig_position;
            fig.PaperPositionMode   = 'auto';
            %--------------------------------------------
            
            % P_{FE}(y)
            subplot(2,2,1)
            hold on
            plot(nod_midline(:,2)',zeros(size(nod_midline,1)),'Color',plot_cont.color_cyan,'Marker','.')
            stem(loaded_nodes_fem_line(:,2)',load_ampl_fem_line,'Color','r','Marker','.')
            set(gca,'XDir','reverse');
            xlabel('$y$','Interpreter',plot_cont.interpreter)
            ylabel('$P_{\mathrm{FE}}(y)$','Interpreter',plot_cont.interpreter)
          
            % P_{FE}(x)
            subplot(2,2,2)
            hold on
            stem(x,zeros(Nx_hs),'Color',plot_cont.color_cyan,'Marker','.')
            stem(x,Pz_xKn_FE(:,3*node_id(2)+DOF_skip),'Color','r','Marker','.')
            set(gca,'XDir','reverse');
            xlabel('$x$','Interpreter',plot_cont.interpreter)
            ylabel('$P_{\mathrm{FE}}(x)$','Interpreter',plot_cont.interpreter)

            
            % P_{FE}(x,y)
            subplot(2,2,3)
            hold on
          
            % loaded_nodes_xy = ones(Nx_hs,1)*loaded_nodes_xy';
            [X,Y] = ndgrid(x,nod_midline(:,2));
            plot3(X,Y,zeros(size(X)),'.','Color',plot_cont.color_cyan)
            [X,Y_l] = ndgrid(x,loaded_nodes_fem_line(:,2));
            stem3(X,Y_l,load_ampl_fem_line_xy,'Color','r','Marker','.')
            set(gca,'XDir','reverse');
            xlabel('$x$','Interpreter',plot_cont.interpreter)
            ylabel('$y_{tunnel}$','Interpreter',plot_cont.interpreter)
            zlabel('$P_{\mathrm{FE}}(x,y)$','Interpreter',plot_cont.interpreter)
            view(3)
            
            %P_{FE}(f)
            subplot(2,2,4)
            stem(f,squeeze(P_f(1,1,:)),'Marker','.')   
            xlabel('$f$','Interpreter',plot_cont.interpreter)
            ylabel('$P_{\mathrm{FE}}(f)$','Interpreter',plot_cont.interpreter)
            xlim([min(f)*1.2 max(f)+0.2])
            ylim([0 1.2])    
            
    end
end
end
