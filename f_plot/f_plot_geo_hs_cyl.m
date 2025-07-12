%__________________________________________________________________________
%
%       f_plot_geo_hs_cyl.m    
%
%       INPUT: 
%       - geo_T1     = Geometry of tunnel/trench 1
%       - dis_itm    = discretization ITM subsystem 
%       - dis_fem_T1 = discretization of FEM subsystem tunnel/trench 1
%       - loading    = load parameter
%       - path       = path variables
%       - plot_cont  = plot control parameters
%       - calc       = calculation parameters
%
%       OUTPUT: 
%       - plot of geometry and load of hs_cyl_1L and reference system hs1L
%
%       DESCRIPTION: 
%       -> y = [-By/2 ... +By/2-dy]. However, in plots XDir plotted reverse with positive axis pointing to the left.
%       -> not all points plotted to fasten plot and save memory
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_plot_geo_hs_cyl)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [] = f_plot_geo_hs_cyl(geo_T1,dis_fem_T1,loading,path,plot_cont,calc)

%% Initialize: Input
% Size of '.' for load points over (x,y)
MarkerSizeLoadxy = 8;

% geometry
x         = geo_T1.x;
y         = geo_T1.y;
By        = geo_T1.By;
Bx        = geo_T1.Bx;

% Tunnel T1
H_T1         = geo_T1.H;
y_Tc_T1      = geo_T1.y_Tc;
radius_T1    = geo_T1.radius;
H_cyl_tot_T1 = geo_T1.H_cyl_tot;


% dis_fem_T1
Edof_T1         = dis_fem_T1.Edof;
Knoten_T1       = dis_fem_T1.Knoten;
elem_coord_T1   = dis_fem_T1.elem_coord;
nod_midline_T1  = dis_fem_T1.nod_midline;

% Loading ITM
P0_ITM_hscyl      = loading.P0_itm_hs;
P_ITM_hscyl_xy    = loading.Pz_ITM_xy;

% Loading FEM
P0_FE          = loading.P0_FE;
P_FE_T1        = loading.P_FE_T1; 
pos_load_x_fem = loading.pos_load_x_fem;
pos_load_y_fem = loading.pos_load_y_fem;
pos_load_z_fem = loading.pos_load_z_fem;
load_dir       = loading.load_dir;

if P0_FE > 0 && strcmpi(loading.distr_xy,'line')
    loaded_nodes_fem_line   = loading.loaded_nodes_fem_line;
    P0_FEM_hscyl_line_x0y   = loading.load_ampl_fem_line;
    P0_FEM_hscyl_line_xy    = loading.load_ampl_fem_line_xy;
end

if P0_FE > 0 && strcmpi(loading.distr_xy,'rect')
    loaded_nodes_fem_block_x = loading.loaded_nodes_fem_block_x;
    loaded_nodes_fem_block_y = loading.loaded_nodes_fem_block_y;
    loaded_nodes_fem_block_z = loading.loaded_nodes_fem_block_z;
    load_ampl_fem_block_y    = loading.load_ampl_fem_block_y;
    load_ampl_fem_block_xy   = loading.load_ampl_fem_block_xy;
end


%% Hscyl - Coordinates of ITM and FEM part
switch calc.system
    case {'trench'}
        % nodes z=0 w/o nodes inside trench
        y_geo_hscyl_z0_T1 = [y(y<(y_Tc_T1-radius_T1)),y(y>y_Tc_T1+radius_T1)];
        z_geo_hscyl_z0_T1 = [zeros(1,size(y(y<-(radius_T1+y_Tc_T1)),2)), zeros(1,size(y(y>radius_T1-y_Tc_T1),2))];
        % Tunnel 1 nodes 
        y_geo_hscyl_T1    = Knoten_T1(:,2).'+ y_Tc_T1;  
        z_geo_hscyl_T1    = Knoten_T1(:,3).';           
        
    case 'tunnel'    
        % nodes z=0
        y_geo_hscyl_z0_T1  = y;                             
        z_geo_hscyl_z0_T1  = [zeros(1,size(y,2))];            
        
        % Tunnel 1 nodes
        y_geo_hscyl_T1  = Knoten_T1(:,2).' + y_Tc_T1;   
        z_geo_hscyl_T1  = Knoten_T1(:,3).'+H_T1;     
end

% All points on (x=0,y,z=0) without tunnel
switch calc.system
    case {'trench'}  
        y_geo_hscyl_z0 = [y(y<=y_Tc_T1-radius_T1),nod_midline_T1(:,2).'+y_Tc_T1,y(y>=y_Tc_T1+radius_T1)];
    case {'tunnel'}
        y_geo_hscyl_z0 = y;
end

% All points on (x,y,z=0) without tunnel
% only each 4th point -> faster
[X_hscyl_OF_WO_TUN,Y_hscyl_OF_WO_TUN] = ndgrid(x(1:4:end),y_geo_hscyl_z0(1:4:end));


%% Hscyl - Load ITM
% Direction of the load x (2), y (1) or z (0)
switch load_dir, case 'x', load_dir=2; case 'y', load_dir=1; case 'z',load_dir=0; end

% Loaded nodes P_ITM_hscyl(x=0,y,z=0)
y_ITM_hscyl_x0yz0  = y(P_ITM_hscyl_xy(x==0,:)~=0);
P_ITM_hscyl_x0yz0  = -0.1*P_ITM_hscyl_xy(x==0,P_ITM_hscyl_xy(x==0,:)~=0);

% Loaded nodes P_ITM(x,y,z=0)
[X,Y] = ndgrid(x,y);
X_ITM_hscyl_xyz0 = X(P_ITM_hscyl_xy~=0);
Y_ITM_hscyl_xyz0 = Y(P_ITM_hscyl_xy~=0);

P_ITM_hscyl_xyz0 = -0.1*P_ITM_hscyl_xy(P_ITM_hscyl_xy~=0);
clear row col

%% Hscyl - Load FEM
% For all FEM load cases
switch calc.system
    case {'trench'}
        if strcmpi(calc.structure,'trench');        H_pl_hscyl_FEM = -0.1;end
    case 'tunnel'
        if strcmpi(calc.structure,'tunnel');  H_pl_hscyl_FEM = H_T1-0.1;end
end

switch load_dir
    case 2; MarkerFEM = '.';%'x'; %x
    case 1; MarkerFEM = '.';%'<'; %y
    case 0; MarkerFEM = '.';%'v'; %z
end
%% Hscyl - Load FEM POINT    
% Loaded nodes FEM POINT of coupled System in y-direction and over xy
if P0_FE > 0 && strcmpi(loading.distr_xy,'point')
    
    x_FEM_hscyl_point     = pos_load_x_fem;
    y_Pfem_hscyl_point_T1 = pos_load_y_fem + y_Tc_T1;
    
    
    P_FEM_hscyl_point = H_pl_hscyl_FEM*P0_FE;
end

%% Hscyl - Load FEM LINE
% Loaded nodes FEM LINE of coupled System in y-direction
if P0_FE > 0 && strcmpi(loading.distr_xy,'line')
    
    y_FEM_hscyl_x0y_line_T1 = loaded_nodes_fem_line(:,2).'+y_Tc_T1;
        
    P_FEM_hscyl_x0y_line = H_pl_hscyl_FEM.*ones(1,size(y_FEM_hscyl_x0y_line_T1,2));
    % P0_FEM_hscyl_line_x0y % from input
    
    [X_FEM_hscyl_line,Y_FEM_hscyl_line_T1] = ndgrid(x,y_FEM_hscyl_x0y_line_T1.');
    
    
    % P0_FEM_hscyl_line_xy  % from input
end

%% Hscyl - Load FEM RECT
if P0_FE > 0 && strcmpi(loading.distr_xy,'rect')
   
    x_FEM_hscyl_xy_rect     = loaded_nodes_fem_block_x;
   
    y_FEM_hscyl_x0y_rect_T1 = loaded_nodes_fem_block_y + y_Tc_T1;
       
    [X_FEM_hscyl_xy_rect_T1,Y_FEM_hscyl_xy_rect_T1] = ndgrid(x_FEM_hscyl_xy_rect,y_FEM_hscyl_x0y_rect_T1);
        
    P_FEM_hscyl_x0y_rect = H_pl_hscyl_FEM*ones(size(load_ampl_fem_block_y)); 
    P_FEM_hscyl_xy_rect  = H_pl_hscyl_FEM*ones(size(load_ampl_fem_block_xy));
   
end


%% Plot of Geometry and Load for Report File -> Without Text in Subfigures
if strcmpi(plot_cont.geo,'yes')
    %--------------------------------------------------------------------------
    fig                     = figure;
    fig.Name                = 'Geometry and Load of system halfspace with cylindrical cavity';
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------------------------------------
    
    %% Hscyl - Load ITM (x=0,y)
    %---------------------------------------
    subplot(2,2,1)
    %---------------------------------------
    figure_title = 'Hscyl - Load ITM';
    hold on
    
    % plot nodes
    plot(y_geo_hscyl_z0_T1,z_geo_hscyl_z0_T1,'.','color',plot_cont.color_cyan);
    plot(y_geo_hscyl_T1,z_geo_hscyl_T1,'.','color',plot_cont.color_cyan);
    

    
    % plot load
    if P0_ITM_hscyl~=0 
        plot(y_ITM_hscyl_x0yz0,P_ITM_hscyl_x0yz0,'.','color','red')
    end
    
    grid on;
    axis([-By/2 By/2 -By/2 By/2]); % quadratic
    pbaspect([1 1 1]);
    
    xlabel('$y$','Interpreter',plot_cont.interpreter)
    ylabel('$z$','Interpreter',plot_cont.interpreter)
    title(figure_title,'Interpreter',plot_cont.interpreter);
    set(gca,'XDir','reverse');
    set(gca,'YDir','reverse');
    clear figurename
    
    %% Hscyl - Load FEM (x=0,y)
    %---------------------------------------
    subplot(2,2,2)
    %---------------------------------------
    figure_title = 'Hscyl - Load FEM';
    hold on
    
    % plot nodes
    plot(y_geo_hscyl_z0_T1,z_geo_hscyl_z0_T1,'.','color',plot_cont.color_cyan);
    plot(y_geo_hscyl_T1,z_geo_hscyl_T1,'.','color',plot_cont.color_cyan);
        
    % plot load - not implemented so far
    if P0_FE > 0
        
        switch loading.distr_xy
            
            %-----------
            case 'point'
            %-----------
            
                if strcmpi(P_FE_T1,'yes');  plot(y_Pfem_hscyl_point_T1,P_FEM_hscyl_point,'LineStyle','none','Marker',MarkerFEM,'color','red','MarkerSize',4);end
                             
            %-----------
            case 'line'
            %-----------
            
            switch pos_load_z_fem
                case {'z0'}
                    if strcmpi(P_FE_T1,'yes');plot(y_FEM_hscyl_x0y_line_T1,P_FEM_hscyl_x0y_line,MarkerFEM,'LineStyle','none','Marker',MarkerFEM,'color','red','MarkerSize',4);end
            end

            %-----------
            case 'rect'
            %-----------
                
                if strcmpi(P_FE_T1,'yes');  plot(y_FEM_hscyl_x0y_rect_T1,P_FEM_hscyl_x0y_rect,'LineStyle','none','Marker',MarkerFEM,'color','red','MarkerSize',4);end %'.','color','red'
        end
        
    end
    
    grid on;
    axis([-By/2 By/2 -By/2 By/2]);
    pbaspect([1 1 1]);
    
    xlabel('$y$','Interpreter',plot_cont.interpreter)
    ylabel('$z$','Interpreter',plot_cont.interpreter)
    title(figure_title,'Interpreter',plot_cont.interpreter);
    set(gca,'XDir','reverse');
    set(gca,'YDir','reverse');
    clear figurename
    
   
    
    %% Hscyl - Load ITM (x,y)
    subplot(2,2,3)
    %---------------------------------------
    figure_title = 'Hscyl - Load ITM';
    hold on
    
    % Plot Points
    plot3(X_hscyl_OF_WO_TUN,Y_hscyl_OF_WO_TUN,zeros(size(X_hscyl_OF_WO_TUN)),'.','color',plot_cont.color_cyan);
    
    % Plot load
    if P0_ITM_hscyl ~= 0
        plot3(X_ITM_hscyl_xyz0,Y_ITM_hscyl_xyz0,P_ITM_hscyl_xyz0,'LineStyle', 'none','Marker','.','Color','red','MarkerSize',MarkerSizeLoadxy);
    end
    
    xlabel('$x$','Interpreter',plot_cont.interpreter)
    ylabel('$y$','Interpreter',plot_cont.interpreter)
    title(figure_title,'Interpreter',plot_cont.interpreter);
    xlim([-Bx/2-Bx/16 +Bx/2+Bx/16])
    ylim([-By/2-By/16 +By/2+By/16])
    pbaspect([1 1 1]);
    view(270, 90);
    set(gca,'XDir','reverse');
    set(gca,'ZDir','reverse');
    
    %% Hscyl - Load FEM (x,y)
    subplot(2,2,4)
    %---------------------------------------
    figure_title = 'Hscyl - Load FEM';
    hold on
    
    % Plot Points
    plot3(X_hscyl_OF_WO_TUN,Y_hscyl_OF_WO_TUN,zeros(size(X_hscyl_OF_WO_TUN)),'.','color',plot_cont.color_cyan);
    
    if P0_FE > 0
    % Plot load
    switch loading.distr_xy
        
        %-----------
        case 'point'
        %-----------
        
        if strcmpi(P_FE_T1,'yes');  plot3(x_FEM_hscyl_point,y_Pfem_hscyl_point_T1,P_FEM_hscyl_point,'LineStyle','none','Marker',MarkerFEM,'color','red','MarkerSize',4);end
                
        %-----------
        case 'line'
        %-----------

        switch pos_load_z_fem
            case {'z0'}
                if strcmpi(P_FE_T1,'yes');plot3(X_FEM_hscyl_line,Y_FEM_hscyl_line_T1,-P0_FEM_hscyl_line_xy,'LineStyle','none','Marker',MarkerFEM,'color','red','MarkerSize',4);end
              
        end   
        
       %-----------
        case 'rect'
       %-----------

            if strcmpi(P_FE_T1,'yes'); plot3(X_FEM_hscyl_xy_rect_T1,Y_FEM_hscyl_xy_rect_T1,P_FEM_hscyl_xy_rect,'Color','red','Marker','.');  end
          
    end
    end
    xlabel('$x$','Interpreter',plot_cont.interpreter)
    ylabel('$y$','Interpreter',plot_cont.interpreter)
    title(figure_title,'Interpreter',plot_cont.interpreter);
    xlim([-Bx/2-Bx/16 +Bx/2+Bx/16])
    ylim([-By/2-By/16 +By/2+By/16])
    pbaspect([1 1 1]);
    view(270, 90);
    set(gca,'XDir','reverse');
    set(gca,'ZDir','reverse');
    
   
%% Plot Undeformed FEM mesh

if strcmpi(plot_cont.geo,'yes')
    %--------------------------------------------------------------------------
    fig                     = figure;
    fig.Name                = 'Undeformed FE mesh';
    fig.Position            = plot_cont.fig_position_2;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------------------------------------
    pl_row = 1;
    pl_col = 1;
    %
    subplot(pl_row,pl_col,1)
    ex       = elem_coord_T1(:,[3,6,9,12]);
    ey       = elem_coord_T1(:,[4,7,10,13]);
    plotpar  = [1 1 0];
    Elnum    = Edof_T1(:,1);
    eldraw2(ex,ey,plotpar)%,Elnum) % from calfem-3.4
    
       
    title({'Undeformed FEM mesh T1'},'Color',[0 0 0], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$z$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    set (gca,'YDir','reverse')
    axis square
    xlim([-max(radius_T1)  +max(radius_T1)])
    ylim([-max(radius_T1)  +max(radius_T1)])
    
   
end

end




