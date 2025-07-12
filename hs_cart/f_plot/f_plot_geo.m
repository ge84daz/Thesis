%__________________________________________________________________________
%
%       f_plot_geo.m    
%
%       INPUT: 
%       path      = pathes
%       calc      = calculation setup
%       plot_cont = plot control parameters
%       loading   = load amplitude, width, position
%       dis_itm   = discretization 
%       geo       = geometry
%
%       OUTPUT: 
%       Pz_kxky_f    = load matrix itm at z=0 in (kx,ky,omega)->(Nx_hs,Ny_hs,Nf);
%       Pz_xy        = load matrix itm at z=0 in (x,y)        ->(Nx_hs,Ny_hs);
%       Pz_kxky      = load matrix itm at z=0 in (kx,ky)      ->(Nx_hs,Ny_hs);
%
%       DESCRIPTION: 
%       - Plot geometry of hs, hs1L
%       - Plot loading distribution
%
%       REMARK: 
%       Original Author: Julian Freisinger  (f_plot_geo_hs1L)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________


function [] = f_plot_geo(path,~,plot_cont,loading,geo)

%% Initialize: Input

x      = geo.x;
y      = geo.y;

Bx     = geo.Bx;
By     = geo.By;

Pz_xy_hs1L = loading.Pz_xy_hs1L;

%% Geometry 

z_geo_z0  = zeros(size(y));

%% Loading
% Location of load for (x,y,z=0)
[x_P_hs1L,y_P_hs1L,~] = find(Pz_xy_hs1L);
z_P_hs1L_xy           = 0.1*ones(size(x_P_hs1L));

% Location and amplitude of load for (x=0,y,z=0)
y_hs1L_z0_x0 = find(Pz_xy_hs1L(x==0,:));

% set z-coord. for load plot (x=0,y,z)
if strcmpi(loading.pos_z_itm_hs1L,'z0')
    z_P_hs1L_x0y = 0.1*ones(size(y_hs1L_z0_x0)); 
end
%% Plots

%-----------------------------------------------------
mainfigurename = 'Hs1L - Geometry and Load distribution';
fig = figure('Name',mainfigurename);
fig.Position = [632,260,1403,1040];
%-----------------------------------------------------
%
subplot(1,2,1)
hold on
% Geometry

plot(y,z_geo_z0,'.','color',plot_cont.color_cyan);
  
% Load
plot(y(y_hs1L_z0_x0),z_P_hs1L_x0y,'.','color','red')

xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel('$z$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
axis([-By/2 By/2 -By/2 By/2])
pbaspect([1 1 1]);


subplot(1,2,2)
hold on
% Geometry
[X,Y] = meshgrid(x,y);
plot3(X,Y,zeros(size(X)),'.','color',plot_cont.color_cyan)
% Load
plot3(x(x_P_hs1L),y(y_P_hs1L),z_P_hs1L_xy,'.','color','red')
xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
xlim([-Bx/2-Bx/16 +Bx/2+Bx/16])
ylim([-By/2-By/16 +By/2+By/16])
pbaspect([1 1 1]);
view([90 90])

sgtitle(['Load type:  ' loading.type],'Interpreter',plot_cont.interpreter)

savename = '\geo_hs1L';
saveas(fig,[path.figures savename],'fig')

if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig);end



%% Output

