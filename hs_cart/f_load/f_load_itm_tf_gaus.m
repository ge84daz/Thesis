%__________________________________________________________________________
%
%       f_load_itm_tf_gaus.m    
%
%       INPUT: 
%       path      = pathes
%       calc      = calculation setup
%       calc_cont = calculation control parameters
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
%       - Create load matrix in (x,y,z=0,omega) domain with gaussian pulse with midpoint 
%         at a given position (pos_x_hs1L,pos_y_hs1L)
%       - FT into (kx,ky,0,omega) -> Pz_kxky_f 
%       - Unit constant frequency spectrum for all freq specified for calc of TF
%       - Spatial load distribution used for loading in x,y and z direction
%
%       REMARK:
%       Adaptions for calc.flex  = 'yes';
%       -> Factors dx/Bx*dy/By*dt/T for FFT omitted
%       -> P0_itm_hs_1L/bx_itm_hs_1L/by_itm_hs_1L/4 
%       => normalized amplitude -> resultant Ptot=1
%
%       REMARK: 
%       Original Author: Julian Freisinger  (f_load_hs1L_itm_tf_rect)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    12-06-2024 - Hicks   
%__________________________________________________________________________



function [loading] = f_load_itm_tf_gaus(path,calc,plot_cont,loading,dis_itm,geo)

% Initialize: Input

P0_hs1L     = loading.P0_itm_hs1L;
bx_hs1L     = loading.bx_itm_hs1L;
by_hs1L     = loading.by_itm_hs1L;
pos_x_hs1L  = loading.pos_x_itm_hs1L;
pos_y_hs1L  = loading.pos_y_itm_hs1L;

dx          = dis_itm.dx;
dy          = dis_itm.dy;
Nx          = dis_itm.Nx_hs;
Ny          = dis_itm.Ny_hs;  
Nf          = dis_itm.Nf;
f           = dis_itm.f;

kx          = dis_itm.kx;
ky          = dis_itm.ky;

x           = geo.x;
y           = geo.y;
Bx          = geo.Bx;
By          = geo.By;

%% Set up load

if P0_hs1L > 0
    
      % Load distribution
    [~,~,xe] = gauspuls(x,10,0.02);
    [~,~,ye] = gauspuls(y,10,0.02);
        
    % Spatial load vector
    Pz_xy = xe.'*ye;

    %% Time dependence -> 3d matrix
    
    % Constant frequency spectrum 
    % -> unit load for each frequency = constant freq spectrum
    P_f = ones(Nx,Ny,Nf);
        
    % Spatial FFT (x,y) -> (kx,ky)
    Pz_kxky = dx/Bx*dy/By*fftshift(fftn(ifftshift(Pz_xy)));
    % Constant frequency spectrum
    Pz_kxky_f = Pz_kxky.*P_f;
    
else
    Pz_kxky_f  = zeros(Nx_hs,Ny_hs,Nf);
    Pz_xy      = zeros(Nx_hs,Ny_hs);
    Pz_kxky    = zeros(Nx_hs,Ny_hs);
end

%% Output

loading.Pz_kxky_f_hs1L  = Pz_kxky_f;
loading.Pz_xy_hs1L      = Pz_xy;
loading.Pz_kxky_hs1L    = Pz_kxky;

% %% Plots

if strcmpi(plot_cont.load,'yes') && P0_hs1L > 0

        %--------------------------------------------
    fig                     = figure;
    fig.Name                = 'HsSph - Load ITM';
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------
    sgtitle('Spatial distribution of load $P(x,y,z=0)$ and $P(k_x,k_y,z=0)$ for $\mathrm{TF}_{u_{\mathrm{soil}}}$','interpreter','latex','FontSize',12)
    %
    % P(x,y,z=0,omega)
    subplot(2,3,1)
    [X,Y] = ndgrid(x(1:4:end),y(1:4:end));
    mesh(X,Y,Pz_xy(1:4:end,1:4:end),'FaceColor',plot_cont.FaceColor);
    grid on
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    zlabel('$P_{z}(x,y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        
    % P(kx,ky,z=0,omega)
    subplot(2,3,2)
    [KX,KY] = ndgrid(kx(1:4:end),ky(1:4:end));
    mesh(KX,KY,real(Pz_kxky(1:4:end,1:4:end)),'FaceColor',plot_cont.FaceColor)
    grid on
    xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    zlabel('$\mathrm{Re} \; P_{z}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        
    % P(f)
    subplot(2,3,3)
    Pplot = squeeze(P_f(Nx/2+1,Ny/2+1,:));
    plot(f,Pplot)
    grid on
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{z}(f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    xlim([min(f)-1 max(f)+1])
        
    % P(x,y,z=0,omega)
    subplot(2,3,4)
    plot(y,Pz_xy(Nx/2+1,:));
    grid on
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{z}(x=0,y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        
    % P(kx=0,ky,z=0,omega)
    subplot(2,3,5)
    Pplot = squeeze(Pz_kxky_f(Nx/2+1,:,1));
    plot(ky,real(Pplot))
    grid on
    xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{z}(k_x=0,k_y,f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    clear Pplot
    
    % P(kx=0,ky=0,z=0,f)
    subplot(2,3,6)
    Pplot = squeeze(Pz_kxky_f(Nx/2+1,Ny/2+1,:));
    stem(f,Pplot,'.')
    grid on
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{z}(k_x=k_y=0,f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    xlim([min(f)-1 max(f)+1])
    clear Pplot



end

end