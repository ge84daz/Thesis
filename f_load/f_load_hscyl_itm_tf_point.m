%__________________________________________________________________________
%
%       f_load_hsycyl_itm_tf_point.m    
%
%       INPUT: 
%       - path          = pathes
%       - plot_cont     = plot control parameter
%       - loading       = load parameters
%       - dis_itm       = discretization
%       - geo           = geometry
%       - nnPi          = loop variable for load position y for flexibility
%
%       OUTPUT: 
%       - loading       = load parameters
%       - Px_ITM, Py_ITM, Pz_ITM  = Total load vector ITM at z=0 with normalized amplitude -> resultant Ptot=1
%       - Pz_ITM_kxky, Pz_ITM_xy  = load vector ITM z=0 of (kx,ky) resp. (x,y)
%
%       DESCRIPTION: 
%       - Same spatial load distribution used for loading in x,y and z direction
%       - set point load at resp. position in Pz_xy
%       - fft2 -> no unit amplitude as for flex calculation
%       - Factors dx/Bx*dy/By for FFT included
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_load_hscyl_itm_tf_point)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [loading] = f_load_hscyl_itm_tf_point(path,plot_cont,calc_cont,loading,dis_itm,geo)

% Initialize: Input
P0_itm_hs  = loading.P0_itm_hs;

dx         = dis_itm.dx;
dy         = dis_itm.dy;
Nx_hs      = dis_itm.Nx_hs;
Ny_hs      = dis_itm.Ny_hs;
kx         = dis_itm.kx;
ky         = dis_itm.ky;
Nf         = dis_itm.Nf;
f          = dis_itm.f;
f_tot      = dis_itm.f_tot;

x          = geo.x;
y          = geo.y;
Bx         = geo.Bx;
By         = geo.By;
nnPiF_ini  = geo.nnPiF_ini;
nnPiF_end  = geo.nnPiF_end;
nnPiF      = 1;

pos_load_x_itm = loading.pos_load_x_itm;
pos_load_y_itm = loading.pos_load_y_itm;

% Creation of load matrix
       
    if P0_itm_hs > 0
        % Predefinition
        P_kxky_f = zeros(Nx_hs,Ny_hs,Nf);
        Px_ITM   = zeros(Nx_hs,3*Ny_hs,Nf);
        Py_ITM   = zeros(Nx_hs,3*Ny_hs,Nf);
        Pz_ITM   = zeros(Nx_hs,3*Ny_hs,Nf);
        P_xy     = zeros(Nx_hs,Ny_hs);
     
        P_xy(x==pos_load_x_itm,y==pos_load_y_itm) = P0_itm_hs;
        P_kxky = dx/Bx*dy/By*ifftshift(fft2(fftshift(P_xy)));
                
        % Predefine constant f spectrum in 3rd dim
        P_f = ones(Nx_hs,Ny_hs,Nf);
        
        % Load const. over f for TF
        P_kxky_f = P_kxky.*P_f;
        
        Px_ITM(:,1:3:3*Ny_hs,:) = P_kxky_f;
        Py_ITM(:,2:3:3*Ny_hs,:) = P_kxky_f;
        Pz_ITM(:,3:3:3*Ny_hs,:) = P_kxky_f;
        
    else
        Px_ITM  = zeros(Nx_hs,3*Ny_hs,Nf,nnPiF_end);
        Py_ITM  = zeros(Nx_hs,3*Ny_hs,Nf,nnPiF_end);
        Pz_ITM  = zeros(Nx_hs,3*Ny_hs,Nf,nnPiF_end);
        
        P_xy    = zeros(Nx_hs,Ny_hs);
        P_kxky  = zeros(Nx_hs,Ny_hs);
    end
    
    %% Output
    
    % For f_displ_hs_cyl
    loading.Px_ITM(:,:,:,nnPiF) = Px_ITM;
    loading.Py_ITM(:,:,:,nnPiF) = Py_ITM;
    loading.Pz_ITM(:,:,:,nnPiF) = Pz_ITM;
     
    
    % For plot_geo
    loading.Pz_ITM_kxky = P_kxky;
    loading.Pz_ITM_xy   = P_xy;
    
    
    %% Plots
    
    if strcmpi(plot_cont.load,'yes') && P0_itm_hs > 0
        
        %% Plot P(x,y), P(kx,ky), P(f)
        %--------------------------------------------
        fig = figure;
        fig.Position = plot_cont.fig_position;
        fig.Name     = 'Hs cyl flex - POINT load';
        %--------------------------------------------
        
        % Plot load over kx,ky RE
        subplot(2,3,1)
        [KX,KY] = ndgrid(kx,ky);
        mesh(KX,KY,real(P_kxky),'FaceColor',plot_cont.FaceColor);
        
        xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        zlabel('$\mathrm{Re} \; P_{z,HsCyl}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        
        % Plot load over kx,ky
        subplot(2,3,2)
        Pplot   = squeeze(P_kxky(Nx_hs/2+1,:));
        plot(ky,real(Pplot));
        
        xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel('$\mathrm{Re} \; P_{z,HsCyl}(k_x=0,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        
        clear Pplot
        
        % Plot P(f)
        subplot(2,3,3)
        Pplot = squeeze(P_kxky_f(Nx_hs/2+1,Ny_hs/2+1,:));
        stem(f,dx*dy*Pplot,'.')
        xlabel('$f$',    'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel('$P(f) \cdot dx \, dy$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        xlim([min(f_tot)-5 max(f_tot)+5])
        
        % Plot of load matrix over (x,y)
        subplot(2,3,4)
        [X,Y] = ndgrid(x,y);
        mesh(X,Y,dx*dy*real(P_xy),'FaceColor',plot_cont.FaceColor)
        xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        zlabel('$\mathrm{Re} \;\; P_{z,HsCyl}(x,y) dx \,dy$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        
        % Plot of load matrix over (x,y)
        subplot(2,3,5)
        Pplot   = dx*dy*squeeze(P_xy(Nx_hs/2+1,:));
        plot(y,real(Pplot));
        xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel('$\mathrm{Re} \; P_{z,HsCyl}(x=0,y) dx \,dy$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        
        % Plot of load matrix over (x,y,z=0)
        subplot(2,3,6)
        [X,Y] = ndgrid(x,y);
        hold on
        plot(X,Y,'.','Color',plot_cont.color_cyan)
        plot(pos_load_x_itm,pos_load_y_itm,'Marker','x','MarkerSize',12,'Color','r');
        
        % xlim([min(x_found)-2 max(x_found)+2])
        % ylim([min(y_found)-2 max(y_found)+2])
        xlabel('$x$','Interpreter','latex');
        ylabel('$y$','Interpreter','latex')
        view([90 -90])
        
     
    end
end

