%__________________________________________________________________________
%
%       f_load_itm_tf_rect.m    
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
%       - Create load matrix in (x,y,z=0,omega) domain with rect block with midpoint 
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
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________



function [loading] = f_load_itm_tf_rect(path,calc,plot_cont,loading,dis_itm,geo)

% Initialize: Input

P0_hs1L     = loading.P0_itm_hs1L;
bx_hs1L     = loading.bx_itm_hs1L;
by_hs1L     = loading.by_itm_hs1L;
pos_x_hs1L  = loading.pos_x_itm_hs1L;
pos_y_hs1L  = loading.pos_y_itm_hs1L;

dx          = dis_itm.dx;
dy          = dis_itm.dy;
Nx_hs       = dis_itm.Nx_hs;
Ny_hs       = dis_itm.Ny_hs;  
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
    
    %% Check if bx and by are multiples of dx
    if mod(bx_hs1L,dx)~= 0; error('Error! The width bx must be a even multiple of dx'); end
    if mod(by_hs1L,dy)~= 0; error('Error! The width by must be a even multiple of dy'); end
    
    %% Spatial load distribution
    
    % Load array vector
    a = zeros(Nx_hs,1);       % x-direction
    c = zeros(Ny_hs,1);       % y-direction
    
    % Number of dx of load block in x resp. y direction
    nbx = bx_hs1L/dx;
    nby = by_hs1L/dy;
    
    % Loaded nodes in load function a symmetric to origin (pos_load_x,pos_load_y)
    sx = Nx_hs/2+1+(pos_x_hs1L/dx)-nbx:Nx_hs/2+1+(pos_x_hs1L/dx)+nbx;
    if size(sx,2)>size(x,2), sx(size(sx,2)) = [];  end
    
    % Loaded nodes in load function c symmetric to origin (pos_load_x,pos_load_y)
    sy = Ny_hs/2+1+(pos_y_hs1L/dy)-nby:Ny_hs/2+1+(pos_y_hs1L/dy)+nby;
    if size(sy,2)>size(y,2), sy(size(sy,2)) = []; end
    
     % Set amplitude
    P0 = P0_hs1L;

    % Load amplitudes of loaded nodes (= load distribution)
    a(sx) = P0;  % in x direction [N/m^2]
    c(sy) = 1;   % in y direction via extrusion function
    
    % Load matrix on hs surf in (x,y,z = 0,omega)
    Pz_xy = a*c.';
    
    % Adapt load: resultant in Pz(kx = 0,ky = 0) Volume of Pz(x,y)
    if bx_hs1L ~= Bx/2 && by_hs1L ~= By/2 % for block load 
        
        Pz_xy(sx(1),sy(1):sy(end))    =  0.5*Pz_xy(sx(1),sy(1):sy(end));   % first x
        Pz_xy(sx(end),sy(1):sy(end))  =  0.5*Pz_xy(sx(end),sy(1):sy(end)); % last x
        
        Pz_xy(sx(1):sx(end),sy(1))    =  0.5*Pz_xy(sx(1):sx(end),sy(1));   % first y
        Pz_xy(sx(1):sx(end),sy(end))  =  0.5*Pz_xy(sx(1):sx(end),sy(end)); % last y
        
    elseif bx_hs1L == Bx/2 && by_hs1L ~= By/2 % for Line load
        
        Pz_xy(sx(1):sx(end),sy(1))    =  0.5*Pz_xy(sx(1):sx(end),sy(1));   % first y
        Pz_xy(sx(1):sx(end),sy(end))  =  0.5*Pz_xy(sx(1):sx(end),sy(end)); % last y
        
    elseif bx_hs1L ~= Bx/2 && by_hs1L == By/2 % for Line load
        
        Pz_xy(sx(1),sy(1):sy(end))    =  0.5*Pz_xy(sx(1),sy(1):sy(end));   % first x
        Pz_xy(sx(end),sy(1):sy(end))  =  0.5*Pz_xy(sx(end),sy(1):sy(end)); % last x
        
    end
    
    %% Time dependence -> 3d matrix
    
    % Constant frequency spectrum 
    % -> unit load for each frequency = constant freq spectrum
    Pz_f = ones(Nx_hs,Ny_hs,Nf);

    % Spatial FFT (x,y) -> (kx,ky)
    Pz_kxky = dx/Bx*dy/By*fftshift(fftn(fftshift(Pz_xy)));
    % Constant frequency spectrum 
    Pz_kxky_f = Pz_kxky.*Pz_f;

    
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

    %-----------------------------------------------------
    mainfigurename = 'Hs1L - ITM rect load P(x,y), P(kx,ky) at z=0';
    fig = figure('Name',mainfigurename);
    fig.Position = [632,260,1403,1040];
    %-----------------------------------------------------
    sgtitle(['Spatial distribution of load $P(x,y,z=0)$ and $P(k_x,k_y,z=0)$ for $\mathrm{TF}_{u_{\mathrm{soil}}}$'],'interpreter','latex','FontSize',12)
%
    % P(x,y,z=0,omega)
    subplot(2,3,1)
    [X,Y] = ndgrid(x(1:4:end),y(1:4:end));
    mesh(X,Y,Pz_xy(1:4:end,1:4:end),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    zlabel('$P_{z}(x,y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)

    % P(kx,ky,z=0,omega)
    subplot(2,3,2)
    [KX,KY] = ndgrid(kx(1:4:end),ky(1:4:end));
    mesh(KX,KY,real(Pz_kxky(1:4:end,1:4:end)),'FaceColor',plot_cont.FaceColor)
    xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    zlabel('$\mathrm{Re} \; P_{z}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)

    % P(f)
    subplot(2,3,3)
    Pplot = squeeze(Pz_f(Nx_hs/2+1,Ny_hs/2+1,:));
    plot(f,Pplot)
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{z}(f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    xlim([min(f)-1 max(f)+1])

    % P(x,y,z=0,omega)
    subplot(2,3,4)
    plot(y,Pz_xy(Nx_hs/2+1,:));
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{z}(x=0,y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)

    % P(kx=0,ky,z=0,omega)
    subplot(2,3,5)
    Pplot = real(squeeze(Pz_kxky_f(Nx_hs/2+1,:,1)));
    %Pplot = squeeze(Pz_kxky_f(Nx_hs/2+1,:,1));
    plot(ky,Pplot)
    xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{z}(k_x=0,k_y,f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    clear Pplot

    % P(kx=0,ky=0,z=0,f)
    subplot(2,3,6)
    Pplot = squeeze(Pz_kxky_f(Nx_hs/2+1,Ny_hs/2+1,:));
    stem(f,Pplot,'.')
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{z}(k_x=k_y=0,f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    xlim([min(f)-1 max(f)+1])
    clear Pplot

    savename = ['\load_xyf_TF'];
    saveas(fig,[path.figures savename],'fig')
    %
    if strcmp(plot_cont.export,'Yes')
        exportfig(gcf,[path.figures savename],'Width',plot_cont.exp_width,'Height',plot_cont.exp_heigth,'Fontmode','fixed','Fontsize',plot_cont.Fontsize,'Color','cmyk','resolution',plot_cont.resolution);
    end

end

end