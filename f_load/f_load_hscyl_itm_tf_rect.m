%__________________________________________________________________________
%
%       f_load_hsycyl_itm_tf_rect.m    
%
%       INPUT: 
%       - path          = pathes
%       - plot_cont     = plot control parameter
%       - loading       = load parameters
%       - dis_itm       = discretization
%       - geo           = geometry
%
%       OUTPUT: 
%       - P_ITM         = Total load vector for ITM subsystem in (kx,ky,f)
%
%       DESCRIPTION: 
%       - Setup rectangular load distribution P(x,y) at z=0 
%       - Adapt load amplitudes at edges and corners such that the sum of the nodal
%         loads is equivalent to the analytic resultant of the rectangular load (volume)
%       - FFT of load P(x,y) → P(kx,ky)
%       - Multiply with unit frequency spectrum in 3rd dimension
%         → unit amplitude P(f) for determination of transfer function (TF)
%       - Assemble load distribution into load matrices 
%         Px_ITM(:,1:3:3*Ny_hs,:) = P_kxkyf;
%         Py_ITM(:,2:3:3*Ny_hs,:) = P_kxkyf;
%         Pz_ITM(:,3:3:3*Ny_hs,:) = P_kxkyf; 
%         → Thereby each Pi_ITM contains the loads in all directions Pzx, Pzy, Pzz at z=0
%         → Load matrix for ITM subsystem of hscyl for rect harmonic load to determine TF
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_load_hscyl_itm_tf_rect)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________


function [loading] = f_load_hscyl_itm_tf_rect(path,plot_cont,loading,dis_itm,geo)

%% Initialize: Input

dx    = dis_itm.dx;
dy    = dis_itm.dy;
Nx_hs = dis_itm.Nx_hs;
Ny_hs = dis_itm.Ny_hs;
Nf    = dis_itm.Nf;
f     = dis_itm.f;
kx    = dis_itm.kx;
ky    = dis_itm.ky;

x     = geo.x;
y     = geo.y;
Bx    = geo.Bx;
By    = geo.By;

P0_itm_hs       = loading.P0_itm_hs;
bx_itm          = loading.bx_itm;
by_itm          = loading.by_itm;
pos_load_x_itm  = loading.pos_load_x_itm;
pos_load_y_itm  = loading.pos_load_y_itm;

%% Creation of load matrix

if P0_itm_hs > 0
        
    % Check if bx and by are multiples of dx
    if mod(bx_itm,dx)~=0; error('Error! The width bx must be a even multiple of dx'); end
    if mod(by_itm,dy)~=0; error('Error! The width by must be a even multiple of dy'); end
    
    % Creation of load matrix

    % Load array vector
    a = zeros(Nx_hs,1);       % x-direction
    c = zeros(Ny_hs,1);       % y-direction
    
    % Number of dx of load block in x resp. y direction
    nbx = bx_itm/dx;
    nby = by_itm/dy;
    
    % Loaded nodes 
    sx = Nx_hs/2+1+(pos_load_x_itm/dx)-nbx:Nx_hs/2+1+(pos_load_x_itm/dx)+nbx; % symmetric to origin (pos_load_x,pos_load_y)  
    sy = Ny_hs/2+1+(pos_load_y_itm/dy)-nby:Ny_hs/2+1+(pos_load_y_itm/dy)+nby; % symmetric to origin (pos_load_x,pos_load_y)
    if size(sx,2)>size(x,2); sx(size(sx,2))=[]; end
    if size(sy,2)>size(y,2); sy(size(sy,2))=[]; end
 
    % Load amplitudes of the loaded nodes 
    a(sx) = P0_itm_hs; % load distribution  in x direction [N/m^2]
    c(sy) = 1;         % Extrusion function in y direction [N/m^2]
    
    % Load matrix on hs in (x,y,z=0,omega)
    P_xy = a*c';
    
    % Adapt load: resultant in Pz(kx=0,ky=0) = Volume of Pz(x,y)
    if bx_itm ~= Bx/2 && by_itm ~= By/2
        
        P_xy(sx(1),sy(1):sy(end))   = 0.5*P_xy(sx(1),sy(1):sy(end));   % first x
        P_xy(sx(end),sy(1):sy(end)) = 0.5*P_xy(sx(end),sy(1):sy(end)); % last x
        
        P_xy(sx(1):sx(end),sy(1))   = 0.5*P_xy(sx(1):sx(end),sy(1));   % first y
        P_xy(sx(1):sx(end),sy(end)) = 0.5*P_xy(sx(1):sx(end),sy(end)); % last y
        
    elseif bx_itm == Bx/2 && by_itm ~= By/2
        
        P_xy(sx(1):sx(end),sy(1))   = 0.5*P_xy(sx(1):sx(end),sy(1));   % first y
        P_xy(sx(1):sx(end),sy(end)) = 0.5*P_xy(sx(1):sx(end),sy(end)); % last y
        
    elseif bx_itm ~= Bx/2 && by_itm == By/2
        
        P_xy(sx(1),sy(1):sy(end))   = 0.5*P_xy(sx(1),sy(1):sy(end));   % first x
        P_xy(sx(end),sy(1):sy(end)) = 0.5*P_xy(sx(end),sy(1):sy(end)); % last x
        
    end
    
    % Fourier Transformation (x,y) -> (kx,ky)
    P_kxky = dx/Bx*dy/By*ifftshift(fft2(fftshift(P_xy)));
    
    % Predefine constant f spectrum in 3rd dim
    P_f = ones(Nx_hs,Ny_hs,Nf);
    
    % Load const. over f for TF
    P_xyf   = P_xy.*P_f;
    P_kxkyf = P_kxky.*P_f;
    
    % Predefinition
    Px_ITM  = zeros(Nx_hs,3*Ny_hs,Nf);
    Py_ITM  = zeros(Nx_hs,3*Ny_hs,Nf);
    Pz_ITM  = zeros(Nx_hs,3*Ny_hs,Nf);
    
    % Sort the vertical loads into the complete load matrix
    Px_ITM(:,1:3:3*Ny_hs,:) = P_kxkyf;
    Py_ITM(:,2:3:3*Ny_hs,:) = P_kxkyf;
    Pz_ITM(:,3:3:3*Ny_hs,:) = P_kxkyf;
    
    % Structure for one omega
    %                   -ky_max             +ky_max
    %               Pzx  Pzy Pzz  | .... | Pzx  Pzy Pzz
    %           -----------------------------------------
    % -kx_max   |                                       |
    %           |                                       |
    %           |                                       |
    % +kx_max   |                                       |
    %           -----------------------------------------
    
    
else

    Px_ITM  = zeros(Nx_hs,3*Ny_hs,Nf);
    Py_ITM  = zeros(Nx_hs,3*Ny_hs,Nf);
    Pz_ITM  = zeros(Nx_hs,3*Ny_hs,Nf);
    
    P_kxkyf = zeros(Nx_hs,Ny_hs,Nf);
    P_xy    = zeros(Nx_hs,Ny_hs);
    
end

%%  Output

% For f_displ_hs_cyl
loading.Px_ITM      = Px_ITM;
loading.Py_ITM      = Py_ITM;
loading.Pz_ITM      = Pz_ITM;

% For plot_geo
loading.Pz_ITM_kxky = P_kxkyf;
loading.Pz_ITM_xy   = P_xy;

%%  Plots
if strcmpi(plot_cont.load,'yes') && P0_itm_hs > 0
    
    %--------------------------------------------
    fig                     = figure;
    fig.Name                = 'HsCyl - Load ITM';
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------
  
    % P_{z,ITM}(x,y)
    subplot(2,2,1)
    hold on
    [X,Y] = ndgrid(x,y);
    Z = zeros(Nx_hs,Ny_hs);
    plot3(X,Y,Z,'Marker','.','Color',plot_cont.color_cyan)
    stem3(X(P_xy~=0),Y(P_xy~=0),P_xy(P_xy~=0),'Marker','.','Color','r');
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    zlabel('$P_{\mathrm{ITM}}(x,y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    view(3)
    
    % P_{ITM}(f)
    subplot(2,2,2)
    Pplot = squeeze(P_xyf(sx(nbx+1),sy(nby+1),:)); % in middle of block load
    stem(f,Pplot,'Marker','.');
    xlabel('$t$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P_{\mathrm{ITM}}(f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    clear Pplot
    ylim([0 P0_itm_hs*1.2])
    xlim([min(f)*1.2 max(f)+0.1])
    
    % P_{ITM}(k_x,k_y)
    subplot(2,2,3)
    [KX,KY] = ndgrid(kx,ky);
    mesh(KX,KY,real(P_kxky),'FaceColor',plot_cont.FaceColor)
    xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    zlabel('$\mathrm{Re} \;\, P_{\mathrm{ITM}}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    clear Pplot
   
    % P_{\mathrm{ITM}}(k_x,f)
    subplot(2,2,4)
    [FF,KX] = ndgrid(f,kx);
    Pplot   = squeeze(P_kxkyf(:,ky==0,:)); % one column for ky = 0
    mesh(FF,KX,real(Pplot).','FaceColor',plot_cont.FaceColor)
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    zlabel('$P_{\mathrm{ITM}}(k_x,f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    clear Pplot
    
   
end
end