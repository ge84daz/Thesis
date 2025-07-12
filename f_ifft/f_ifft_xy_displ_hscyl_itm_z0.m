%__________________________________________________________________________
%
%        f_ifft_xy_displ_hscyl_itm_z0.m    
%
%       INPUT: 
%       calc system   = current calc system
%       displ_dir     = ux,uy,uz
%       load_dir      = Px,Py,Pz
%       plot_cont     = plot control parameter
%       displ         = displ from ITM-FEM calc for hscyl
%       dis_itm       = ITM discretization
%       geo           = geometry parameter
%       calc_cont     = calc control parameter
%
%       OUTPUT: 
%       -uiPi_xyf_hscyl_itm_z0 = displ at z=0 for hscyl from itm solution in u(x,y,f) in displ dir i due to load Pi
%
%       DESCRIPTION: 
%       - Select displ ui_kxkyf_ges from solution of hscyl for current system and load dir Pi
%       - Conj. complex expansion for f>0
%       - Determine IFFT at z=0 from (kx,ky,f) to (x,y,f) for uiPi
% 
%       REMARKS: 
%       - displ. at z=0 calculated for ALL (x,y) at virtual hs surface
%       - BUT: actually equilibrium only on the physically existing system
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_ifft_xy_displ_hscyl_itm_z0)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________

function [displ] = f_ifft_xy_displ_hscyl_itm_z0(displ_dir,load_dir,plot_cont,displ,dis_itm,geo,calc_cont,calc)

%%  Initialize: Input

Nx_hs  = dis_itm.Nx_hs;
Ny_hs  = dis_itm.Ny_hs;
dx     = dis_itm.dx;
dy     = dis_itm.dy;
f_tot  = dis_itm.f_tot;
Nf     = dis_itm.Nf;
Nf_tot = dis_itm.Nf_tot;
fev    = dis_itm.fev;

Bx    = geo.Bx;
By    = geo.By;
x     = geo.x;
y     = geo.y;

u_xyf_hs_z0 = zeros(Nx_hs,Ny_hs,Nf_tot);


%% Select displacements u_kx_ky_omega_ges for current system and load direction
u_kxkyf_tot = getfield(displ.hscyl,calc.system,'tot',['u' load_dir '_kxkyf_hscyl_tot']);
    
% Select displ. on hs surface z=0 depending on displ_dir
switch displ_dir
    case 'ux'; u_kxkyf_hscyl_z0 = squeeze(u_kxkyf_tot(:,1:3:3*Ny_hs,:));
    case 'uy'; u_kxkyf_hscyl_z0 = squeeze(u_kxkyf_tot(:,2:3:3*Ny_hs,:));
    case 'uz'; u_kxkyf_hscyl_z0 = squeeze(u_kxkyf_tot(:,3:3:3*Ny_hs,:));
end

%% Conj. complex expansion for f>0 (vgl. homogeneous halfspace)
% Expand matrix for f > 0 => TF(kx,ky,z0,f_all)
% flipdim cp. Müller 1989 p27-29 and Hackenberg 2016 p.18
u_kxkyf_hscyl_z0(:,:,Nf+1:Nf_tot)             = conj(u_kxkyf_hscyl_z0(:,:,1:Nf-1));                        % Conj complex
u_kxkyf_hscyl_z0(2:Nx_hs,2:Ny_hs,Nf+1:Nf_tot) = flip(u_kxkyf_hscyl_z0(2:Nx_hs,2:Ny_hs,Nf+1:Nf_tot),1);     % flip dimensions
u_kxkyf_hscyl_z0(2:Nx_hs,2:Ny_hs,Nf+1:Nf_tot) = flip(u_kxkyf_hscyl_z0(2:Nx_hs,2:Ny_hs,Nf+1:Nf_tot),2);
u_kxkyf_hscyl_z0(1,2:Ny_hs,Nf+1:Nf_tot)       = flip(u_kxkyf_hscyl_z0(1,2:Ny_hs,Nf+1:Nf_tot),2);
u_kxkyf_hscyl_z0(2:Nx_hs,1,Nf+1:Nf_tot)       = flip(u_kxkyf_hscyl_z0(2:Nx_hs,1,Nf+1:Nf_tot),1);
u_kxkyf_hscyl_z0(1:Nx_hs,1:Ny_hs,Nf+1:Nf_tot) = flip(u_kxkyf_hscyl_z0(1:Nx_hs,1:Ny_hs,Nf+1:Nf_tot),3);


%% IFFT u(z=0)
u_xkyf_hs_z0             = Bx/dx*fftshift(ifft(ifftshift(u_kxkyf_hscyl_z0),[],1));  % kx to x (-> over columns)
u_xyf_hs_z0(:,:,:) = By/dy*fftshift(ifft(ifftshift(u_xkyf_hs_z0), [],2));  % ky to y (-> over rows)
%   u_xyt_hs_z0              = T/dt*fftshift (ifft(ifftshift(u_xyf_hs_z0),  [],3));  % omega to t (-> over 3rd dim)

%%  Output

varname1 = [displ_dir load_dir '_xyf_hscyl_itm_z0'];
displ.hscyl.(calc.system).itm.(varname1) = u_xyf_hs_z0;

%%  Plots

if strcmpi(plot_cont.displ_trans,'yes')
    
    % Plots the displacements on the complete surface for all x,y at z=0.
    % Actually only displacements on surface z=0 outside the trench can be
    % calculated with ITM. The displ at z=0 within the trench have to
    % be taken out of the FEM subsystem. However, displ at z=0 calculated for all
    % x,y only due to superposition procedure of fundamental solutions.
    
    fig             = figure;
    fig.Name        = ['ITM-FEM - Displ. uz from ITM solution: u(x,y,z=0,fev) - ' calc.system];
    fig.Position    = plot_cont.fig_position;
    %
    [X,Y] = ndgrid(x,y);
    mesh(X,Y,abs(u_xyf_hs_z0(:,:,f_tot==fev)),'FaceColor',plot_cont.FaceColor,'EdgeColor','interp');
    %
    if strcmpi(plot_cont.map,'yes');     colormap(plot_cont.colormap);   end
    if strcmpi(plot_cont.shading,'yes'); shading(plot_cont.shading_typ); end
    %
    xlim([-Bx/2 Bx/2-dx]);
    ylim([-By/2 By/2-dy]);
    
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$|u_z(x,y,z=0,f=' num2str(fev) ')|$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    
    if strcmpi(plot_cont.leg,'yes')
        leg=legend('$|u_z(x,y,z=0,f=f_{ev}|$','Location',plot_cont.loc_ledg,'Interpreter',plot_cont.interpreter);
        set(leg,'Box',plot_cont.leg_box );
    end
end
end
