%__________________________________________________________________________
%
%        f_ifft_xy_displ_hscyl_fem.m    
%
%       INPUT: 
%       - calc.system   = current calc system
%       - load_dir      = Px,Py,Pz
%       - calc_cont     = calc control parameter
%       - plot_cont     = plot control parameter
%       - displ         = displ from ITM-FEM calc for hscyl
%       - dis_itm       = ITM discretization
%       - dis_fem       = FEM discretization
%       - geo           = geometry parameter
%
%       OUTPUT: 
%       - u_xyf_fem_xyz         = all FEM displ ux,uy,uz of cyl (boundary + inside) 
%       - u_kxyf_fem_xyz        = all FEM displ ux,uy,uz of cyl (boundary + inside) 
%       - u_x_midline_f_fem_x   = displ. ux on cyl midline
%       - u_y_midline_f_fem_y   = displ. uy on cyl midline
%       - u_z_midline_f_fem_z   = displ. uz on cyl midline
%
%       DESCRIPTION: 
%       - Select displ ui_kxkyf_ges from solution of hscyl for fem subsystem and load dir Pi
%       - Inverse Transformation of displ on the cyl boundary
%       - Conj. complex expansion for f>0 for displ inside FEM subsystem
%       - IFFT displ inside FEM subsystem (kx,y,f) to (x,y,f) for uxyz_f due to Pi
%       - Determine IFFT at z=0 from (kx,ky,f) to (x,y,f) for uiPi
%       - load_dir input referes to the load direction of Pi of the ITM subsystem
%       - load_dir for FEM is given at beginning in input
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_ifft_xy_displ_hscyl_fem)
%       Modified:   Tom Hicks                 
%       Date:       25-09-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________


function [displ] = f_ifft_xy_displ_hscyl_fem(load_dir,tun_sel,calc,calc_cont,plot_cont,displ,dis_itm,dis_fem,geo)
   
% Initialize: Input

Nx_hs    = dis_itm.Nx_hs;
Ny_hs    = dis_itm.Ny_hs;
Nphi_t   = dis_itm.Nphi_t;

Nf       = dis_itm.Nf;
Nf_tot   = dis_itm.Nf_tot;
f_tot    = dis_itm.f_tot;
dx       = dis_itm.dx;
dy       = dis_itm.dy;
fev      = dis_itm.fev;

Bx       = geo.Bx;
By       = geo.By;
x        = geo.x;


DOFS            = dis_fem.DOFS;
T12             = dis_fem.T12;
num_nod         = dis_fem.num_nod;
num_nod_in      = dis_fem.num_nod_in;
nod_midline     = dis_fem.nod_midline;
num_nod_midline = dis_fem.num_nod_midline;



% Predefinitions
ux_xyf_fem = zeros(Nx_hs,num_nod,Nf_tot);
uy_xyf_fem = zeros(Nx_hs,num_nod,Nf_tot);
uz_xyf_fem = zeros(Nx_hs,num_nod,Nf_tot);
ux_x_midline_f_fem  = zeros(Nx_hs,num_nod_midline,Nf_tot);
uy_x_midline_f_fem  = zeros(Nx_hs,num_nod_midline,Nf_tot);
uz_x_midline_f_fem  = zeros(Nx_hs,num_nod_midline,Nf_tot);
u_xyf_fem_xyz_nnPiF  = zeros(Nx_hs,size(DOFS,1)*3,Nf_tot);
u_kxyf_fem_xyz_nnPiF = zeros(Nx_hs,size(DOFS,1)*3,Nf_tot);

    
%% Select displacements u_kx_ky_omega_ges for current system and load direction
u_kxkyf_tot = getfield(displ.hscyl,calc.system,'tot',['u' load_dir '_kxkyf_hscyl_tot']);
    
%% Inverse Transformation of displ on the cyl boundary
% u_fem_Gamma(kx,r,n,omega) -> u_fem_Gamma(kx,y,z,omega)
% Choosing all DOF at the boundary of the trench and inside
% Choosing Displacements for all Nx_hs, Nodes of FEM and all frequencies
u_kxyf_fem_all = squeeze(u_kxkyf_tot(:,  3*(Ny_hs)+1:3*(Ny_hs+Nphi_t+num_nod_in),:));
       
% Choosing only boundary elements for transformation from (kx,r,n,omega) to (kx,y,z,omega)
for nnf = 1:Nf
    % Select displ on cyl boundary
    u_rand_r_n = u_kxyf_fem_all(:,1:3*Nphi_t,nnf).';
        
    % Transformation of boundary displacements (kx,r,n,omega) to (kx,y,z,omega)
    u_rand_y_z = T12*u_rand_r_n;
        
    % Replacing of old boundary displacements in global displacement vector
    %      Gamma             Omega
    % [ (kx,y,z,omega)   (kx,y,z,omega)  ]
    u_kxyf_fem_all(:,1:3*Nphi_t,nnf) = u_rand_y_z.';
      
end
    
    %% Conj Complex expansion of displ inside FEM subsystem
    
    switch calc_cont.mesh_type
        case 'full_cyl'
            u_kxyf_fem_xyz = u_kxyf_fem_all;
        case 'half_cyl'
            % DOFs of not existing nodes are skipped
            DOF_hc_xyz     = sortrows([DOFS(:,2); DOFS(:,3); DOFS(:,4)]);
            u_kxyf_fem_xyz = u_kxyf_fem_all(:,DOF_hc_xyz,:);
    end
    
    % Expand matrix for f > 0 => TF(kx,ky,z0,f_all)
    % flipdim cp. Müller 1989 p27-29 and Hackenberg 2016 p.18
    % Remarks:
    % - Before only calculated for negative f
    % - Static excitation not included here
    % - Here: Evaluation for FEM mesh, therefore columns not sorted for ky but
    %   for node IDs. Thus also mirroring of matrix members only in 1st and 3rd
    %   dimension. y doesn't have to be flipped because it wasn't transformed.
    u_kxyf_fem_xyz(1:Nx_hs,:,Nf+1:Nf_tot) = conj(u_kxyf_fem_xyz(1:Nx_hs,:,1:Nf-1));
    u_kxyf_fem_xyz(2:Nx_hs,:,Nf+1:Nf_tot) = flip(u_kxyf_fem_xyz(2:Nx_hs,:,Nf+1:Nf_tot),1); % flip over kx for 2:Nx_hs
    u_kxyf_fem_xyz(2:Nx_hs,:,Nf+1:Nf_tot) = flip(u_kxyf_fem_xyz(2:Nx_hs,:,Nf+1:Nf_tot),3); % flip over Nt for 2:Nx_hs
    u_kxyf_fem_xyz(1,:,Nf+1:Nf_tot)       = flip(u_kxyf_fem_xyz(1,      :,Nf+1:Nf_tot),3); % flip over Nt for -kxmax
    
       
    %% IFFT displ inside FEM subsystem
        % hscyl: ITM load P(kx,ky) but IFFT in FEM subsystem only w.r.t. x
        % ITM load hscyl: 
        % flex 'yes' -> P_kxky =             fft(P_xy) -> u_xy =       dy/By*ifft(u_kxy)
        %            -> factor dy/By needs to be included in u_xy explicitly
        %            -> no factor Bx/dx as no dx/Bx in forward transformation
        % flex 'no'  -> P_kxky = dx/Bx*dy/By*fft(P_xy) -> u_xy = Bx/dx*ifft(u_kxy)
        %            -> factor dy/By is included in u_xy implicitly
        %            -> actually also factor By/dy in ifft, but not applied
        %               as dy/By is needed in FEM displ 
        %                          
        % for calc_cont.flex = 'no' factors in FFT of load dx/Bx used -> Bx/dx for IFFT
        u_xyf_fem_xyz   = Bx/dx*fftshift(ifft(ifftshift(u_kxyf_fem_xyz),[],1)); % IFFT kx → x
    
    % all displ within fem in u(x,y,z,f) & u(kx,y,z,f) for all load positions
    % Tunnel 1
    u_kxyf_fem_xyz_nnPiF(:,:,:) = u_kxyf_fem_xyz;
    u_xyf_fem_xyz_nnPiF(:,:,:)  = u_xyf_fem_xyz;
    
    % Select displacements ux, uy, uz
    % Tunnel 1
    ux_xyf_fem(:,:,:) = u_xyf_fem_xyz(:,1:3:end,:); % ux
    uy_xyf_fem(:,:,:) = u_xyf_fem_xyz(:,2:3:end,:); % uy
    uz_xyf_fem(:,:,:) = u_xyf_fem_xyz(:,3:3:end,:); % uz
        
    % Select displacements on midline of cyl. cavity
    % sort into 3D array w.r.t.(kx,nodmidline,f)
    % Tunnel 1
    ux_x_midline_f_fem(:,:,:) = ux_xyf_fem(:,nod_midline(:,1),:); % ux(x,y,z,omega)
    uy_x_midline_f_fem(:,:,:) = uy_xyf_fem(:,nod_midline(:,1),:); % uy(x,y,z,omega)
    uz_x_midline_f_fem(:,:,:) = uz_xyf_fem(:,nod_midline(:,1),:); % uz(x,y,z,omega)


%% Output
% Tunnel 1
varname_fem_all_xyz = {['uxyz' load_dir '_xyf_hscyl_fem_all'],['uxyz' load_dir '_kxyf_hscyl_fem_all']};
varname_fem_all     = {['ux' load_dir '_xyf_hscyl_fem_all'],    ['uy' load_dir '_xyf_hscyl_fem_all'],    ['uz' load_dir '_xyf_hscyl_fem_all']};
varname_fem         = {['ux' load_dir '_xyf_hscyl_fem_midline'],['uy' load_dir '_xyf_hscyl_fem_midline'],['uz' load_dir '_xyf_hscyl_fem_midline']};

displ.hscyl.(calc.system).fem.(tun_sel).(varname_fem_all{1}) = ux_xyf_fem;
displ.hscyl.(calc.system).fem.(tun_sel).(varname_fem_all{2}) = uy_xyf_fem;
displ.hscyl.(calc.system).fem.(tun_sel).(varname_fem_all{3}) = uz_xyf_fem;

displ.hscyl.(calc.system).fem.(tun_sel).(varname_fem{1}) = ux_x_midline_f_fem;
displ.hscyl.(calc.system).fem.(tun_sel).(varname_fem{2}) = uy_x_midline_f_fem;
displ.hscyl.(calc.system).fem.(tun_sel).(varname_fem{3}) = uz_x_midline_f_fem;

% displ.hscyl.(calc.system).fem.(tun_sel).(varname_fem_all_xyz{1}) = u_xyf_fem_xyz_nnPiF;
displ.hscyl.(calc.system).fem.(tun_sel).(varname_fem_all_xyz{2}) = u_kxyf_fem_xyz_nnPiF;


%% Plots
if strcmpi(plot_cont.displ_trans,'yes')
    
    % Plot of the disp in original domain
    %----------------------------------------------------------------------------
    fig                     = figure;
    fig.Name                = ['ItmFem - uz only FEM in (x,y,z,f) - ' calc.system];
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %----------------------------------------------------------------------------
    %
    subplot(1,2,1)
    [X,Y] = ndgrid(x,nod_midline(:,2));
    surf(X,Y,abs(uz_x_midline_f_fem(:,:,f_tot==-fev)),'FaceColor',plot_cont.FaceColor);
    if strcmpi(plot_cont.map,'yes');     colormap(plot_cont.colormap);   end
    if strcmpi(plot_cont.shading,'yes'); shading(plot_cont.shading_typ); end
    %
    title('$u_z$ within tunnel over $x$ and $y$','Interpreter',plot_cont.interpreter)
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel('$|u_z(x,y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    %
    subplot(1,2,2)
    plot(nod_midline(:,2),abs(uz_x_midline_f_fem(size(x,2)/2+1,:,f_tot==-fev)),'color',plot_cont.color_cyan,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth);
    
    xlabel('$y$','Interpreter',plot_cont.interpreter)
    ylabel('$|u_z(y)|$','Interpreter',plot_cont.interpreter)
    
    % Set Figure and axis handles
    set(gca,'XDir','reverse');
    set(gca,'FontSize',plot_cont.Fontsize)                                  % Fontsize axis
    set(findall(gca, 'Type', 'Line'),'LineWidth',plot_cont.LineWidth);      % Linewidth all lines
    
end
end

