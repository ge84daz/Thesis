%__________________________________________________________________________
%
%       f_ifft_xy_displ_TF.m    
%
%       INPUT: 
%       displ = Displ. hs1L for each triple (kx,ky,omega) in x,y,z direction due to Px,Py,Pz
%              (ux_Px,uy_Px,uz_Px, ux_Py,uy_Py,uz_Py, ux_Pz,uy_Pz,uz_Pz)
%
%       OUTPUT: 
%       u_x_y_omega_z0 = Displ. in (x,y,z1=0,omega) at halfspace surface
%       u_x_y_omega_zH = Displ. in (x,y,z2=0,omega) at interface layer-halfspace
%
%       DESCRIPTION: 
%       Spatial IFFT of displ for constant f spectrum from (kx,ky,z,omega) to (x,y,z,omega)
%
%       - hs, hs1L  -> calc displ of hs or hs1L at z=0 and z=H
%
%       - For interp xy large RAM needed -> remedy: loop over frequencies and
%       IFFT kx,ky -> x,y for each f seperately and sort into matrix.
%
%       REMARK: 
%       Original Author: Julian Freisinger  (f_ifft_xy_displ_TF)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________


function [displ] = f_ifft_xy_displ_TF(displ_dir,load_dir,displ,dis_itm,geo,loading,calc)

%%  Initialize: Input

Nx_hs  = dis_itm.Nx_hs;
Ny_hs  = dis_itm.Ny_hs;
Nf     = dis_itm.Nf;
Nf_tot = dis_itm.Nf_tot;
f_tot  = dis_itm.f_tot;
dx     = dis_itm.dx;
dy     = dis_itm.dy;

Bx     = geo.Bx;
By     = geo.By;

% load displ
varname1 = [displ_dir load_dir '_kxkyf_z0'];

u_kxkyf_z0 = displ.TF.(varname1);

%% Conj complex expansion and flip u(z1=0)
% Expand matrix for f > 0 => TF(kx,ky,z0,f_all)
% flip cp. Müller 1989 p27-29 and Hackenberg 2016 p.18

u_kxkyf_z0(:,:,Nf+1:Nf_tot)             = conj(u_kxkyf_z0(:,:,1:Nf-1));                        % Conj complex
u_kxkyf_z0(2:Nx_hs,2:Ny_hs,Nf+1:Nf_tot) = flip(u_kxkyf_z0(2:Nx_hs,2:Ny_hs,Nf+1:Nf_tot),1);     % flip dimensions
u_kxkyf_z0(2:Nx_hs,2:Ny_hs,Nf+1:Nf_tot) = flip(u_kxkyf_z0(2:Nx_hs,2:Ny_hs,Nf+1:Nf_tot),2);
u_kxkyf_z0(1,2:Ny_hs,Nf+1:Nf_tot)       = flip(u_kxkyf_z0(1,2:Ny_hs,Nf+1:Nf_tot),2);
u_kxkyf_z0(2:Nx_hs,1,Nf+1:Nf_tot)       = flip(u_kxkyf_z0(2:Nx_hs,1,Nf+1:Nf_tot),1);
u_kxkyf_z0(1:Nx_hs,1:Ny_hs,Nf+1:Nf_tot) = flip(u_kxkyf_z0(1:Nx_hs,1:Ny_hs,Nf+1:Nf_tot),3);



%% IFFT (kx,ky)→(x,y) z=0
u_xkyf_z0 = Bx/dx*fftshift(ifft(ifftshift(u_kxkyf_z0,1),[],1),1);  % IFFT kx    -> x (over columns)
u_xyf_z0  = By/dy*fftshift(ifft(ifftshift(u_xkyf_z0,2), [],2),2);  % IFFT ky    -> y (over rows)



%%  Output

varname1 = [displ_dir load_dir '_xyf_z0'];
varname2 = [displ_dir load_dir '_kxkyf_z0'];

displ.TF.(varname1) = u_xyf_z0;
displ.TF.(varname2) = u_kxkyf_z0;

end

