%__________________________________________________________________________
%
%       f_displ_itm.m    
%
%       INPUT: 
%
%       OUTPUT: 
%        displ = Displ. hs1L for each triple (kx,ky,omega) in x,y,z direction due to Px,Py,Pz
%                (ux_Px,uy_Px,uz_Px, ux_Py,uy_Py,uz_Py, ux_Pz,uy_Pz,uz_Pz) at z1=0 
%
%       DESCRIPTION: 
%       - Calculate displacements with the fundamental ITM solution for each triple
%         (kx,ky,omega) in x,y,z direction due to Px,Py,Pz 
%          ->(ux_Px,uy_Px,uz_Px, ux_Py,uy_Py,uz_Py, ux_Pz,uy_Pz,uz_Pz) 
%          at z1=0 (=hs surf)           -> uz_0 
%
%       - For f=0 static solution is calculated
%       - static solution not defined for kx=ky=0 -> adapted by small wavenumber shift
%
%       REMARK: 
%       Original Author: Julian Freisinger  (f_displ_itm)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    24-09-2023 - Hicks   
%__________________________________________________________________________



function [displ] = f_displ_itm(dis_itm,~,mat,loading,calc,calc_cont)


%% Predefinition

Nx_hs  = dis_itm.Nx_hs;
Ny_hs  = dis_itm.Ny_hs;
Nf     = dis_itm.Nf;
kx     = dis_itm.kx;
ky     = dis_itm.ky;
f      = dis_itm.f;
omega  = dis_itm.omega;

Pz_kxky_f  = loading.Pz_kxky_f_hs1L;
pos_load_z = loading.pos_z_itm_hs1L; % z=z0 

% prestore
u_z0    = zeros(Nx_hs,Ny_hs,Nf,3);


%% Set variables for loops
switch calc.kx
    case 'kx=0';   nx_start = Nx_hs/2+1;  nx_end   = Nx_hs/2+1;
    case 'kx=all'; nx_start = 1;          nx_end   = Nx_hs;
end

%% Calculation for each triple (kx,ky,omega) 
% (ux_Px,uy_Px,uz_Px,ux_Py,uy_Py,uz_Py,ux_Pz,uy_Pz,uz_Pz)
% - frequency sorting [-fmax ... -f1] with f1 <=0 
Nf_neg = Nf-1;
switch calc_cont.mode
    case 'serial'
            for nnf = 1:Nf-1 % for all frequencies of freq
                for nx = nx_start:nx_end
                    for ny = 1:Ny_hs
                        [u_z0(nx,ny,nnf,:)] = f_displ_hs(kx(nx),ky(ny),omega(nnf),Pz_kxky_f(nx,ny,nnf),mat);
                    end
                end
            end


    case 'parallel'
            % waitbar
            global numCompleted
            numCompleted = 0;               % counter variable for waitbar
            numMax       = dis_itm.Nx_hs;   % Maximum number of iterations in parfor loop
            wait = waitbar(0,'Progress:','Name','Calculate stiffness and displ');
            Nf_calc = Nf;
            D = parallel.pool.DataQueue;
            afterEach(D, @(x) updateWaitbar(numMax,wait,x,Nf_calc));

        for nnf = 1:Nf-1 % for all frequencies of freq
            parfor (nx = nx_start:nx_end,calc_cont.M)
                for ny = 1:Ny_hs
                [u_z0(nx,ny,nnf,:)] = f_displ_hs_paralell(kx(nx),ky(ny),omega(nnf),Pz_kxky_f(nx,ny,nnf),mat);
                end
                send(D, nx)
            end
        end
        close(wait);
        delete(wait)
        set(0,'defaultTextInterpreter','latex'); 
end

    

%% Output
% z0
displ.TF.uxPz_kxkyf_z0 = u_z0(:,:,:,1);
displ.TF.uyPz_kxkyf_z0 = u_z0(:,:,:,2);
displ.TF.uzPz_kxkyf_z0 = u_z0(:,:,:,3);

end
