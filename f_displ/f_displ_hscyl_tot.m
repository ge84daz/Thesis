%__________________________________________________________________________
%
%        f_displ_hscyl_tot.m    
%
%       INPUT: 
%       - dis_itm       = Discretization ITM
%       - dis_fem       = Discretization FEM
%       - geo           = Geometry
%       - loading       = Load (amplitudes, distribution,...)
%       - mat           = Material Parameters
%       - calc_cont     = Calculation control parameters
%       - calc          = Calculation parameters
%
%       OUTPUT: 
%       - displ.tot   for different systems tunnel, trench, slit, hs_cyl_1L
%       - dis_fem.T   Transformation matrices for FEM part -> used in IFFT_fem_in
%
%       calc.system: TUNNEL -> u_kx_ky_omega_ges_xxx:
%       1.) Column 1                                           : 3* Ny_hs                                                      
%       2.) Column 3*(Ny_hs)+1                                 : 3*(Ny_hs+Nphi_t_T1)                                           
%       3.) Column 3*(Ny_hs+Nphi_t_T1)+1                       : 3*(Ny_hs+Nphi_t_T1+num_nod_in_T1)                             
%
% 
%       1.) Displ. on ground surface z=0 (LAMBDA_1)
%       2.) Displ. on cyl. coupling surface T1 (GAMMA_1)
%       3.) Displ. within cyl. inclusion T1 (=FEM domain) (OMEGA_1)
%
%       DESCRIPTION: 
%       Calculation of the displacements of the coupled ITM-FEM system in 
%       the wavenumber frequency domain (kx,ky,omega) for each combination 
%       of (kx,Omega) in each case for all ky.
%
%       1.) Transformation Matrix T
%       2.) Stiffness Matrix FEM 
%       3.) Stiffness Matrix ITM
%       4.) Assemble total stiffness
%       5.) Assemble total load vector
%       6.) Solve System
% 
%       REMARKS: 
%       - Generally omega-loop only over negative frequencies -> nnf=1:Nf
%       - For all frequencise dis_itm.freq
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_displ_hscyl_tot)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [displ,dis_fem_T1] = f_displ_hscyl_tot(dis_itm,dis_fem_T1,geo_T1,loading,mat,displ,trunc_T1,calc_cont,calc)

%% Initialize: Input 

% dis_itm
Nx_hs      = dis_itm.Nx_hs;
Ny_hs      = dis_itm.Ny_hs;
Nphi_t_T1  = dis_itm.Nphi_t;
Nf         = dis_itm.Nf;
omega      = dis_itm.omega;
kx         = dis_itm.kx;

% Tunnel 1
H_T1          = geo_T1.H;
radius_T1     = geo_T1.radius;
Edof_T1       = dis_fem_T1.Edof;
NDOF_T1       = dis_fem_T1.NDOF;
num_nod_in_T1 = dis_fem_T1.num_nod_in;
div           = dis_fem_T1.div; 

% Parameters for passing in f_stiff_fem_mex
P_FE    = loading.P_FE;
% Px_ITM  = loading.Px_ITM;
% Py_ITM  = loading.Py_ITM;
Pz_ITM  = loading.Pz_ITM;
P_FE_T1 = loading.P_FE_T1;

calc_struct = calc.structure;

% Tunnel T1
num_elem_T1             = dis_fem_T1.num_elem;
elem_coord_T1           = dis_fem_T1.elem_coord;

% material for passing in f_stiff_fem_mex
D_fem_soil_T1   = mat.fem_soil_T1.D;
rho_fem_soil_T1 = mat.fem_soil_T1.rho;

%% Predefinitions 
switch calc.system
    case 'tunnel'
            % uPx_kxkyf_tot_tun = complex(zeros(Nx_hs,3*(Ny_hs+Nphi_t_T1+num_nod_in_T1),Nf));
            % uPy_kxkyf_tot_tun = complex(zeros(Nx_hs,3*(Ny_hs+Nphi_t_T1+num_nod_in_T1),Nf));
            uPz_kxkyf_tot_tun = complex(zeros(Nx_hs,3*(Ny_hs+Nphi_t_T1+num_nod_in_T1),Nf));  
    case 'trench'
            % uPx_kxkyf_tot_tr = complex(zeros(Nx_hs,3*(Ny_hs+Nphi_t_T1+num_nod_in_T1),Nf));
            % uPy_kxkyf_tot_tr = complex(zeros(Nx_hs,3*(Ny_hs+Nphi_t_T1+num_nod_in_T1),Nf));
            uPz_kxkyf_tot_tr = complex(zeros(Nx_hs,3*(Ny_hs+Nphi_t_T1+num_nod_in_T1),Nf));
end

switch calc_cont.kx
    case 'kx=0'                         % 2.5D case
        nx_start   = Nx_hs/2+1;
        nx_end     = Nx_hs/2+1;

    case 'kx=all'                       % 3D case
        nx_start = 1;
        nx_end      = length(kx);
end

% Loop parameter for frequency 
% till now only consider one frequency for harmonic load case
nf_start = 1;
nf_end   = 1;
Nf_calc  = 1;

%% Waitbar
global waitbar_system numCompleted
waitbar_system = 'ItmFem';
numCompleted   = 0;        % counter variable for waitbar
numMax         = Nx_hs*Nf_calc;    % Maximum number of iterations in parfor loop

% Initialize waitbar
wait = waitbar(0,'Progress:','Name','Calculate stiffness and displ ItmFem');

% Setup data queue
D = parallel.pool.DataQueue;
afterEach(D, @(x) updateWaitbar(numMax,wait,x));


%%  Calculate Transformation Matrix T -> same for all (kx, omega)
[dis_fem_T1] = f_trans_matrix_fem(dis_itm,dis_fem_T1,geo_T1,calc_cont); % Tunnel T1
   
warning('off','all')
%%  Calculate stiffness, assemble, solve for displacements - SERIAL
switch calc_cont.mode 
    case 'serial'
    for nnf = nf_start:nf_end % for all f of freq without f=0
        
        % parfor (nnx = nx_start:nx_end,calc_cont.M) % - no parallel computation implemented so far
        for nnx = nx_start:nx_end 
            omega_calc = omega(nnf);
    
            if omega_calc/2/pi <0                                   % f<0
                % calc displ for neg f
            elseif omega_calc /2/pi>=-0.1 && omega_calc/2/pi <=+0.1 % f~0
                % As no static solution is included 
                % -> approximate f=0 by solution for very low frequency
                omega_calc = -0.5*2*pi;
            elseif omega_calc > 0.1                                 % f>0
                % calc displ for neg f and take conj complex of result
                omega_calc = -omega_calc;
            end
    
            % Calculate stiffness matrix K_fem for one combination (kx,Omega)
            switch calc.system
            case {'tunnel'} 
                K_fem_tunnel_T1 = f_stiff_fem_cyl('tunnel',calc_struct,num_elem_T1,Edof_T1,div,NDOF_T1,elem_coord_T1,D_fem_soil_T1,rho_fem_soil_T1,H_T1,radius_T1,omega_calc,kx(nnx));
            case 'trench'
                K_fem_trench_T1 = f_stiff_fem_cyl('trench',calc_struct,num_elem_T1,Edof_T1,div,NDOF_T1,elem_coord_T1,D_fem_soil_T1,rho_fem_soil_T1,H_T1,radius_T1,omega_calc,kx(nnx));
            end
    
      
            % Transform FEM stiffness matrix (kx,r,n,omega)
            switch calc.system % switch due to different variable names
                case {'tunnel'}
                    % Tunnel T1
                    K_fem_pol_tunnel_T1 = f_trans_stiff_fem(K_fem_tunnel_T1,dis_fem_T1);
                case 'trench'
                    % Trench T1
                    K_fem_pol_trench_T1 = f_trans_stiff_fem(K_fem_trench_T1,dis_fem_T1);
            end
    
            % Calculation of ITM stiffnessmatrix hs_cyl
            [K_itm_hscyl] = f_stiff_itm(calc,dis_itm,geo_T1,mat,trunc_T1,omega_calc,nnx);
                    
            % Assembling of total stiffness matrix 
            switch calc.system
                case 'tunnel'   
                    K_GES_hs_cyl_tun    = f_assem_stiff_hs_cyl(K_fem_pol_tunnel_T1,K_itm_hscyl,dis_itm,dis_fem_T1,calc);
                case 'trench'
                    K_GES_hs_cyl_trench = f_assem_stiff_hs_cyl(K_fem_pol_trench_T1,K_itm_hscyl,dis_itm,dis_fem_T1,calc);
            end
    
            % Assembling of total stiffness load 
            Pz_GES_hs_cyl = f_assem_load_hs_cyl(P_FE_T1,P_FE(nnx,:,nnf),Pz_ITM(nnx,:,nnf),dis_itm,dis_fem_T1,calc);
    
            % Solution of system
            switch calc.system
                case 'tunnel'   
                    [uPz_kxkyf_tot_tun(nnx,:,nnf)] = f_displ_hscyl_tot_solve(K_GES_hs_cyl_tun,Pz_GES_hs_cyl,omega_calc);
                case 'trench'                                                     
                    [uPz_kxkyf_tot_tr(nnx,:,nnf)] = f_displ_hscyl_tot_solve(K_GES_hs_cyl_trench,Pz_GES_hs_cyl,omega_calc);         
            end
            send(D, nnx)
         end % end wavenumber kx loop 
    end % end frequency loop
    close(wait)
    delete(wait)
    clear numCompleted numMax D

%%  Calculate stiffness, assemble, solve for displacements - SERIAL
    case 'parallel'
    for nnf = nf_start:nf_end % for all f of freq without f=0
        omega_nnf = omega(nnf);
        parfor (nnx = nx_start:nx_end,calc_cont.M) 
            omega_calc = omega_nnf;
            % prestore
            K_fem_tunnel_T1     = zeros(NDOF_T1);
            K_fem_trench_T1     = zeros(NDOF_T1);
            K_fem_pol_tunnel_T1 = zeros(NDOF_T1);
            K_fem_pol_trench_T1 = zeros(NDOF_T1);
            K_GES_hs_cyl_tun    = zeros(3*(Ny_hs+Nphi_t_T1+num_nod_in_T1));
            K_GES_hs_cyl_trench = zeros(3*(Ny_hs+Nphi_t_T1+num_nod_in_T1));
                
            if omega_calc/2/pi <0                                   % f<0
                % calc displ for neg f
            elseif omega_calc /2/pi>=-0.1 && omega_calc/2/pi <=+0.1 % f~0
                % As no static solution is included 
                % -> approximate f=0 by solution for very low frequency
                omega_calc = -0.5*2*pi;
            elseif omega_calc > 0.1                                 % f>0
                % calc displ for neg f and take conj complex of result
                omega_calc = -omega_calc;
            end
    
            % Calculate stiffness matrix K_fem for one combination (kx,Omega)
            switch calc.system
            case {'tunnel'} 
                K_fem_tunnel_T1 = f_stiff_fem_cyl('tunnel',calc_struct,num_elem_T1,Edof_T1,div,NDOF_T1,elem_coord_T1,D_fem_soil_T1,rho_fem_soil_T1,H_T1,radius_T1,omega_calc,kx(nnx));
            case 'trench'
                K_fem_trench_T1 = f_stiff_fem_cyl('trench',calc_struct,num_elem_T1,Edof_T1,div,NDOF_T1,elem_coord_T1,D_fem_soil_T1,rho_fem_soil_T1,H_T1,radius_T1,omega_calc,kx(nnx));
            end

            % Transform FEM stiffness matrix (kx,r,n,omega)
            switch calc.system % switch due to different variable names
                case {'tunnel'}
                    % Tunnel T1
                    K_fem_pol_tunnel_T1 = f_trans_stiff_fem(K_fem_tunnel_T1,dis_fem_T1);
                case 'trench'
                    % Trench T1
                    K_fem_pol_trench_T1 = f_trans_stiff_fem(K_fem_trench_T1,dis_fem_T1);
            end
    
            % Calculation of ITM stiffnessmatrix hs_cyl
            [K_itm_hscyl] = f_stiff_itm(calc,dis_itm,geo_T1,mat,trunc_T1,omega_calc,nnx);
                    
            % Assembling of total stiffness matrix 
            switch calc.system
                case 'tunnel'   
                    K_GES_hs_cyl_tun    = f_assem_stiff_hs_cyl(K_fem_pol_tunnel_T1,K_itm_hscyl,dis_itm,dis_fem_T1,calc);
                case 'trench'
                    K_GES_hs_cyl_trench = f_assem_stiff_hs_cyl(K_fem_pol_trench_T1,K_itm_hscyl,dis_itm,dis_fem_T1,calc);
            end
    
            % Assembling of total stiffness load 
            Pz_GES_hs_cyl = f_assem_load_hs_cyl(P_FE_T1,P_FE(nnx,:,nnf),Pz_ITM(nnx,:,nnf),dis_itm,dis_fem_T1,calc);
    
            % Solution of system
            switch calc.system
                case 'tunnel'   
                    [uPz_kxkyf_tot_tun(nnx,:,nnf)] = f_displ_hscyl_tot_solve(K_GES_hs_cyl_tun,Pz_GES_hs_cyl,omega_calc);
                case 'trench'                                                     
                    [uPz_kxkyf_tot_tr(nnx,:,nnf)] = f_displ_hscyl_tot_solve(K_GES_hs_cyl_trench,Pz_GES_hs_cyl,omega_calc);         
            end
            send(D, nnx)
         end % end wavenumber kx loop 
     end % end frequency loop
    close(wait)
    delete(wait)
    clear numCompleted numMax D
end

warning('on','all')

% Gather results
switch calc.system
    case {'tunnel'}
        varname = {'uPx_kxkyf_hscyl_tot','uPy_kxkyf_hscyl_tot','uPz_kxkyf_hscyl_tot'};
        displ.hscyl.(calc.system).tot.(varname{3}) = uPz_kxkyf_tot_tun;
   case {'trench'}
        varname = {'uPx_kxkyf_hscyl_tot','uPy_kxkyf_hscyl_tot','uPz_kxkyf_hscyl_tot'};
        displ.hscyl.(calc.system).tot.(varname{3}) = uPz_kxkyf_tot_tr;
end


end % end function

