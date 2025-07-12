%__________________________________________________________________________
%
%        f_stiff_itm.m    
%
%       INPUT: 
%       - dis_itm     = Discretization ITM
%       - geo         = Geometry
%       - mat         = Material properties
%       - trunc       = trunction criteria
%       - nomega,nx   = index of current freq omega and wavenumber kx
%
%       OUTPUT: 
%       - K_itm       = stiffness matrix of ITM system hs with cyl. cavity
%
%       DESCRIPTION: 
%       K_itm =
%       |          |        |      ky(1)   ...    ky(n)     |     n_t(1)    ...    n_t(n)     |
%       |          |--------|---------------|---------------|----------------|----------------|
%       |          |        | Pzx  Pzy  Pzz | Pzx  Pzy  Pzz |  Prr Prphi Prx |  Prr Prphi Prx |
%       |----------|--------|---------------|---------------|----------------|----------------|
%       |   ky (1) |u_x     | ...  ...  ...   ...  ...  ... |  ...  ...  ...    ...  ...  ... | 
%       |          |u_y     | ...  ...  ...   ...  ...  ... |  ...  ...  ...    ...  ...  ... | 
%       |          |u_z     | ...  ...  ...   ...  ...  ... |  ...  ...  ...    ...  ...  ... | 
%       |     :    |   :    |       :               :       |        :                :       |
%       |   ky(n)  |   :    | ...  ...  ...   ...  ...  ... |  ...  ...  ...    ...  ...  ... | 
%       |----------|--------|-------------------------------|---------------- ----------------|
%       |  n_t(1)  |u_x     | ...  ...  ...   ...  ...  ... |  ...  ...  ...    ...  ...  ... | 
%       |          |u_r     | ...  ...  ...   ...  ...  ... |  ...  ...  ...    ...  ...  ... | 
%       |          |u_phi   | ...  ...  ...   ...  ...  ... |  ...  ...  ...    ...  ...  ... | 
%       |     :    |   :    |       :               :       |        :                :       |
%       |  n_t(n)  |   :    | ...  ...  ...   ...  ...  ... |  ...  ...  ...    ...  ...  ... | 
%       |----------|--------|-------------------------------|---------------- ----------------|
%
%
%       REMARK:
%       1.) Truncation criteria for Fourier terms
%       2.) Stress/displ at CAVITY  - unit loads at HS SURF
%       3.) Stress/displ at HS SURF - unit loads at CAVITY
%       4.) Stress/displ at HS SURF - unit loads at HS SURF
%       5.) Stress/displ at CAVITY  - unit loads at CAVITY
%       6.) Flexibility matrix of coupled ITM-FEM system
%       7.) Resorting and stiffness matrix
%
%       REMARKS:
%       Compared to original function of Frühe 2010 last_verformung_ofl_tun
%       - Without expanding the matrices with splines in case of trench
%       - without additional point at By/2+1 and corr. averaging of outmost points
%       - Without inclusion of influence of repetition sections
%       - Corresponds to the function f_last_verf_ofl_tun_wo_spline_wo_wh
%
%       Direct calculation of stiffness matrix leads to more instable results of
%       final displacements uz due to bad condition of the stress and
%       displ matrices. The calculation of the stiffness K_itm via the
%       flexibility leads to more stable results for uz, as somehow the inversion
%       of the better conditioned N_ITM leads to a more exact K_itm.
%       
%       REMARK: 
%       Original Author: Julian Freisinger (f_stiff_itm)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________



function [K_ITM] = f_stiff_itm(~,dis_itm,geo_T1,mat,trunc_T1,omega,nx)


%% 0.) Initialize: Input

% Paramters equal for T1 
ky          = dis_itm.ky;               % total discretization of hs surf
ky0         = dis_itm.ky0;              % increment of ky
Ny_hs       = dis_itm.Ny_hs;            % total amount of points at hs surf
Nphi_t_T1   = dis_itm.Nphi_t;           % total amount of points at cyl bound

kx          = dis_itm.kx(nx);           % step of kx

% get material parameters of halfspace soil
cs2         = mat.hs.cs;
cp2         = mat.hs.cp;
mue2        = mat.hs.mu;
lambda2     = mat.hs.lambda;

% nodes at hs surf
% OF = Oberflache
N_OF        = geo_T1.N_OF;              % number of points at hs surf
geo_tun_tr  = geo_T1.geo;
r0_T1       = geo_T1.radius;            % radius of tunnel T1
r_OF_T1     = geo_T1.r_OF;              % r-coord. of points at hs surf
phi_OF_T1   = geo_T1.phi_OF;            % phi-coord. of points at hs surf
y_OF_ges_T1 = geo_T1.y_OF_ges;          % y-coord. of points at hs surf

% nodes at cyl bound. T1
% LR = Lochrand
N_LR_T1      = geo_T1.N_LR;             % number of points at cyl bound
phi_LR_T1    = geo_T1.phi_LR;           % phi-coord. of points at cyl bound
y_LR_T1      = geo_T1.y_LR;             % y-coord. of points at cyl bound
z_LR_T1      = geo_T1.z_LR;             % z-coord. of points at cyl bound

y_Tc_T1      = geo_T1.y_Tc;             % y-coord. of tunnel center

% truncation criteria
delta_wh_T1    = trunc_T1.delta_wh;     % truncation criteria for repetition section
abstand_wh_T1  = trunc_T1.abstand_wh;   % distance from tunnel center to edge of next section
delta_T1       = trunc_T1.delta;        % truncation criteria for cyl bound - hs surf
abstand_T1     = trunc_T1.abstand;      % distance from tunnel upper boundary to hs surf (Ueberdeckung)

% inform matlab to use bessel function
coder.extrinsic('besselh');

%% 1.) Truncation criteria for Fourier terms
% Determine number of Fourier terms that have an influence at surface/tunnel
% compare Diss. Frühe 2010 p.78

% Predefinitions
ks          = omega/cs2;
kbeta       = sqrt(ks^2-kx^2);
kbeta_r0_T1 = kbeta*r0_T1;

%__________________________________________________________________________
% Influence of load tunnel to REPETITION SECTIONS
%__________________________________________________________________________
% Number of Fourier terms at cavity boundary that have an influence on the
% repetition section in y-direction
% For tunnel T1
kbeta_r_wh_T1  = kbeta*abstand_wh_T1;
N_t_wh_T1      = 0; % preset order of the Hankel function = 0

% using Hankel functions of first kind with order n = N_t_wh_T1 
% check till what order (kx,n) has to be taken into account
while abs(besselh(N_t_wh_T1,kbeta_r_wh_T1))/abs(besselh(N_t_wh_T1,kbeta_r0_T1))>delta_wh_T1 && N_t_wh_T1<Nphi_t_T1/2
    N_t_wh_T1 = N_t_wh_T1+1;
end

%__________________________________________________________________________
% Inflence of load tunnel|surface to surface|tunnel displ
%__________________________________________________________________________

% Influnce of load tunnel|surface to surface|tunnel tractions/displ
% -----------------------------------------------------------------
switch geo_tun_tr
    case 'tunnel'
        
        %% Number of Fourier members at cavity that have an influence at the hs surface
        % For tunnel T1
        kbeta_r_T1 = kbeta*(r0_T1+abstand_T1);
        N_t_T1     = N_t_wh_T1;
        
        while abs(besselh(N_t_T1,kbeta_r_T1))/abs(besselh(N_t_T1,kbeta_r0_T1))>delta_T1 && N_t_T1<Nphi_t_T1/2
            N_t_T1 = N_t_T1+1;
        end   

        %% Number of Fourier members at hs surf that have an influence at the cavity
        %  e^-(lambda2)    -> amplitude of wave in depth direction attenuated until cavity boundary is reached
        %  lambda2 > 1.3kR -> no influence of Rayleigh wave onto the cavity
        
        % For tunnel T1
        N_h_T1 = 0;
        ky1_T1 = 0;
        kR_T1  = 1.2*abs(ks); % approximation for Rayleigh wavenumber
        while (abs(exp(-sqrt(kx^2+ky1_T1^2-ks^2)*abstand_T1))>delta_T1 || sqrt(kx^2+ky1_T1^2)<1.3*kR_T1) && N_h_T1<Ny_hs/2
            %    while abs(exp(-sqrt(kx^2+ky1^2-ks^2)*abstand))>delta && N_h<Nmax_h/2
            N_h_T1 = N_h_T1+1;
            ky1_T1 = N_h_T1*ky0;
        end
        
    otherwise % trench,slit
        N_t_T1 = Nphi_t_T1/2;
        N_h_T1 = Ny_hs/2;
end

%% 2.) Stress/displ at CAVITY - unit loads at HS SURF
%_________________________________________________________________________________________________
% Determination of stress/displ at the CAVITY (y,z,kx,omega) due to unit load cases at hs surf
%_________________________________________________________________________________________________

% Tunnel: For all nodes on cavity boundary
% Tunnel T1
% LR = Lochrand - cylindrical boundary
% kart = cartesian coordinates
stress_LR_kart_T1 = complex(zeros(9,N_LR_T1*3*2*N_h_T1));
for s = 1:2*N_h_T1  % for all relevant series members N_h in y-direction at hs surf
    n_h = (Ny_hs/2-N_h_T1)+s;
    stress_LR_kart_T1(:,N_LR_T1*3*(s-1)+1:N_LR_T1*3*s) = f_stress_cyl_load_surf(kx,ky(n_h),omega,y_LR_T1,z_LR_T1,mue2,lambda2,cp2,cs2);
end
% stress_LR_kart_T1 =
%            |     N_h(1)    |     N_h(2)    |     N_h(3)     |  ...  
% -------------------------------------------------------------------
%            | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx  |  ...    
% -------------------------------------------------------------------
% u_x        | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    
% u_y        | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    
% u_z        | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    
% sigma_xx   | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    
% sigma_yy   | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    
% sigma_zz   | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    
% sigma_xy   | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    
% sigma_zy   | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    
% sigma_zx   | ...  ...  ...   ...  ...  ...    ...  ...  ... |  ...    

%________________________________________________________________________________________
% Transform stress/displ at CAVITY form cartesian into polar coordinates (phi,r,kx,omega)
%________________________________________________________________________________________

% Tunnel: For all nodes on cavity boundary (one after one)
% Tunnel T1
stress_LR_pol_temp_T1 = complex(zeros(6*Nphi_t_T1,3*2*N_h_T1));
for n_t=1:N_LR_T1  % for all r,phi on cavity boundary
    stress_LR_pol_temp_T1(6*(n_t-1)+1:6*n_t,:) = f_trans_kart_pol(stress_LR_kart_T1(:,n_t:N_LR_T1:N_LR_T1*3*2*N_h_T1),phi_LR_T1(n_t));
end
clear stress_LR_kart_T1 
% stress_LR_pol_T1 =
%           |           |     N_h(1)    |     N_h(2)    |     N_h(3)    |  ...  
%           |-----------|---------------|---------------|---------------|-----
%           |           | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx |  ...  
%|----------|-----------|---------------|---------------|---------------|-----
%|  phi(1)  |u_x        | ...  ...  ... | ...  ...  ... | ...  ...  ... |  ...  
%|          |u_r        | ...  ...  ... | ...  ...  ... | ...  ...  ... |  ...  
%|          |u_phi      | ...  ...  ... | ...  ...  ... | ...  ...  ... |  ...  
%|          |sigma_r    | ...  ...  ... | ...  ...  ... | ...  ...  ... |  ...  
%|          |sigma_rphi | ...  ...  ... | ...  ...  ... | ...  ...  ... |  ...  
%|          |sigma_rx   | ...  ...  ... | ...  ...  ... | ...  ...  ... |  ...  
%|----------|-----------|---------------|---------------|---------------|-----
%|  phi(2)  |u_x        | ...  ...  ... | ...  ...  ... | ...  ...  ... |  ...  
%|    :     |     :     |       :       |       :       |       :


%______________________________________________________________________________
% Sorting of stress/displ at CAVITY into matrix stress_LR_pol (phi,r,kx,omega)
%______________________________________________________________________________
if strcmpi(geo_tun_tr,'trench')
    % Sorting of stress/displ at CAVITY into matrix stress_LR_pol.
    % Top/bottom: nodes of trench within hs
    % Middle:     nodes above hs surf -> set to zero
    %
    % XXXXXXXXX y_OF [0 ... -radius]
    % 000000000 y_OF [-radius ... 0]
    % 000000000 y_OF [0 ... +radius]
    % XXXXXXXXX y_OF [+radius ... 0]
    %
    % Trench T1
    stress_LR_pol_T1 = complex(zeros(6*Nphi_t_T1,3*2*N_h_T1));
    stress_LR_pol_T1(1:6*(N_LR_T1)/2,:) = stress_LR_pol_temp_T1(1:6*(N_LR_T1)/2,:);
    stress_LR_pol_T1(6*(N_LR_T1)/2+6*(Nphi_t_T1-N_LR_T1+1)+1:6*Nphi_t_T1,:) = stress_LR_pol_temp_T1(6*(N_LR_T1)/2+1:6*(N_LR_T1-1),:);

else % 'tunnel'
    % Tunnel 1
    stress_LR_pol_T1 = stress_LR_pol_temp_T1;
end


%__________________________________________________________________________
% Develop stress/displ at CAVITY into Fourier series along boundary (n_t,r,kx,omega)
%__________________________________________________________________________
% Tunnel T1
stress_LR_T1_POF = complex(zeros(6*Nphi_t_T1,3*Ny_hs));
for s=1:6  % for all stress/displ
    stress_LR_T1_POF(s:6:6*Nphi_t_T1,3*(Ny_hs/2-N_h_T1)+1:3*(Ny_hs/2+N_h_T1)) = fftshift(1/Nphi_t_T1*fft(stress_LR_pol_T1(s:6:6*Nphi_t_T1,:)),1);
end
clear stress_LR_pol_T1
% stress_LR_T1_POF =
%           |           |     ky(1)     |     ky(2)     |     ky(3)      |  ...  
%           |-----------|---------------|---------------|----------------|-----
%           |           | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx  |  ...  
%|----------|-----------|---------------|---------------|----------------|-----
%|  n_t(1)  |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%|          |u_r        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%|          |u_phi      | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%|          |sigma_r    | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%|          |sigma_rphi | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%|          |sigma_rx   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%|----------|-----------|---------------|---------------|----------------|-----
%|  n_t(2)  |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%|    :     |     :     |       :       |       :       |        :

% The ROWS of the matrix contain:
% stress/displ (u_x, u_r, u_phi, sigma_r, sigma_rphi, sigma_rx) one after another
% for all Fourier members n_t (from -Nmax_t/2 to Nmax_t/2-1) at cavity
% boundary due to one unit load case at hs surf
%
% Over the COLUMNS the matrix contains the stress/displ for the
% load cases (Pzz,Pzy,Pzx) consecutively for all wavenumbers ky
% (-kymax ... +kymax) in y-dir at the hs surf

%% 3.) Stress/displ at HS SURF - unit loads at CAVITY
%__________________________________________________________________________
% Determination of stress/displ at the HS SURF (phi,r,kx,omega) due to unit
% load cases at the CAVITY
%__________________________________________________________________________

% Tunnel: For all nodes on hs surf
% Tunnel T1
% OF = surface of hs (N_OF = points at surface)
stress_OF_pol_T1 = complex(zeros(9,N_OF*3*2*N_t_T1));
for s=1:2*N_t_T1  % for all relevant Fourier members n_t at cavity
    nn = -N_t_T1+s-1;
    stress_OF_pol_T1(:,N_OF*3*(s-1)+1:N_OF*3*s) = f_stress_surf_load_cyl(nn,kx,omega,r_OF_T1,phi_OF_T1,mue2, cp2, cs2, r0_T1);
end
% stress_OF_pol_T1 =
%               |      n_t(1)     |      n_t(2)     |      n_t(3)      |  ...  |
% --------------|-----------------|-----------------|------------------|-------|
%               | Prr  Prphi  Prx | Prr  Prphi  Prx | Prr  Prphi  Prx  |  ...  |
% --------------|-----------------|-----------------|------------------|-------|
% u_x           | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
% u_r           | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
% u_phi         | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
% sigma_xx      | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
% sigma_rr      | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
% sigma_phiphi  | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
% sigma_rx      | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
% sigma_rphi    | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
% sigma_phix    | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |

%__________________________________________________________________________
% Transform stress/displ at HS SURF into cartesian coord (y,z,kx,omega)
%__________________________________________________________________________

% Tunnel: For all nodes on hs surf
% Tunnel T1
stress_OF_kart2_T1 = complex(zeros(6*(Ny_hs),3*2*N_t_T1));
for n_h=1:N_OF  % for all loactions y on hs surf
    stress_OF_kart2_T1(6*(n_h-1)+1:6*n_h,:)=f_trans_pol_kart(stress_OF_pol_T1(:,n_h:N_OF:N_OF*3*2*N_t_T1),phi_OF_T1(n_h));
end
clear stress_OF_pol_T1
% stress_OF_kart_T1 =
%           |           |      n_t(1)     |      n_t(2)     |      n_t(3)     |  ...  
%           |-----------|-----------------|-----------------|-----------------|
%           |           | Prr  Prphi  Prx | Prr  Prphi  Prx | Prr  Prphi  Prx |  ...  
%|----------|-----------|-----------------|-----------------|-----------------|
%|  y(1)    |u_x        | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |u_y        | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |u_z        | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |sigma_zz   | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |sigma_zy   | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |sigma_zx   | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|----------|-----------|-----------------|-----------------|-----------------|-----
%|  y(2)    |u_x        | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|    :     |     :     |        :        |        :        |        :



%__________________________________________________________________________
% Sorting of stress/displ at HS SURF into matrix stress_OF_kart (y,z,kx,omega)
%__________________________________________________________________________
if strcmp(geo_tun_tr,'trench')
    % Sorting of stress/displ at HS SURF into matrix stress_OF_kart
    % including all location at hs surf y_OF = [0...By/2-dy, -By/2...-dy]
    %
    % stresses & displ should be zero within trench
    % -> sort them at rows of matrix, where y_OF has values that are
    %    located outside of the trench. As in total 6 displ & stresses, blow
    %    up y_OF_ges by 6 and then search for values outside trench
    % -> this gives directly the rows for sorting the displ/stress into the
    %    matrix stress_OF_kart_T1 or/and stress_OF_kart_T2
    % -> For example y_Tc =0
    %    000000000 y_OF [0 ... +radius-dy]
    %    XXXXXXXXX y_OF [+radius ... By/2-dy]
    %    XXXXXXXXX y_OF [-By/2 ... -radius]
    %    000000000 y_OF [-radius+dy ... -dy]
    
    % Trench T1
    for n=1:1:6;temp1(n:6:6*Ny_hs) = y_OF_ges_T1; end
    sort_T1 = temp1<=+y_Tc_T1-r0_T1 |temp1 >=+y_Tc_T1+r0_T1;
    
    stress_OF_kart_T1 = complex(zeros(6*Ny_hs,3*2*N_t_T1));
    stress_OF_kart_T1(sort_T1,:) = stress_OF_kart2_T1(1:6*N_OF,:);

else % 'tunnel'   
    % Tunnel 1
    stress_OF_kart_T1 = stress_OF_kart2_T1;
end

%__________________________________________________________________________________
% Develope stress/displ at HS SURF into Fourier series along hs surf (ky,z,kx,omega)
%__________________________________________________________________________________

% Tunnel T1
stress_OF_PLR_T1 = complex(zeros(6*Ny_hs,3*Nphi_t_T1));
for s=1:6  % for all stress/displ
    stress_OF_PLR_T1(s:6:6*Ny_hs,3*(Nphi_t_T1/2-N_t_T1)+1:3*(Nphi_t_T1/2+N_t_T1)) = fftshift(1/Ny_hs*fft(stress_OF_kart_T1(s:6:6*Ny_hs,:)),1);
end
clear stress_OF_kart_T1
% stress_OF_PLR_T1 =
%           |           |      n_t(1)     |      n_t(2)     |      n_t(3)     |  ...  
%           |-----------|-----------------|-----------------|-----------------|
%           |           | Prr  Prphi  Prx | Prr  Prphi  Prx | Prr  Prphi  Prx |  ...  
%|----------|-----------|-----------------|-----------------|-----------------|
%|  ky(1)   |u_x        | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |u_y        | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |u_z        | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |sigma_zz   | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |sigma_zy   | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|          |sigma_zx   | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|----------|-----------|-----------------|-----------------|-----------------|-----
%|  ky(2)   |u_x        | ...   ...   ... | ...   ...   ... | ...   ...   ... |  ...  
%|    :     |     :     |        :        |        :        |        :

% The ROWS of the matrix contain:
% stress/displ (u_x, u_y, u_z, sigma_zz, sigma_zy, sigma_zx) one after another
% for all wavenumbers ky (from -Ny_hs/2 to +Ny_hs/2-1) at hs surface due
% to one unit load case at cyl boundary
%
% Over the COLUMNS the matrix contains the stress/displ for the
% load cases (Prr, Prphi, Prx) consecutively for all Fourier members 
% n_t(from -Nmax_t/2 to Nmax_t/2-1) at cyl boundary


%% 4.) Stress/displ at HS SURF - unit loads at HS SURF
%__________________________________________________________________________
% Determination of stress/displ at the HS SURF (ky,z,kx,omega) due to unit
% load cases at the HS SURF
%__________________________________________________________________________

% original formulation
stress_OF_POF    = complex(zeros(6*Ny_hs,3*Ny_hs));
for n_h = 1:Ny_hs   % for all wavenumbers ky_ev in y-dir at surface
    stress_OF_POF(6*(n_h-1)+1:6*n_h,3*(n_h-1)+1:3*n_h) = f_stress_surf_load_surf(kx,ky(n_h),omega,mue2, lambda2, cp2, cs2);
end

% stress_OF_POF =
%           |           |     ky(1)     |     ky(2)     |     ky(3)      |  ...  |
%           |-----------|---------------|---------------|----------------|-------|
%           |           | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx  |  ...  |
%|----------|-----------|---------------|---------------|----------------|-------|
%|   ky(1)  |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  |
%|          |u_y        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  |
%|          |u_z        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  |
%|          |sigma_zz   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  |
%|          |sigma_zy   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  |
%|          |sigma_zx   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  |
%|----------|-----------|---------------|---------------|----------------|-------|
%|  ky(2)   |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  |

% The ROWS of the matrix contain:
% stress/displ (u_x, u_y, u_z, sigma_zz, sigma_zy, sigma_zx) one after another
% for all wavenumbers ky (from -Ny_hs/2 to +Ny_hs/2-1) at hs surface due
% to one unit load case at hs surf
%
% Over the COLUMNS the matrix contains the stress/displ for the
% load cases (Pzz,Pzy,Pzx) consecutively for all wavenumbers ky
% (-kymax ... +kymax) in y-dir at the hs surf

%% 5.) Stress/displ at CAVITY  - unit loads at CAVITY
%__________________________________________________________________________
% Determination of stress/displ at the CAVITY (n_t,r,kx,omega) due to unit
% load cases at the CAVITY
%__________________________________________________________________________
% Tunnel T1
stress_LR_PLR_T1  =complex(zeros(6*Nphi_t_T1,3*Nphi_t_T1));
for s=1:Nphi_t_T1  % for all Fourier members along cavity boundary
    n_t = -Nphi_t_T1/2+s-1;
    stress_LR_PLR_T1(6*(s-1)+1:6*s,3*(s-1)+1:3*s) = f_stress_cyl_load_cyl(n_t,kx,omega,mue2, cp2, cs2,r0_T1);
end
% stress_LR_PLR_T1 =
% |          |           |     n_t(1)    |     n_t(2)    |     n_t(3)     |  ...  
% |          |-----------|---------------|---------------|----------------|-------
% |          |           | Prr Prphi Prx | Prr Prphi Prx | Prr Prphi Prx  |  ...  
% |----------|-----------|---------------|---------------|----------------|-------
% |   n_t(1) |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
% |          |u_r        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
% |          |u_phi      | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
% |          |sigma_rr   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
% |          |sigma_rphi | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
% |          |sigma_rx   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
% |----------|-----------|---------------|---------------|----------------|-------
% |   n_t(2) |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%     :           :              :               :                :           : 
% The ROWS of the matrix contain:
% stress/displ (u_x, u_r, u_phi, sigma_rr, sigma_rphi, sigma_rx) one after another
% for all Fourier members n_t(from -Nmax_t/2 to Nmax_t/2-1) at cyl bound
% due to one unit load case at cyl bound
%
% Over the COLUMNS the matrix contains the stress/displ for the
% load cases (Prr, Prphi, Prx) consecutively for all Fourier members 
% n_t(from -Nmax_t/2 to Nmax_t/2-1) at cyl boundary

%% 7.) Flexibility matrix of coupled ITM-FEM system

% Assembling of total coupling matrix (stresses and displ)
coupling = [stress_OF_POF     stress_OF_PLR_T1
            stress_LR_T1_POF  stress_LR_PLR_T1];
clear stress_OF_POF stress_OF_PLR_T1 stress_LR_T1_POF stress_LR_PLR_T1 

% Adapt sorting for 2nd tunnel via Nphi_t_T2
Nphi_t_T2 = 0; % as if no 2nd tunnel exists

% Seperate into stress and displ matrices
stress     = complex(zeros(3*((Ny_hs)+(Nphi_t_T1 + Nphi_t_T2))));
displ      = complex(zeros(3*((Ny_hs)+(Nphi_t_T1 + Nphi_t_T2))));

for s=1:3
    stress(    s:3:3*(Ny_hs+Nphi_t_T1+Nphi_t_T2),:) = coupling(s+3:6:6*(Ny_hs+Nphi_t_T1+Nphi_t_T2),:);
    displ(s:3:3*(Ny_hs+Nphi_t_T1+Nphi_t_T2),:) = coupling(s:6:6*  (Ny_hs+Nphi_t_T1+Nphi_t_T2),:);
end
% clear coupling
% Flexibility Matrix for all unit load cases
N_itm_us = -displ/stress;
clear displ stress


%% 8.) Resorting and stiffness matrix
% Remark: Structure of Matrix:
% Frühe:    Pzz,Pzy,Pzx         Ground surface
%           Prr, Prphi, Prx     Tunnel surface
%
% Needed:   Pzx,Pzy,Pzz
%           Prx,Prr,Prphi
% -> resortation needed

N_ITM = complex(eye(3*(Ny_hs+Nphi_t_T1+Nphi_t_T2)));

N_ITM(:,1:3:3*Ny_hs) = N_itm_us(:,3:3:3*Ny_hs);
N_ITM(:,2:3:3*Ny_hs) = N_itm_us(:,2:3:3*Ny_hs);
N_ITM(:,3:3:3*Ny_hs) = N_itm_us(:,1:3:3*Ny_hs);

N_ITM(:,3*Ny_hs+1:3:3*(Ny_hs+Nphi_t_T1+Nphi_t_T2)) = N_itm_us(:,3*Ny_hs+3:3:3*(Ny_hs+Nphi_t_T1+Nphi_t_T2));
N_ITM(:,3*Ny_hs+2:3:3*(Ny_hs+Nphi_t_T1+Nphi_t_T2)) = N_itm_us(:,3*Ny_hs+1:3:3*(Ny_hs+Nphi_t_T1+Nphi_t_T2));
N_ITM(:,3*Ny_hs+3:3:3*(Ny_hs+Nphi_t_T1+Nphi_t_T2)) = N_itm_us(:,3*Ny_hs+2:3:3*(Ny_hs+Nphi_t_T1+Nphi_t_T2));
clear N_itm_us

K_ITM = eye(3*(Ny_hs+Nphi_t_T1+Nphi_t_T2))/N_ITM;

end % end function






