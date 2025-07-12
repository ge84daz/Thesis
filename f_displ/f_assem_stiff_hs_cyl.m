%__________________________________________________________________________
%
%        f_assem_stiff_hs_cyl.m    
%
%       INPUT: 
%       - Each stiffness matrix is given for one combination (kx,omega)
%       - K_fem_pol_T1 = Stiffness of FEM system tunnel T1
%       - K_itm_hs_cyl = Stiffness of ITM system hs_cyl with one tunnel
%
%       OUTPUT: 
%       - K_GES = Total stiffness of complete ITM-FEM system for one (kx,Omega)
%
%       DESCRIPTION: 
%   	- Assembles the total stiffness matrix of the complete ITM-FEM System.
%       - tunnel, trench, slit  -> stiffness matrix K_GES has same size
%       - hs_cyl:
%       ______________________________________________________________________
%       K_itm_hs_cyl_LL    K_itm_hs_cyl_LG1                     0
%       K_itm_hs_cyl_G1L   K_itm_hs_cyl_G1G1 + K_fem_pol_G1G1   K_fem_pol_G1O1
%       0                  K_fem_pol_O1G1                       K_fem_pol_O1O1
%       ______________________________________________________________________
%       with:
%               L = Lambda  - halfspace surface
%               G = Gamma1  - cylindrical boundary of tunnel T1
%               O = Omega   - interior fem domain
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_assem_stiff_hs_cyl)
%       Modified:   Tom Hicks                 
%       Date:       25-09-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________

function [K_GES] = f_assem_stiff_hs_cyl(varargin)

%% Initialize: Input

% Input arguments
K_fem_pol_T1    = varargin{1};
K_itm_hs_cyl    = varargin{2};
dis_itm         = varargin{3};
dis_fem_T1      = varargin{4};
calc            = varargin{5};

% Variable declaration
% Tunnel T1
Ny_hs         = dis_itm.Ny_hs;          % number of points at hs surf
Nphi_t_T1     = dis_itm.Nphi_t;         % number of points at cyl bound
num_nod_in_T1 = dis_fem_T1.num_nod_in;  % number of nodes interior FEM domain

%% Stiffness members (all)
%   K_GES = 
%   ______________________________________________________________________
%   K_itm_hs_cyl_LL    K_itm_hs_cyl_LG1                     0
%   K_itm_hs_cyl_G1L   K_itm_hs_cyl_G1G1 + K_fem_pol_G1G1   K_fem_pol_G1O1
%   0                  K_fem_pol_O1G1                       K_fem_pol_O1O1
%   ______________________________________________________________________
K11 = K_itm_hs_cyl(1:3*Ny_hs,1:3*Ny_hs);                        
K12 = K_itm_hs_cyl(1:3*Ny_hs,3*Ny_hs+1:3*(Ny_hs+Nphi_t_T1));   
K13 = zeros(3*Ny_hs,3*num_nod_in_T1);

K21 = K_itm_hs_cyl(3*Ny_hs+1:3*(Ny_hs+Nphi_t_T1),1:3*Ny_hs);
K22 = K_itm_hs_cyl(3*Ny_hs+1:3*(Ny_hs+Nphi_t_T1),3*Ny_hs+1:3*(Ny_hs+Nphi_t_T1)) + K_fem_pol_T1(1:3*Nphi_t_T1,1:3*Nphi_t_T1);
K23 = K_fem_pol_T1(1:3*Nphi_t_T1,3*Nphi_t_T1+1:3*(Nphi_t_T1+num_nod_in_T1));

K31 = zeros(3*num_nod_in_T1,3*Ny_hs);
K32 = K_fem_pol_T1(3*Nphi_t_T1+1:3*(Nphi_t_T1+num_nod_in_T1),1:3*Nphi_t_T1);
K33 = K_fem_pol_T1(3*Nphi_t_T1+1:3*(Nphi_t_T1+num_nod_in_T1),3*Nphi_t_T1+1:3*(Nphi_t_T1+num_nod_in_T1));

%% Assemble K_Ges
K_GES = [ K11  K12  K13
          K21  K22  K23
          K31  K32  K33 ];
   


end
