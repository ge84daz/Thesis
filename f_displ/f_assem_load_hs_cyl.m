%__________________________________________________________________________
%
%        f_assem_load_hs_cyl.m    
%
%       INPUT: 
%       - P_FE          = FEM load vector for one combination (kx,Omega)
%       - P_ITM         = ITM load vector for one combination (kx,Omega)
%       - Ny_hs         = Number of Fourier members on ground surface in y
%       - Nphi_t        = Number of Fourier members on tunnel surface
%       - num_nod_in    = Number of nodes inside Tunnel
%
%       OUTPUT: 
%       - P_GES         = Total load vector of the complete ITM-FEM System
%                         for one combination (kx,Omega)
%
%       DESCRIPTION: 
%       - Assembles the total load vecotr of the ITM-FEM hs cyl 
%       - tunnel -> stiffness matrix K_GES has same size
%   
%       REMARK: 
%       Original Author: Julian Freisinger (f_assem_load_hs_cyl)
%       Modified:   Tom Hicks                 
%       Date:       25-09-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________


function [P_GES] = f_assem_load_hs_cyl(P_FE_T1,P_FE,P_ITM,dis_itm,dis_fem_T1,calc)

%  Initialize: Input
Ny_hs          = dis_itm.Ny_hs;    % number of points at hs surf
Nphi_t_T1      = dis_itm.Nphi_t;   % number of points at cyl bound

% Tunnel 1
num_nod_in_T1  = dis_fem_T1.num_nod_in; % number of points interior FEM domain

% Assemble
switch calc.system
    case {'tunnel','trench'}
        P_GES = zeros(3*(Ny_hs+Nphi_t_T1+num_nod_in_T1),1);
        % Load at halfspace surface z=0
        P_GES(1:3*Ny_hs) = P_ITM;
        
        % Load within FEM Tunnel T1
        if strcmpi(P_FE_T1,'yes')
        P_GES(3*Ny_hs+1:3*(Ny_hs+Nphi_t_T1+num_nod_in_T1)) = P_FE;
        end  
end
end