%__________________________________________________________________________
%
%        f_trans_pol_kart.m    
%
%       INPUT: 
%       - stress/displacements at one point of the halfspace (stress_pol)
%       - angle phi of the corresponding point (phi)
%
%       OUTPUT: 
%       - stress_pol  = transformed stress/displ from polar to cartesian coordinates 
%
%       DESCRIPTION: 
%       |u_x          | 
%       |u_r          |       |u_x       |
%       |u_phi        |       |u_y       | 
%       |sigma_xx     |       |u_z       | 
%       |sigma_rr     | ----> |sigma_zz  | 
%       |sigma_phiphi |       |sigma_zy  | 
%       |sigma_rx     |       |sigma_zx  | 
%       |sigma_rphi   |            
%       |sigma_phix   | 
%          [9x1]               [6x1]
%
%       REMARK: 
%       Original Author: Georg Frühe (trans_pol_kart)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [stress_kart] = f_trans_pol_kart(stress_pol,phi)

% transformation matrix [6x9] - Diss frühe p.22
T=[ 1   0          0          0   0                   0                   0          0                       0
    0  -sin(phi)  -cos(phi)   0   0                   0                   0          0                       0
    0   cos(phi)  -sin(phi)   0   0                   0                   0          0                       0
    0   0          0          0   cos(phi)^2          sin(phi)^2          0         -2*sin(phi)*cos(phi)     0
    0   0          0          0  -sin(phi)*cos(phi)   sin(phi)*cos(phi)   0          sin(phi)^2-cos(phi)^2   0
    0   0          0          0   0                   0                   cos(phi)   0                      -sin(phi)];

% transformation of stress and displacements
stress_kart = T * stress_pol;





