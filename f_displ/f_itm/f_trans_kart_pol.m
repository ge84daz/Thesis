%__________________________________________________________________________
%
%        f_trans_kart_pol.m    
%
%       INPUT: 
%       - stress/displacements at one point of the cylinderical surface (stress_kart)
%       - angle phi of the corresponding point (phi)
%
%       OUTPUT: 
%       - stress_pol  = transformed stress/displ from cartesian to polar coordinates 
%
%       DESCRIPTION: 
%       |u_x        | 
%       |u_y        |       |u_x       |
%       |u_z        |       |u_r       | 
%       |sigma_xx   |       |u_phi     | 
%       |sigma_yy   | ----> |sigma_r   | 
%       |sigma_zz   |       |sigma_rphi| 
%       |sigma_xy   |       |sigma_rx  | 
%       |sigma_zy   |            
%       |sigma_zx   | 
%          [9x1]               [6x1]
%
%       REMARK: 
%       Original Author: Georg Frühe (trans_kart_pol)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [stress_pol] = f_trans_kart_pol(stress_kart,phi)

% transformation matrix [6x9] - Diss Frühe p.22
T = [1   0          0          0   0                   0                   0           0                       0
     0  -sin(phi)   cos(phi)   0   0                   0                   0           0                       0
     0  -cos(phi)  -sin(phi)   0   0                   0                   0           0                       0
     0   0          0          0   sin(phi)^2          cos(phi)^2          0          -2*sin(phi)*cos(phi)     0
     0   0          0          0   sin(phi)*cos(phi)  -sin(phi)*cos(phi)   0           sin(phi)^2-cos(phi)^2   0
     0   0          0          0   0                   0                  -sin(phi)    0                       cos(phi)];


% transformation of stress and displacements
stress_pol = T * stress_kart;

















