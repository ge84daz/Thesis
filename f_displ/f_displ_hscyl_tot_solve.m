%__________________________________________________________________________
%
%        f_displ_hs_cyl_tot.m    
%
%       INPUT: 
%       - K_GES = total dynamic stiffness matrix
%       - P_GES = total load vector of calculated system
%       - omega = omega discretization
%
%       OUTPUT: 
%       - u_GES = total displacements due to given load case
%
%       DESCRIPTION: 
%       - u_GES = K_GES\P_GES; = u_GES*K_GES = P_GES - solve for u_GES
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_displ_hs_cyl_tot)
%       Modified:   Tom Hicks                 
%       Date:       25-09-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________

function [u_GES] = f_displ_hscyl_tot_solve(K_GES,P_GES,omega_calc)

if     omega_calc < 0;  u_GES = K_GES\P_GES;       % negative f and f=0 -> f=-0.5*2*pi
elseif omega_calc > 0;  u_GES = conj(K_GES\P_GES); % postive  f
end