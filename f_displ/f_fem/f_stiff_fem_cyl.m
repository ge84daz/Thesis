%__________________________________________________________________________
%
%        f_stiff_fem_cyl.m    
%
%       INPUT: 
%       - num_elem                      = Number of elements
%       - elem_coord                    = Elements with coordinates
%       - Edof                          = Mapping Matrix
%       - H                             = Depth of tunnel center on z-axis
%       - omega(nomega)                 = Frequency
%       - kx(nx)                        = Inkrement in wavenumber domain in x
%       - NDOF                          = Number of DOFs
%       - mat                           = material properties
%       
%       OUTPUT: 
%       - K_fem                         = Stiffness Matrix in (kx,y,z,omega) for one combination (kx,omega).
%       
%       DESCRIPTION: 
%       The function calculates the stiffness matrix K_fem for one
%       combination (kx,Omega). Only the element stiffnesses for the Elements
%       that are located under the ground surface are calculated.
%
%       Assemling direct in loop not over assem function
%       -> used for MEX file generation
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_stiff_fem_cyl)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [K_fem] = f_stiff_fem_cyl(calc_sys,calc_struct,num_elem,Edof,div,NDOF,elem_coord,D_fem_soil,rho_fem_soil,H,radius,Omega,Kx)

% FEM Stiffness matrix
% Definition of Gauß Points
nPoints             = 2;
[wp,zeta,eta,ngp]   = f_gauss_points_2d(nPoints);
ep                  = [2 1];

% Predefinition of K_fem
K_fem = complex(eye(NDOF));

% Calculation of element stiffness and assemling
switch calc_sys
case 'trench'
    % Calculate the Element Stiffnesses for all Elements that are located under
    % the ground surface -> z=>0 (for hc mesh directly fullfilled)
    for s=1:num_elem
        ey  = elem_coord(s,[3 6 9 12]);
        ez  = elem_coord(s,[4 7 10 13]);
        Ke = f_planequad(ey,ez,ep,ngp,wp,zeta,eta,Kx,Omega,D_fem_soil,rho_fem_soil);   % soil fem
        K_fem (Edof(s,2:end),Edof(s,2:end)) = K_fem (Edof(s,2:end),Edof(s,2:end)) + Ke;
    end


case 'tunnel'
% Calculate the Element Stiffnesses for all Elements
    for s=1:num_elem
        ey = elem_coord(s,[3 6 9 12]);
        ez = elem_coord(s,[4 7 10 13]);
        % => Only soil fem used as in case of homog. filled tunnel.
        Ke = f_planequad(ey,ez,ep,ngp,wp,zeta,eta,Kx,Omega,D_fem_soil,rho_fem_soil); % soil fem 
        K_fem (Edof(s,2:end),Edof(s,2:end)) = K_fem (Edof(s,2:end),Edof(s,2:end)) + Ke;
    end

end
end



