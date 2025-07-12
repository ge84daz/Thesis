%__________________________________________________________________________
%
%        f_stress_cyl_load_surf.m    
%
%       INPUT: 
%       - N_t = relevant Fourier member at hs surf
%       - kx, ky, omega
%
%       OUTPUT: 
%       - the functions calculates the following matrix for one given N_t(.):
%       
%       stress_LR_kart_T1_temp =
%       |           |                    N_t(.)                    |
%       |-----------|---------------------|----------------|-------|
%       |           |         Pzz         |      Pzy       |  Pzx  |
%       |           | z(1) z(2) ...  z(n) | z(1) ...  z(n) |  ...  |
%       |-----------|---------------------|----------------|-------|
%       |u_x        | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |u_y        | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |u_z        | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_xx   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_yy   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_zz   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_xy   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_zy   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_zx   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%
%       where as z(n) are the points on the cylindrical surface.
%       the resulting matrix is then stored in the total matrix such that:
%
%        stress_LR_kart_T1 =
%       |           |     N_t(1)    |     N_t(2)    |     N_t(3)     |  ...    
%       |-----------|---------------|---------------|----------------|------
%       |           | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx  |  ... 
%       |-----------|---------------|---------------|----------------|------
%       |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%       |u_y        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%       |u_z        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%       |sigma_xx   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%       |sigma_yy   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%       |sigma_zz   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%       |sigma_xy   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%       |sigma_zy   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%       |sigma_zx   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...    
%
%       DESCRIPTION: 
%       The function determines the displ/stress at the cylindrical
%       surface due to unit loads (Pzz, Pzy, Pzx) applied at the halfspace
%       surface
%
%       REMARK: 
%       Original Author: Georg Frühe (span_tun_last_ofl)
%       Modified:   Tom Hicks                 
%       Date:       22-09-2023 
%       Changed:    22-09-2023 - Hicks   
%__________________________________________________________________________

function [stress_LR_kart_T1_temp] = f_stress_cyl_load_surf(kx,ky,omega,y,z,mue2,lambda2,cp2,cs2)

% Predefinitions
K1 = complex(zeros(3,3));

% material parameters of halfspace soil
% wavenumbers for p- and s-waves
kp = omega/cp2;
ks = omega/cs2;

% Substitution
kr_2 = kx^2+ky^2;

% exponential coefficients
lambda1_2 = sqrt(kr_2-kp^2);
lambda2_2 = sqrt(kr_2-ks^2);

% K - Matrix for hom. halfspace - Diss Freisinger p. 241 - mistake within K(6.1)
K=[ 1i*kx                          0                       lambda2_2              % u_x
    1i*ky                          -lambda2_2              0                      % u_y
   -lambda1_2                      -1i*ky                  1i*kx                  % u_z
   -(2*kx^2+lambda2/mue2*kp^2)     0                       2*1i*kx*lambda2_2      % sigma_xx
   -(2*ky^2+lambda2/mue2*kp^2)     -2*1i*ky*lambda2_2      0                      % sigma_yy
    2*kr_2-ks^2                    2*1i*ky*lambda2_2       -2*1i*kx*lambda2_2     % sigma_zz
   -2*kx*ky                        -1i*kx*lambda2_2        1i*ky*lambda2_2        % sigma_xy
   -2*1i*ky*lambda1_2              lambda2_2^2+ky^2       -kx*ky                  % sigma_yz
   -2*1i*kx*lambda1_2              kx*ky                  -(lambda2_2^2+kx^2)];   % sigma_zx

K(4:9,:)=K(4:9,:)*mue2;

% Unit load states
Pzz=[-1,0,0].';  % Pzz
Pzy=[0,-1,0].';  % Pzy
Pzx=[0,0,-1].';  % Pzx

% Solve system of equations for coefficients for all unit load states
% Boundary condition is given by unit load states
K1(:,:) = K([6 8 9],:);
C_Pzz   = K1\Pzz;
C_Pzy   = K1\Pzy;
C_Pzx   = K1\Pzx;

% exponential functions in z-direction for all unit load states
% z = z-coordiantes of points on cylindrical surface
% with unknowns A2, Bx2 and By2 for each unit load state at each position z 
expo_z = [C_Pzz(1)*exp(-lambda1_2*z)  C_Pzy(1)*exp(-lambda1_2*z)  C_Pzx(1)*exp(-lambda1_2*z)
          C_Pzz(2)*exp(-lambda2_2*z)  C_Pzy(2)*exp(-lambda2_2*z)  C_Pzx(2)*exp(-lambda2_2*z)   
          C_Pzz(3)*exp(-lambda2_2*z)  C_Pzy(3)*exp(-lambda2_2*z)  C_Pzx(3)*exp(-lambda2_2*z)]; 
%        expo_z =
%       |     |         Pzz         |      Pzy       |  Pzx  |
%       |----------------------------------------------------|
%       |     | z(1) z(2) ...  z(n) | z(1) ...  z(n) |  ...  |
%       |----------------------------------------------------|
%       |A2   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |Bx2  | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |By2  | ...  ...  ...  ...  | ...  ...  ...  |  ...  |


% FOURIER-transformed Stress and Displacements in image space (ky,z,kx,omega)
sigmaFT = K*expo_z;


% exponential functions for transformation in y-direction
% y = y-coordiantes of points on cylindrical surface
% exponential decrease of wave amplitudes with increasing distance
expo_y = [ones(9,1)*exp(1i*ky*y)  ones(9,1)*exp(1i*ky*y)  ones(9,1)*exp(1i*ky*y)];

% Stress and Displacements in original space (y,z,kx,omega)
stress_LR_kart_T1_temp = sigmaFT.*expo_y;

end % function

