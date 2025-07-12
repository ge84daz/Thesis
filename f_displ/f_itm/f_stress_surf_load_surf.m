%__________________________________________________________________________
%
%        f_stress_surf_load_surf.m    
%
%       INPUT: 
%       - kx, ky, omega
%
%       OUTPUT: 
%       - the functions calculates the following matrix for one given wavenumber ky:
%       
%       stress_OF_POF_temp =
%       |         |          ky(.)        |
%       |---------|-------|-------|-------|
%       |         |  Pzz  |  Pzy  |  Pzx  |
%       |---------|-------|-------|-------|
%       |u_x      | ...   |  ...  |  ...  |
%       |u_y      | ...   |  ...  |  ...  |
%       |u_z      | ...   |  ...  |  ...  |
%       |sigma_zz | ...   |  ...  |  ...  |
%       |sigma_zy | ...   |  ...  |  ...  |
%       |sigma_zx | ...   |  ...  |  ...  |
%
%       the resulting matrix is then stored in the total matrix such that:
%
%        stress_OF_POF =
%       |          |           |     ky(1)     |     ky(2)     |     ky(3)      |  ...  
%       |          |-----------|---------------|---------------|----------------|-------
%       |          |           | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx | Pzz  Pzy  Pzx  |  ...  
%       |----------|-----------|---------------|---------------|----------------|-------
%       |   ky(1)  |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |u_y        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |u_z        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |sigma_zz   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |sigma_zy   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |sigma_zx   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |----------|-----------|---------------|---------------|----------------|-------
%       |  ky(2)   |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%           :           :              :               :                :           : 
%
%       DESCRIPTION: 
%       The function determines the displ/stress at the halfspace
%       surface due to unit loads (Pzz, Pzy, Pzx) applied at the halfspace
%       surface
%
%       REMARK: 
%       Original Author: Georg Frühe (span_ofl_last_ofl)
%       Modified:   Tom Hicks                 
%       Date:       22-09-2023 
%       Changed:    22-09-2023 - Hicks   
%__________________________________________________________________________


function [stress_OF_POF_temp] = f_stress_surf_load_surf(kx, ky, omega, mue2, lambda2, cp2, cs2)

% Predefinitions
stress_OF_POF_temp = complex(zeros(6,3));
K1    = complex(zeros(3,3));
K2    = complex(zeros(3,3));

% Wavenumbers for p- and s-wave
kp=omega/cp2;
ks=omega/cs2;

% Substitution
kr_2=kx^2+ky^2;

% exponential coefficients
lambda1_2=sqrt(kr_2-kp^2);
lambda2_2=sqrt(kr_2-ks^2);

% Evaluate K-Matrix (note mistake in K(6.1) due to Diss Fruehe)
K=[ 1i*kx                           0                      lambda2_2             % u_x
    1i*ky                          -lambda2_2              0                     % u_y
   -lambda1_2                     -1i*ky                   1i*kx                 % u_z
   -(2*kx^2+lambda2/mue2*kp^2)     0                      2*1i*kx*lambda2_2      % sigma_x
   -(2*ky^2+lambda2/mue2*kp^2)    -2*1i*ky*lambda2_2       0                     % sigma_y
    2*kr_2-ks^2                    2*1i*ky*lambda2_2      -2*1i*kx*lambda2_2     % sigma_z
   -2*kx*ky                       -1i*kx*lambda2_2         1i*ky*lambda2_2       % sigma_xy
   -2*1i*ky*lambda1_2               lambda2_2^2+ky^2      -kx*ky                 % sigma_yz
   -2*1i*kx*lambda1_2               kx*ky                 -(lambda2_2^2+kx^2)];  % sigma_zx
K(4:9,:)=K(4:9,:)*mue2;

% unit load states
Pzz=[-1,0,0].';
Pzy=[0,-1,0].';
Pzx=[0,0,-1].';

% Solve for unknown coefficients
K1(:,:)=K([6 8 9],:);
C_Pzz=K1\Pzz;
C_Pzy=K1\Pzy;
C_Pzx=K1\Pzx;

% FOURIER-transformed stress/displ in (ky,z,kx,omega) space
K2(:,:)=K([1 2 3],:);

stress_OF_POF_temp(1:3,:)=[K2*C_Pzz  K2*C_Pzy  K2*C_Pzx];

% stress equals unit load states
stress_OF_POF_temp(4:6,:)=-eye(3);
