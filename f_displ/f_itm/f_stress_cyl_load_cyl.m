%__________________________________________________________________________
%
%        f_stress_cyl_load_cyl.m    
%
%       INPUT: 
%       - n_t = series member on cyl bound
%       - kx  = wavenumber in x-direction
%
%       OUTPUT: 
%       - the functions calculates the following matrix for one given series member:
%       
%       stress_LR_PLR_temp =
%       |           |           n_t(.)        |
%       |-----------|-------|---------|-------|
%       |           |  Prr  |  Prphi  |  Prx  |
%       |-----------|-------|---------|-------|
%       |u_x        | ...   |   ...   |  ...  |
%       |u_r        | ...   |   ...   |  ...  |
%       |u_phi      | ...   |   ...   |  ...  |
%       |sigma_rr   | ...   |   ...   |  ...  |
%       |sigma_rphi | ...   |   ...   |  ...  |
%       |sigma_tx   | ...   |   ...   |  ...  |
%
%       the resulting matrix is then stored in the total matrix such that:
%
%        stress_LR_PLR =
%       |          |           |     n_t(1)    |     n_t(2)    |     n_t(3)     |  ...  
%       |          |-----------|---------------|---------------|----------------|-------
%       |          |           | Prr Prphi Prx | Prr Prphi Prx | Prr Prphi Prx  |  ...  
%       |----------|-----------|---------------|---------------|----------------|-------
%       |   n_t(1) |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |u_r        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |u_phi      | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |sigma_rr   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |sigma_rphi | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |          |sigma_rx   | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%       |----------|-----------|---------------|---------------|----------------|-------
%       |   n_t(2) |u_x        | ...  ...  ... | ...  ...  ... |  ...  ...  ... |  ...  
%           :           :              :               :                :           : 
%
%       DESCRIPTION: 
%       The function determines the displ/stress at the cyl bound
%       due to unit loads (Prr, Prphi, Prx) applied at the cyl bound
%       Note:  Diss Frühe      - sigma = K*C and u = H*C
%              Diss Freisinger - sigma = S*C and u = U*C 
%              and K matrix = [H K]
%
%       REMARK: 
%       Original Author: Georg Frühe (span_tun_last_tun)
%       Modified:   Tom Hicks                 
%       Date:       22-09-2023 
%       Changed:    22-09-2023 - Hicks   
%__________________________________________________________________________


function [stress_LR_PLR_temp] = f_stress_cyl_load_cyl(n,kx,omega,mue2, cp2, cs2,r0)

% Predefinitions
stress_LR_PLR_temp  = complex(zeros(6,3));
H                   = complex(zeros(3,3));
K                   = complex(zeros(3,3));
coder.extrinsic('besselh');

% wavenumbers for p- and s-wave for soil of hs
kp = omega/cp2;
ks = omega/cs2;

% dimensionless parameters as input for hankel functions
kalpha = sqrt(kp^2-kx^2);
kbeta  = sqrt(ks^2-kx^2);

% Precompile H-matrix (relation between displacements and wave amplitudes)
besselh_n_kalpha_r0  = complex(zeros(1,1));
besselh_n1_kalpha_r0 = complex(zeros(1,1));
besselh_n_kbeta_r0   = complex(zeros(1,1));
besselh_n1_kebta_r0  = complex(zeros(1,1));

% get hankel function of first order for n and n+1 (needed for K)
% hankel function of second order not needed - diss frühe p.36
besselh_n_kalpha_r0  = besselh(n,kalpha*r0);
besselh_n1_kalpha_r0 = besselh(n+1,kalpha*r0);
besselh_n_kbeta_r0   = besselh(n,kbeta*r0);
besselh_n1_kebta_r0  = besselh(n+1,kbeta*r0);

% for u_x (Element H11, H12, H13) (Appendix Diss Frühe p. 126 or Diss Freisinger p.248)
H(1,1) = 1i*kx*besselh_n_kalpha_r0;
H(1,2) = 0;
H(1,3) = kbeta^2*besselh_n_kbeta_r0;
% for u_r (Element H21, H22, H23)
H(2,1) = n/r0*besselh_n_kalpha_r0-kalpha*besselh_n1_kalpha_r0;
H(2,2) = 1i*n/r0*besselh_n_kbeta_r0;
H(2,3) = 1i*kx*n/r0*besselh_n_kbeta_r0-1i*kx*kbeta*besselh_n1_kebta_r0;
% for u_phi (Element H31, H32, H33)
H(3,1) = 1i*n/r0*besselh_n_kalpha_r0;
H(3,2) = -n/r0*besselh_n_kbeta_r0+kbeta*besselh_n1_kebta_r0;
H(3,3) = -kx*n/r0*besselh_n_kbeta_r0;

% Evaluation of the part of the K-matrix which is corresponding to the implied boundary conditions 
% see diss freisinger p.245 component for sigma_rr
% for sigma_rr (Element K21, K22, K23)
K(1,1) = ((n^2-n)/r0^2+kx^2-1/2*ks^2)*besselh_n_kalpha_r0+kalpha/r0*besselh_n1_kalpha_r0;
K(1,2) = 1i*(n^2-n)/r0^2*besselh_n_kbeta_r0-1i*kbeta*n/r0*besselh_n1_kebta_r0;
K(1,3) = (-1i*kx*kbeta^2+1i*kx*(n^2-n)/r0^2)*besselh_n_kbeta_r0+1i*kx*kbeta/r0*besselh_n1_kebta_r0;
% for sigma_rphi (Element K51, K52, K53)
K(2,1) = 1i*(n^2-n)/r0^2*besselh_n_kalpha_r0-1i*kalpha*n/r0*besselh_n1_kalpha_r0;
K(2,2) = (1/2*kbeta^2-(n^2-n)/r0^2)*besselh_n_kbeta_r0-kbeta/r0*besselh_n1_kebta_r0;
K(2,3) = -kx*(n^2-n)/r0^2*besselh_n_kbeta_r0+kbeta*kx*n/r0*besselh_n1_kebta_r0;
% for sigma_rx (Element K41, K42, K43)
K(3,1) = 1i*kx*n/r0*besselh_n_kalpha_r0-1i*kalpha*kx*besselh_n1_kalpha_r0;
K(3,2) = -1/2*kx*n/r0*besselh_n_kbeta_r0;
K(3,3) = 1/2*n*(kbeta^2-kx^2)/r0*besselh_n_kbeta_r0-1/2*kbeta*(kbeta^2-kx^2)*besselh_n1_kebta_r0;

% factor 2*mue within formulations - see diss freisinger p.245 ff.
K = K*2*mue2;

% unit load states
Prr   = [-1,0,0].';
Prphi = [0,-1,0].';
Prx   = [0,0,-1].';

% Konditionierung und Lösung des Gleichungssystems nach den Koeffizienten
% für alle Lastzustände Prr=1, Prphi=1 und Prx=1
cond=abs(besselh_n_kbeta_r0/besselh_n_kalpha_r0);
K(:,1)     = cond*K(:,1);
C_Prr      = (K\Prr);
C_Prphi    = (K\Prphi);
C_Prx      = (K\Prx);
C_Prr(1)   = C_Prr(1)  *cond;
C_Prphi(1) = C_Prphi(1)*cond;
C_Prx(1)   = C_Prx(1)  *cond;

% FOURIER-transformed stresses and displacements in (n,r,kx,omega)
stress_LR_PLR_temp(1:3,:)=[H*C_Prr  H*C_Prphi  H*C_Prx];
% stresses are corresponding to the unit load case = -1
stress_LR_PLR_temp(4:6,:)=-eye(3);
