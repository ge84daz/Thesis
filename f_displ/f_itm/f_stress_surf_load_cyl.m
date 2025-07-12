%__________________________________________________________________________
%
%        f_stress_surf_load_cyl.m    
%
%       INPUT: 
%       - n_t = relevant Fourier members of cyl boundary
%       - kx, ky, omega
%       - r0 = radius of cylinder
%       - r = r-coordiante of points on hs surf
%
%       OUTPUT: 
%       - the functions calculates the following matrix for one given n_t(.):
%       
%       stress_OF_pol_T1_temp =
%       |             |                    n_t(.)                    |
%       |------------------------------------------------------------|
%       |             |         Prr         |     Prphi      |  Prx  |
%       |             | y(1) y(2) ...  y(n) | y(1) ...  y(n) |  ...  |
%       |------------------------------------------------------------|
%       |u_x          | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |u_r          | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |u_phi        | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_xx     | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_rr     | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_phiphi | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_rx     | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_rphi   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%       |sigma_phix   | ...  ...  ...  ...  | ...  ...  ...  |  ...  |
%
%       whereas y(.) are all points on the halfspace surface
%       the resulting matrix is then stored in the total matrix such that:
%
%       stress_OF_pol_T1 =
%                     |      n_t(1)     |      n_t(2)     |      n_t(3)      |  ...  |
%        -------------|-----------------|-----------------|------------------|-------|
%                     | Prr  Prphi  Prx | Prr  Prphi  Prx | Prr  Prphi  Prx  |  ...  |
%       --------------|-----------------|-----------------|------------------|-------|
%       u_x           | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%       u_r           | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%       u_phi         | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%       sigma_xx      | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%       sigma_rr      | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%       sigma_phiphi  | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%       sigma_rx      | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%       sigma_rphi    | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%       sigmaphix     | ...   ...   ... | ...   ...   ... | ...   ...   ...  |  ...  |
%
%       DESCRIPTION: 
%       The function determines the displ/stress at the hs surf
%       due to unit loads (Prr, Prphi, Prx) applied at the cyl boundary
%       Note:  Diss Frühe      - sigma = K*C and u = H*C
%              Diss Freisinger - sigma = S*C and u = U*C 
%              and K matrix = [H K]
%
%       REMARK: 
%       Original Author: Georg Frühe (span_ofl_last_tun)
%       Modified:   Tom Hicks                 
%       Date:       22-09-2023 
%       Changed:    22-09-2023 - Hicks   
%__________________________________________________________________________

function [stress_OF_pol_T1_temp] = f_stress_surf_load_cyl(nn,kx,omega,r,phi,mue2, cp2, cs2, r0)

% For MEX file generation -> besselh not supported until now
coder.extrinsic('besselh');

% Note: bessel function of third order = Hankel function of first and second order

% material parameters of halfspace
% wavenumbers for p- and s-wave
kp = omega/cp2;
ks = omega/cs2;

% dimensionless parameters
% input parameters for bessel functions
kalpha = sqrt(kp^2-kx^2); % Diss Freisinger: kalpha = k1
kbeta  = sqrt(ks^2-kx^2); % Diss Freisinger: kbeta  = k2

% Evaluation of the part of the K-matrix which is corresponding to the implied boundary conditions 
% prestore
K     = complex(zeros(3,3));


% get hankel function of first order for n and n+1 (needed for K)
% hankel function of second order not needed - diss frühe p.36
besselh_n_kalpha_r0  = besselh(nn,kalpha*r0);
besselh_n1_kalpha_r0 = besselh(nn+1,kalpha*r0);
besselh_n_kbeta_r0   = besselh(nn,kbeta*r0);
besselh_n1_kbeta_r0  = besselh(nn+1,kbeta*r0);

% see diss freisinger p.245 component for sigma_rr
% für sigma_rr (Elemente K21, K22, K23)
K(1,1) = ((nn^2-nn)/r0^2+kx^2-1/2*ks^2)*besselh_n_kalpha_r0 + kalpha/r0*besselh_n1_kalpha_r0;
K(1,2) = 1i*(nn^2-nn)/r0^2*besselh_n_kbeta_r0-1i*kbeta*nn/r0*besselh_n1_kbeta_r0;
K(1,3) = 1i*kx*(-kbeta^2+(nn^2-nn)/r0^2)*besselh_n_kbeta_r0+1i*kx*kbeta/r0*besselh_n1_kbeta_r0;

% see diss freisinger p.247 component for sigma_rphi
% für sigma_rphi (Elemente K51, K52, K53)
K(2,1) = 1i*(nn^2-nn)/r0^2*besselh_n_kalpha_r0-1i*kalpha*nn/r0*besselh_n1_kalpha_r0;
K(2,2) = (1/2*kbeta^2-(nn^2-nn)/r0^2)*besselh_n_kbeta_r0-kbeta/r0*besselh_n1_kbeta_r0;
K(2,3) = -kx*(nn^2-nn)/r0^2*besselh_n_kbeta_r0+kbeta*kx*nn/r0*besselh_n1_kbeta_r0;

% see diss freisinger p.246 component for sigma_rx
% für sigma_rx (Elemente K41, K42, K43)
K(3,1) = 1i*kx*nn/r0*besselh_n_kalpha_r0-1i*kalpha*kx*besselh_n1_kalpha_r0;
K(3,2) = -1/2*kx*nn/r0*besselh_n_kbeta_r0;
K(3,3) = 1/2*nn*(kbeta^2-kx^2)/r0*besselh_n_kbeta_r0-1/2*kbeta*(kbeta^2-kx^2)*besselh_n1_kbeta_r0;

% factor 2*mue within formulations - see diss freisinger p.245 ff.
K = K*2*mue2;

% unit load states
Prr   = [-1,0,0].';
Prphi = [0,-1,0].';
Prx   = [0,0,-1].';

% Konditionierung und Lösung des Gleichungssystems nach den Koeffizienten
% für alle Lastzustände P_rr=1, P_rphi=1 und P_rx=1
cond       = abs(besselh_n_kbeta_r0/besselh_n_kalpha_r0); % why? 
K(:,1)     = cond*K(:,1);
C_Prr      = (K\Prr).';
C_Prphi    = (K\Prphi).';
C_Prx      = (K\Prx).';
C_Prr(1)   = C_Prr(1)*cond;
C_Prphi(1) = C_Prphi(1)*cond;
C_Prx(1)   = C_Prx(1)*cond;

% values of the hankel functions for all points r,phi at hs surf
H_n_a  = complex(zeros(1,size(r,2)));
H_n1_a = complex(zeros(1,size(r,2)));
H_n_b  = complex(zeros(1,size(r,2)));
H_n1_b = complex(zeros(1,size(r,2)));

% precompute - needed for matrix H 
H_n_a    = besselh(nn,kalpha*r);
H_n_a_r  = 1./r.*H_n_a;
H_n_a_r2 = 1./r.^2.*H_n_a; 
H_n1_a   = besselh(nn+1,kalpha*r);
H_n1_a_r = 1./r.*H_n1_a;
H_n_b    = besselh(nn,kbeta*r);
H_n_b_r  = 1./r.*H_n_b;
H_n_b_r2 = 1./r.^2.*H_n_b;
H_n1_b   = besselh(nn+1,kbeta*r);
H_n1_b_r = 1./r.*H_n1_b;

% Elements of H-Matrix (Appendix Diss Frühe p. 126 or Diss Freisinger p.248)
% for u_x (Elements H11, H12, H13)
H1=[1i*kx*H_n_a;
    0*r;
    kbeta^2*H_n_b];
% for u_r (Elements H21, H22, H23)
H2=[nn*H_n_a_r-kalpha*H_n1_a;
    1i*nn*H_n_b_r;
    1i*kx*nn*H_n_b_r-1i*kx*kbeta*H_n1_b];
% for u_phi (Elements H31, H32, H33)
H3=[1i*nn*H_n_a_r;
    -nn*H_n_b_r+kbeta*H_n1_b;
    -kx*nn*H_n_b_r];

% Evaluation of total K-Matrix
% note that second order hankel is not needed as describing incoming waves
% from infinity (comparison homogeneous halfspace with A1, Bx1, By1 coefficients)
% (see Diss Freisinger p.248)
% for sigma_xx (Elements K11, K12, K13)
K1=[(kalpha^2-1/2*ks^2)*H_n_a;
    0*r;
    1i*kx*kbeta^2*H_n_b]*2*mue2;
% for sigma_rr (Elements K21, K22, K23)
K2=[(nn^2-nn)*H_n_a_r2+(kx^2-1/2*ks^2)*H_n_a+kalpha*H_n1_a_r;
    1i*(nn^2-nn)*H_n_b_r2-1i*kbeta*nn*H_n1_b_r;
    -1i*kx*kbeta^2*H_n_b+1i*kx*(nn^2-nn)*H_n_b_r2+1i*kx*kbeta*H_n1_b_r]*2*mue2;
% for sigma_phiphi (Elements K31, K32, K33)
K3=[-(nn^2-nn)*H_n_a_r2+(kalpha^2+kx^2-1/2*ks^2)*H_n_a-kalpha*H_n1_a_r;
    -1i*(nn^2-nn)*H_n_b_r2+1i*kbeta*nn*H_n1_b_r;
    -1i*kx*(nn^2-nn)*H_n_b_r2-1i*kbeta*kx*H_n1_b_r]*2*mue2;
% for sigma_rx (Elements K41, K42, K43)
K4=[1i*kx*nn*H_n_a_r-1i*kalpha*kx*H_n1_a;
    -1/2*kx*nn*H_n_b_r;
    1/2*nn*(kbeta^2-kx^2)*H_n_b_r-1/2*kbeta*(kbeta^2-kx^2)*H_n1_b]*2*mue2;
% for sigma_rphi (Elements K51, K52, K53)
K5=[1i*(nn^2-nn)*H_n_a_r2-1i*kalpha*nn*H_n1_a_r;
    1/2*kbeta^2*H_n_b-(nn^2-nn)*H_n_b_r2-kbeta*H_n1_b_r;
    -kx*(nn^2-nn)*H_n_b_r2+kbeta*kx*nn*H_n1_b_r]*2*mue2;
% for sigma_phix (Elements K61, K62, K63)
K6=[-kx*nn*H_n_a_r;
    -1i/2*kx*nn*H_n_b_r+1i/2*kbeta*kx*H_n1_b;
    1i/2*nn*(kbeta^2-kx^2)*H_n_b_r]*2*mue2;

% FOURIER-transformed stresses and displacments in (n,r,kx,omega)
% for all unit load states P_rr=1, P_rphi=1 und P_rx=1
sigmaFT=[C_Prr*H1   C_Prphi*H1   C_Prx*H1     % u_x    
         C_Prr*H2   C_Prphi*H2   C_Prx*H2     % u_r
         C_Prr*H3   C_Prphi*H3   C_Prx*H3     % u_phi
         C_Prr*K1   C_Prphi*K1   C_Prx*K1     % sigma_x
         C_Prr*K2   C_Prphi*K2   C_Prx*K2     % sigma_r
         C_Prr*K3   C_Prphi*K3   C_Prx*K3     % sigma_phi
         C_Prr*K4   C_Prphi*K4   C_Prx*K4     % tau_rx
         C_Prr*K5   C_Prphi*K5   C_Prx*K5     % tau_rphi
         C_Prr*K6   C_Prphi*K6   C_Prx*K6];   % tau_phix

% exponential functions for transformation in phi-direction (Umfangsrichtung)
% phi = phi-coordiantes of points on hs surf
% exponential decrease of wave amplitudes with increasing distance
expo_phi = [ones(9,1)*exp(1i*nn*phi)  ones(9,1)*exp(1i*nn*phi)  ones(9,1)*exp(1i*nn*phi)];

% stresses and displacements in original space (phi,r,kx,omega-Raum)
stress_OF_pol_T1_temp = sigmaFT.*expo_phi;

