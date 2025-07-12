%__________________________________________________________________________
%
%       f_displ_hs.m    
%
%       INPUT: 
%       kx,ky                     = Inkrements in transformed domain
%       omega                     = Frequency
%       Pz                        = Vertical Load for one Triple (kx,ky,omega)
%       mue1,lambda1,cp1,cs1      = Material Parameter of soil
%       z                         = discretization over depth
%
%
%       OUTPUT: 
%       u_z0      = uz at halfspace surface z=0 in (kx,ky,z,omega) 
%       u_zH_1L   = not evaluated
%       uzPz_zvar = uz within halfspace over depth z in (kx,ky,z,omega) 
%
%       DESCRIPTION: 
%       The function calculates the stresses and displ. at z=0 or over z
%       for one Triple (kx,ky,omega) of a homogeneouse halfspace.
%       ux,uy,uz for each load direction Px,Py,Pz are returned.
%
%       REMARK: 
%       Original Author: Julian Freisinger  (f_displ_hs)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    24-09-2023 - Hicks   
%__________________________________________________________________________

function [u_z0] = f_displ_hs_paralell(varargin)

% Standard input
kx      = varargin{1};
ky      = varargin{2};
omega   = varargin{3};
P_hs1L  = varargin{4};
mat     = varargin{5};

% evaluation over depth if z given as input
if nargin == 6
   z = varargin{6}; 
end

mu     = mat.hs.mu;
lambda = mat.hs.lambda;
cp     = mat.hs.cp;
cs     = mat.hs.cs;

% no case study because damping with -2Di in material and only calculated for
% negative frequencies
kp = omega/cp;
ks = omega/cs;

kr2 = kx^2+ky^2;
lambda1 = sqrt(kr2-kp^2);
lambda2 = sqrt(kr2-ks^2);


%% approach 1 - Determine unknowns A2 Bx2 By2
% line 1-3: matrix S - Diss Freisinger Appendix
% line 4-6: matrix U - Diss Frühe eq. 3.17 mit kr^2 = kx^2 + ky^2
% with corrected value for K[4,1] - Fehler in Diss Frühe eq. 3.17
K = [
    1i*kx                               0                   lambda2             % ux
    1i*ky                               -lambda2            0                   % uy
    -lambda1                            -1i*ky              1i*kx               % uz
    2*lambda1^2-lambda/mu*kp^2          2*1i*ky*lambda2     -2*1i*kx*lambda2    % sigma zz
    -2*1i*ky*lambda1                    lambda2^2+ky^2      -kx*ky              % sigma zy
    -2*1i*kx*lambda1                    kx*ky               -(lambda2^2+kx^2)]; % sigma zx

K(4:6,:) = mu*K(4:6,:);

% Unit load states for each wavenumber (only load in z-direction)
Pz = [-P_hs1L;0;0];
% Determine unknown coefficients C = [A2 Bx2 By2]
Cz = K(4:6,:)\Pz;
% Determine displacements at z=0
u_z0 = zeros(3,1);
% ux,uy,uz due to Pz at z=0 (u = U*C)
u_z0(1:3) = K(1:3,:)*Cz(1:3); 




end







