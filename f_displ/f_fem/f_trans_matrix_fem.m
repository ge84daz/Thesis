%__________________________________________________________________________
%
%        f_trans_matrix_fem.m    
%
%       INPUT: 
%       - Nphi_t            = number of points/fouriermembers on cyl. surface
%       - phi               = circumferential angle
%       - mesh_type         = half_cyl or full_cyl
%
%       OUTPUT: 
%       - T12, T12_inv      = Transformation matrix for stiffness matrix of FEM substructure
%
%       DESCRIPTION: 
%       For the coupling of the ITM and the FEM subsystems the stiffnesses on the
%       cylindrical coupling surface have to be described in the same reference
%       system (here: ITM reference system).
%
%       | 1/ds T^-1 K_GammaGamma_FE T     1/ds T^-1 K_GammaOmega_FE | |u(kx,n,r=R,omega)|
%       |           K_OmegaGamma_FE T               K_OmegaOmgea_FE | |u(kx,y,z,omega)  |
% 
%       The function sets up the transformation matrices T12 and T12_inv 
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_trans_matrix_fem)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [dis_fem] = f_trans_matrix_fem(dis_itm,dis_fem,geo,calc_cont)

% Initialize
Nphi_t      = dis_itm.Nphi_t;
radius      = geo.radius;

switch calc_cont.mesh_type
    case 'full_cyl'
        phi = dis_fem.phi;
    case 'half_cyl'
        phi = geo.phi_LR_ges.';
end

%% Transformation of stiffnessmatrix in polar Basis (like ITM System)
% Transformationmatrix for cartesian -> polar (T1)
T1 = zeros(3*Nphi_t,3*Nphi_t);
for sc1 = 1:Nphi_t
    T1((sc1-1)*3+1:(sc1-1)*3+3,(sc1-1)*3+1:(sc1-1)*3+3) = ...
      [1       0              0          
       0  -sin(phi(sc1))  -cos(phi(sc1))   
       0   cos(phi(sc1))  -sin(phi(sc1))];
end

% Transformation polar -> Fouriercoefficients (T2)
T2 = complex(zeros(3*Nphi_t,3*Nphi_t));
n  = -Nphi_t/2:1:Nphi_t/2-1;

for sc2=1:Nphi_t
    T2((sc2-1)*3+1,1:3:3*Nphi_t) = exp(1i*n*phi(sc2));
    T2((sc2-1)*3+2,2:3:3*Nphi_t) = exp(1i*n*phi(sc2));
    T2((sc2-1)*3+3,3:3:3*Nphi_t) = exp(1i*n*phi(sc2));
end

T12     = T1*T2;
T12_inv = T12\eye(3*Nphi_t,3*Nphi_t);
ds      = ((2*pi*radius)/Nphi_t);

%% Output
dis_fem.T12     = T12;
dis_fem.T12_inv = T12_inv;
dis_fem.ds      = ds;

end