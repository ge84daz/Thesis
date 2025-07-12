%__________________________________________________________________________
%
%        f_trans_stiff_fem.m    
%
%       INPUT: 
%       K_fem             = cartesian FEM stiffnessmatrix
%       T12, T12_inv      = Transformation matrices
%       radius            = radius of cyl. FEM subsystem
%
%       OUTPUT: 
%       - K_fem_elem_pol    = stiffness matrix of FEM substructure
%
%       DESCRIPTION: 
%       For the coupling of the ITM and the FEM subsystems the stiffnesses on the
%       cylindrical coupling surface have to be described in the same reference
%       system (here: ITM reference system).
%
%       | 1/ds T^-1 K_GammaGamma_FE T     1/ds T^-1 K_GammaOmega_FE | |u(kx,n,r=R,omega)|
%       |           K_OmegaGamma_FE T               K_OmegaOmgea_FE | |u(kx,y,z,omega)  |
% 
%       The function sets up the transformation matrices and transforms the
%       respective members of K_fem -> K_fem_pol.
% 
%       REMARK: 
%       Original Author: Julian Freisinger (f_trans_stiff_fem)
%       Modified:   Tom Hicks                 
%       Date:       21-09-2023 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________


function [K_fem_pol]=f_trans_stiff_fem(K_fem,dis_fem)

% Initialize
T12        = dis_fem.T12;
T12_inv    = dis_fem.T12_inv;
ds         = dis_fem.ds;
num_nod_bd = dis_fem.num_nod_bd;

%  Transformation of stiffnessmatrix in polar Basis (like ITM System)
K_fem_pol = [(T12_inv)*K_fem(1:3*num_nod_bd,1:3*num_nod_bd)    *(T12)    (T12_inv)*K_fem(1:3*num_nod_bd,3*num_nod_bd+1:end);
                       K_fem(3*num_nod_bd+1:end,1:3*num_nod_bd)*(T12)              K_fem(3*num_nod_bd+1:end,3*num_nod_bd+1:end)];
        
% Distrubution of nodal loads over element length = Adaption on stresses in ITM stiffnessmatrix
K_fem_pol(1:3*num_nod_bd,:) = (1/ds)*K_fem_pol(1:3*num_nod_bd,:);

end
