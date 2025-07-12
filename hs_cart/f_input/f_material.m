%__________________________________________________________________________
%
%       f_material.m    
%
%       INPUT: 
%       mat           = Object including all material properties
%       calc          = string defining which material to be used
%
%       OUTPUT: 
%       mat           = Object including all material properties
%
%       DESCRIPTION: 
%       Define all material properties for different materials used in model.
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_material)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________


function [mat]=f_material(mat,calc)


% 1.)   Material definition
mat.def = cell(1,1);

% 1.1.) Soil
mn = 1;
mat.def{mn}.name  = 'soft_v01';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 260e5*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]- E(1+isign(w)*2D) (Frühe eq. 2.6)
mat.def{mn}.nu    = 0.3;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

mn = mn +1;
mat.def{mn}.name  = 'stiff_v01';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 260e6*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.3;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 1600;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%% 2.)   Material assignment

% Find material number corresponding to material name input
for nmat = 1:mn
    if strcmpi(mat.def{1,nmat}.name, calc.mat_hs);             mn_hs             = nmat; end
end 
    
% assign material
mat.hs          = mat.def{1,mn_hs};

clear mn
end


%% 3.)   Computation of secondary material parameter

function [mat] = f_det_mat_par(mat,mn)
% input
E   = mat.def{mn}.E;
nu  = mat.def{mn}.nu;
rho = mat.def{mn}.rho;

% calc mat param
G      = E./(2*(1+nu));
mu     = G;
lambda = 2*G.*nu./(1-2*nu);
cp     = sqrt((lambda+2*mu)./rho);
cs     = sqrt(mu./rho);
cr     = (0.87+1.12*nu)/(1+nu)*cs;

% output
mat.def{mn}.G      = G;
mat.def{mn}.mu     = mu;
mat.def{mn}.lambda = lambda ;
mat.def{mn}.cp     = cp;
mat.def{mn}.cs     = cs;
mat.def{mn}.cr     = cr;
end

function [D]   = hooke(ptype,E,v)


if ptype==1
    Dm=E/(1-v^2)*[1  v   0;
        v  1   0;
        0  0 (1-v)/2];
elseif ptype==2
    Dm=E/(1+v)/(1-2*v)*[1-v  v    v        0;
        v  1-v   v        0;
        v   v   1-v       0;
        0   0    0   (1-2*v)/2];
elseif ptype==3
    Dm=E/(1+v)/(1-2*v)*[1-v  v    v        0;
        v  1-v   v        0;
        v   v   1-v       0;
        0   0    0   (1-2*v)/2];
elseif ptype==4
    Dm=E/(1+v)/(1-2*v)*[1-v  v    v    0    0    0;
        v  1-v   v    0    0    0;
        v   v   1-v   0    0    0;
        0   0    0 (1-2*v)/2    0    0;
        0   0    0    0 (1-2*v)/2    0;
        0   0    0    0    0  (1-2*v)/2];
else
    error('Error ! Check first argument, ptype=1,2,3 or 4 allowed')
    return
end

D = Dm;

end
