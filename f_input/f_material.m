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
%       - Define all material properties for different materials used in model.
%       - Loss factor eta = 2*Damping ratio
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_material)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________

function [mat] = f_material(mat,calc)

%% 1.)   Material definition
mat.def = cell(1,25);

%% 1.1.) Soil
mn = 1;
mat.def{mn}.name  = 'soft_v01';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 260e5*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.3;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'soft_v02';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 1.04e8*(1-mat.def{mn}.eta*1i);           % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.3;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'stiff_v01';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 260e6*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.3;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 1600;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'stiff_v02';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 260e6*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.3;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'stiff_v03';                             % = 'comp_mh_tun_IL'
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 6e7*(1-mat.def{mn}.eta*1i);              % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.3;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'stiff_v04';                             % comp D.V.Jones mat 1
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 269e6*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.257;                                   % Poisson Ratio [-]
mat.def{mn}.rho   = 1550;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'stiff_v05';                             % comp Bian
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 188e6*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.3;                                   % Poisson Ratio [-]
mat.def{mn}.rho   = 1800;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param


%% 1.2.) Concrete
%
mn = mn +1;
mat.def{mn}.name  = 'concrete_v01';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 34e9*(1-mat.def{mn}.eta*1i);             % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_v02';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 9.6e7*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_v03';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 3.84e8*(1-mat.def{mn}.eta*1i);           % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_v04';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 1.54e9*(1-mat.def{mn}.eta*1i);           % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2000;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_v05';                          % as used in ma_wagner, diss mh tun shell
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 34e9*(1-mat.def{mn}.eta*1i);             % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2600;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_v06';                          % as used Dijckmans
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 30e9*(1-mat.def{mn}.eta*1i);             % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2400;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_stiff';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 3e15*(1-mat.def{mn}.eta*1i);             % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.3;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2600;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_stiff_m0';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 3e15*(1-mat.def{mn}.eta*1i);             % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 0.1;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_c50_60';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 37e9*(1-mat.def{mn}.eta*1i);             % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2500;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_v07';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.02;                                  % Loss factor [-]
mat.def{mn}.E     = 3e10*(1-mat.def{mn}.eta*1i);             % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2500;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'concrete_v07';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 3.6e10*(1-mat.def{mn}.eta*1i);             % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.2;                                     % Poisson Ratio [-]
mat.def{mn}.rho   = 2400;                                    % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param


%% 1.4.) Sylomer
%
mn = mn +1;
mat.def{mn}.name  = 'sylomer_SR11';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 0.125;                                   % Loss factor [-]
mat.def{mn}.E     = 0.06*1e6*(1-mat.def{mn}.eta*1i);         % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.40;                                    % Poisson Ratio [-]
mat.def{mn}.rho   = 600;                                     % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'sylomer_mh';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 1*0.07;                                  % Loss factor [-]
mat.def{mn}.E     = 0.9e6*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.30;                                    % Poisson Ratio [-]
mat.def{mn}.rho   = 620;                                     % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param

%
mn = mn +1;
mat.def{mn}.name  = 'sylomer_jf';
mat.def{mn}.num   = mn;
mat.def{mn}.eta   = 2*0.05;                                  % Loss factor [-]
mat.def{mn}.E     = 5.5e5*(1-mat.def{mn}.eta*1i);            % Complex E-Modulus [N/m^2]
mat.def{mn}.nu    = 0.36;                                    % Poisson Ratio [-]
mat.def{mn}.rho   = 620;                                     % Density [kg/m^2]
mat.def{mn}.D     = hooke(4, mat.def{mn}.E, mat.def{mn}.nu); % Elasticity Matrix
[mat]             = f_det_mat_par(mat,mn);                   % determine other mat param



% write material name in 2nd row of mat.def
for nmat = 1:mn
    mat.def{2,nmat} = mat.def{1,nmat}.name;
end

%% 2.)   Material assignment

% Find material number corr. to material name input
for nmat = 1:mn
    if strcmpi(mat.def{1,nmat}.name, calc.mat_hs);             mn_hs             = nmat; end
    if strcmpi(mat.def{1,nmat}.name, calc.mat_fem_soil_T1);    mn_fem_soil_T1    = nmat; end
end

% assign material
mat.hs          = mat.def{1,mn_hs};
mat.fem_soil_T1 = mat.def{1,mn_fem_soil_T1};

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

