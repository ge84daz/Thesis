%__________________________________________________________________________
%
%        f_comb_displ_tr_f_z0.m    
%
%       INPUT: 
%       -displ_dir     = Displ. direction
%       -load_dir      = Displ. due to load on load_dir
%       -displ         = Displacements
%       -dis_itm       = Discretization ITM
%       -dis_fem       = Discretization FEM
%       -geo           = Geometry
%       -loading       = Load (amplitudes, distribution,...)
%       -mat           = Material Parameters
%       -plot_cont     = plot control parameters
%       -path          = path parameters
%       -calc          = Calculation parameters
%       -fev           = Evaluation frequency
% 
%
%       OUTPUT: 
%       - u_xyf_fem_xyz         = all FEM displ ux,uy,uz of cyl (boundary + inside) 
%       - u_kxyf_fem_xyz        = all FEM displ ux,uy,uz of cyl (boundary + inside) 
%       - u_x_midline_f_fem_x   = displ. ux on cyl midline
%       - u_y_midline_f_fem_y   = displ. uy on cyl midline
%       - u_z_midline_f_fem_z   = displ. uz on cyl midline
%
%       DESCRIPTION: 
%       The function combines the solution parts on the halfspace surface of the
%       ITM and FEM System for a halfspace with a cylindrical trench/slit and
%       displayes them in different plots, also compared to a reference solution
%       calculated using a layered halfspace with one layer.
%       Furthermore the Amplitude Reduction Factors AR are calculated for
%       different setups.
%       - only displ. in frequency domain plotted for negative frequency f
%       - if harmonic time dependent displ shall be plotted:
%         Time of maximum displ t_max_itm_fem -> equal for trench
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_comb_displ_tr_f_z0)
%       Modified:   Tom Hicks                 
%       Date:       09-11-2021 
%       Changed:    02-12-2024 - Hicks   
%__________________________________________________________________________


function [displ,geo_T1] = f_comb_displ_tr_f_z0(displ_dir,load_dir,displ,geo_T1,~,dis_fem_T1,~)

%% Initialize: Input

% Select displ hscyl
% ITM
u_xyf_hscyl_itm_tr_z0 = getfield(displ.hscyl,'trench','itm',   [displ_dir load_dir '_xyf_hscyl_itm_z0']);
% Tunnel 1
u_xyf_hscyl_fem_all_tr_T1 = getfield(displ.hscyl,'trench','fem',   'T1',[displ_dir load_dir '_xyf_hscyl_fem_all']);
    
% Geo
y = geo_T1.y;

% Tunnel 1
r0_T1          = geo_T1.radius;
y_Tc_T1        = geo_T1.y_Tc;
nod_midline_T1 = dis_fem_T1.nod_midline;

% y-coordinate for different sampling in ITM & FEM
y_comb = [y(y<=-r0_T1+y_Tc_T1)  nod_midline_T1(:,2).'+y_Tc_T1  y(y>=r0_T1+y_Tc_T1)];

% displ tunnel/trench midline → line trough center of cyl
u_xyf_fem_midline_tr_T1 = u_xyf_hscyl_fem_all_tr_T1(:,nod_midline_T1(:,1),:,:);

    
%% Combination of the solution parts ITM-FEM-ITM - Displ
% ITM(x<r) - FEM nod_middline_in - ITM(x>r)
u_xyf_comb_tr_z0   = [ u_xyf_hscyl_itm_tr_z0(:,y<=-r0_T1+y_Tc_T1,:,:)   u_xyf_fem_midline_tr_T1(:,:,:,:)     u_xyf_hscyl_itm_tr_z0(:,y>=r0_T1+y_Tc_T1,:,:)];      % (x,   y,z=0,-Omega)


%% Output
geo_T1.y_comb = y_comb;

varname1 = [displ_dir load_dir '_xyf_comb_tr_z0'];
displ.hscyl.trench.itmfem.(varname1) = u_xyf_comb_tr_z0;

end



