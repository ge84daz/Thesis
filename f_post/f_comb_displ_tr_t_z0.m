%__________________________________________________________________________
%
%        f_comb_displ_tr_t_z0.m    
%
%       INPUT: 
%       - displ_dir     = Displ. direction
%       - load_dir      = Displ. due to load on load_dir
%       - displ         = Displacements
%       - dis_itm       = Discretization ITM
%       - dis_fem       = Discretization FEM
%       - geo           = Geometry
%       - loading       = Load (amplitudes, distribution,...)
%       - mat           = Material Parameters
%       - plot_cont     = plot control parameters
%       - path          = path parameters
%       - calc          = Calculation parameters
%       - fev           = Evaluation frequency
%       - nnPiF         = load position w.r.t. y (in case of flexibility)
%
%       OUTPUT:
%
%       DESCRIPTION: 
%       The function combines the solution parts on the halfspace surface of the
%       ITM and FEM System for a halfspace with a cylindrical trench/slit and
%       displayes them in different plots, also compared to a reference solution
%       calculated using a layered halfspace with one layer.
%       Furthermore the Amplitude Reduction Factors AR are calculated for
%       different setups.
%
%       - only displ. in frequency domain plotted for negative frequency f
%       - if harmonic time dependent displ shall be plotted:
%         Time of maximum displ t_max_itm_fem -> equal for trench, slit, hs1L
%       - For calc.system = 'slit'
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_ifft_t_displ_hscyl_harm)
%       Modified:   Tom Hicks                 
%       Date:       07-11-2019 
%       Changed:    04-12-2024 - Hicks   
%__________________________________________________________________________


function [displ,geo_T1] = f_comb_displ_tr_t_z0(displ_dir,load_dir,displ,~,geo_T1,dis_fem_T1,~)

% Initialize: Input
% Select displ hscyl

% ITM
u_xyt_hscyl_tr_z0                   = getfield(displ.hscyl,'trench','itm_Pt',[displ_dir load_dir '_xyt_hscyl_itm_z0']);
u_xyf_hscyl_itm_tr_z0               = getfield(displ.hscyl,'trench','itm_Pt',[displ_dir load_dir '_xyf_hscyl_itm_z0']);
u_xyf_hscyl_itm_tr_z0_TF_nint       = getfield(displ.hscyl,'trench','itm_Pt',[displ_dir load_dir '_xyf_hscyl_itm_z0_TF_nint']);

% Tunnel 1
u_xyt_hscyl_fem_all_tr_T1            = getfield(displ.hscyl,'trench','fem_Pt','T1',[displ_dir load_dir '_xyt_hscyl_fem_all']);
u_xyf_hscyl_fem_all_tr_T1            = getfield(displ.hscyl,'trench','fem_Pt','T1',[displ_dir load_dir '_xyf_hscyl_fem_all']);
u_xyf_hscyl_fem_all_tr_TF_nint_T1    = getfield(displ.hscyl,'trench','fem_Pt','T1',[displ_dir load_dir '_xyf_hscyl_fem_all_TF_nint']);

% Geo
y = geo_T1.y;

% Tunnel 1
r0_T1          = geo_T1.radius;
y_Tc_T1        = geo_T1.y_Tc;
nod_midline_T1 = dis_fem_T1.nod_midline;

% y-coordinate for different sampling in ITM & FEM
y_comb = [y(y<=-r0_T1+y_Tc_T1)  nod_midline_T1(:,2).'+y_Tc_T1  y(y>=r0_T1+y_Tc_T1)];

% displ tunnel/trench midline → line trough center of cyl
u_xyt_fem_midline_tr_T1         = u_xyt_hscyl_fem_all_tr_T1(:,nod_midline_T1(:,1),:,:);
u_xyf_fem_midline_tr_T1         = u_xyf_hscyl_fem_all_tr_T1(:,nod_midline_T1(:,1),:,:);
u_xyf_fem_midline_tr_TF_nint_T1 = u_xyf_hscyl_fem_all_tr_TF_nint_T1(:,nod_midline_T1(:,1),:,:);


%% Combination of the solution parts ITM-FEM-ITM - Displ
% ITM(x<r) - FEM nod_middline_in - ITM(x>r)
u_xyt_comb_tr_z0           = [ u_xyt_hscyl_tr_z0(:,y<=-r0_T1+y_Tc_T1,:,:)              u_xyt_fem_midline_tr_T1(:,:,:,:)          u_xyt_hscyl_tr_z0(:,y>=r0_T1+y_Tc_T1,:,:)];      % (x,y,z=0,t)
u_xyf_comb_tr_z0           = [ u_xyf_hscyl_itm_tr_z0(:,y<=-r0_T1+y_Tc_T1,:,:)          u_xyf_fem_midline_tr_T1(:,:,:,:)          u_xyf_hscyl_itm_tr_z0(:,y>=r0_T1+y_Tc_T1,:,:)];  % (x,y,z=0,f_tot)
u_xyf_comb_tr_z0_TF_nint   = [ u_xyf_hscyl_itm_tr_z0_TF_nint(:,y<=-r0_T1+y_Tc_T1,:,:)  u_xyf_fem_midline_tr_TF_nint_T1(:,:,:,:)  u_xyf_hscyl_itm_tr_z0_TF_nint(:,y>=r0_T1+y_Tc_T1,:,:)];  % (x,y,z=0,f_tot)

%% Output
geo_T1.y_comb = y_comb;

varname1 = [displ_dir load_dir '_xyt_comb_tr_z0'];
varname2 = [displ_dir load_dir '_xyf_comb_tr_z0'];
varname3 = [displ_dir load_dir '_xyf_comb_tr_z0_TF_nint'];

displ.hscyl.trench.itmfem_Pt.(varname1) = u_xyt_comb_tr_z0;
displ.hscyl.trench.itmfem_Pt.(varname2) = u_xyf_comb_tr_z0;
displ.hscyl.trench.itmfem_Pt.(varname3) = u_xyf_comb_tr_z0_TF_nint;

end



