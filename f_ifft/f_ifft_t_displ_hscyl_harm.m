%__________________________________________________________________________
%
%        f_ifft_t_displ_hscyl_harm.m    
%
%       INPUT: 
%       - displ_dir     = cur. eval. displ dir
%       - load_dir      = cur. eval. load dir
%       - displ         = displ. from ITM-FEM calc
%       - dis_itm       = discr. ITM
%       - dis_fem       = discr. FEM
%       - geo           = system geom
%       - loading       = loading
%       - mat           = material param
%       - calc          = calculation param
%       - plot_cont     = plot control param
%       - calc_cont     = calculation control param
%       - path          = pathes
%
%       OUTPUT:
%       - u_xyf_hscyl_itm_z0    = freq spectrum of displ u(f) = TF(f)*P(f)
%       - u_kxkyf_hscyl_itm_z0  = freq spectrum of displ u(f) = TF(f)*P(f)
%       - u_xyf_hscyl_itm_z0    = time resp of displ u(t)
% 
%       all with dimensions u(x,y,f,nnPiF) resp. u(x,y,t,nnPiF)
%
%       DESCRIPTION: 
%       - Goal: determination of harmonic response of hscyl
%       - Apply excitation spectrum P(f) of harmonic load on TF
%           - possibly interpolate freq spectrum of P(f),T(f),u(f)
%           - determine displ u(f) = TF(f)*P(f)
%       - IFFT: u(f) -> u(t)
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_ifft_t_displ_hscyl_harm)
%       Modified:   Tom Hicks                 
%       Date:       25-09-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________

function [displ] = f_ifft_t_displ_hscyl_harm(displ_dir,load_dir,displ,dis_itm,dis_fem,geo,loading,calc,calc_cont)

% Initialize: Input

Ny_hs     = dis_itm.Ny_hs;
f_tot     = dis_itm.f_tot;
Nf        = dis_itm.Nf;
Nt_harm   = dis_itm.Nt_harm;

Pf_tot     = loading.Pf_tot;
loc_f_tot  = loading.loc_f_tot;    % location of f_calc in f_tot
loc_f_harm = loading.loc_f_harm;   % loaction of f_calc in f_harm

x_eval_xy = geo.x_eval_xy;
y_eval_xy = geo.y_eval_xy;
x         = geo.x;
y         = geo.y;
xev       = geo.x(abs(geo.x)<=x_eval_xy);
yev       = geo.y(abs(geo.y)<=y_eval_xy);
num_nod     = dis_fem.num_nod;

% Select displ hscyl u(x,y,f,nnPiF)
u_kxkyf_hscyl_tot       = getfield(displ.hscyl,calc.system,'tot',['u' load_dir '_kxkyf_hscyl_tot']);
u_xyf_hscyl_itm_z0      = getfield(displ.hscyl,calc.system,'itm',[displ_dir load_dir '_xyf_hscyl_itm_z0']);
u_xyf_hscyl_fem_all_T1  = getfield(displ.hscyl,calc.system,'fem','T1',[displ_dir load_dir '_xyf_hscyl_fem_all']);

% non interpolated -> f=f_tot
u_kxkyf_hscyl_itm_z0_TF = u_kxkyf_hscyl_tot(:,3:3:3*Ny_hs,:,:);
u_xyf_hscyl_itm_z0_TF   = u_xyf_hscyl_itm_z0(abs(x)<=x_eval_xy,abs(y)<=y_eval_xy,:,:);

% Tunnel 1
u_xyf_hscyl_fem_all_TF_T1 = u_xyf_hscyl_fem_all_T1(abs(x)<=x_eval_xy,:,:,:);


%% Apply excitation spectrum P(f) of harmonic load on TF
% Set Im(u_xyf_z0(f=0))=0
u_xyf_hscyl_itm_z0_TF(:,:,f_tot==0,:)     = real(u_xyf_hscyl_itm_z0_TF(:,:,f_tot==0,:));
% Tunnel 1
u_xyf_hscyl_fem_all_TF_T1(:,:,f_tot==0,:) = real(u_xyf_hscyl_fem_all_TF_T1(:,:,f_tot==0,:));

%  -> u(x,y,z,f) = TF(x,y,z,f).*P(f) by
%  multiplication of P(f) in 3rd dim on all displ for each (x,y)

% u(x,y,f,nnPiF)
u_kxkyf_hscyl_itm_z0_Pf_nint = u_kxkyf_hscyl_itm_z0_TF.*Pf_tot(1,1,1:Nf); % with f = f_tot
u_xyf_hscyl_itm_z0_Pf_nint   = u_xyf_hscyl_itm_z0_TF.*Pf_tot;         % with f = f_tot

% Tunnel 1
u_xyf_hscyl_fem_all_Pf_nint_T1 = u_xyf_hscyl_fem_all_TF_T1.*Pf_tot;    % with f = f_tot

%% IFFT f -> t
% Frequency/time response u(f), u(t)
dim4 = 1;
u_xyf_hscyl_itm_z0_Pf_harm = zeros(length(xev),length(yev),Nt_harm,dim4);

% Tunnel 1
u_xyf_hscyl_fem_all_Pf_harm_T1 = zeros(length(xev),num_nod,Nt_harm,dim4);

% sort results for f_calc in uf with higher resolution for IFFT
u_kxkyf_hscyl_itm_z0_Pf_harm = u_kxkyf_hscyl_itm_z0_Pf_nint;
u_xyf_hscyl_itm_z0_Pf_harm(:,:,loc_f_harm,:) = u_xyf_hscyl_itm_z0_Pf_nint(:,:,loc_f_tot,:);

% Tunnel 1
u_xyf_hscyl_fem_all_Pf_harm_T1(:,:,loc_f_harm,:) = u_xyf_hscyl_fem_all_Pf_nint_T1(:,:,loc_f_tot,:);

% IFFT
%  Multipy with Nt_harm as Pf not from forward FFT from
%  Pt, but here directly introduced in f domain
u_xyt_hscyl_itm_z0_Pf_harm     = Nt_harm*fftshift(ifft(ifftshift(u_xyf_hscyl_itm_z0_Pf_harm),[],3));

% Tunnel 1
u_xyt_hscyl_fem_all_Pf_harm_T1 = Nt_harm*fftshift(ifft(ifftshift(u_xyf_hscyl_fem_all_Pf_harm_T1),[],3));



%% Output

varname1 = [displ_dir load_dir '_kxkyf_hscyl_itm_z0'];
varname2 = [displ_dir load_dir '_xyf_hscyl_itm_z0'];
varname3 = [displ_dir load_dir '_xyt_hscyl_itm_z0'];
varname4 = [displ_dir load_dir '_xyf_hscyl_itm_z0_TF_nint'];
varname5 = [displ_dir load_dir '_xyf_hscyl_fem_all_TF_nint'];
varname6 = [displ_dir load_dir '_xyf_hscyl_fem_all'];
varname7 = [displ_dir load_dir '_xyt_hscyl_fem_all'];

displ.hscyl.(calc.system).itm_Pt.(varname1) = u_kxkyf_hscyl_itm_z0_Pf_harm;
displ.hscyl.(calc.system).itm_Pt.(varname2) = u_xyf_hscyl_itm_z0_Pf_harm;
displ.hscyl.(calc.system).itm_Pt.(varname3) = u_xyt_hscyl_itm_z0_Pf_harm;
displ.hscyl.(calc.system).itm_Pt.(varname4) = u_xyf_hscyl_itm_z0_TF;

% Tunnel 1
displ.hscyl.(calc.system).fem_Pt.T1.(varname5) = u_xyf_hscyl_fem_all_TF_T1;
displ.hscyl.(calc.system).fem_Pt.T1.(varname6) = u_xyf_hscyl_fem_all_Pf_harm_T1;
displ.hscyl.(calc.system).fem_Pt.T1.(varname7) = u_xyt_hscyl_fem_all_Pf_harm_T1;


end



