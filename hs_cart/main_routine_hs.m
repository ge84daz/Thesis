%__________________________________________________________________________
%
%       main_routine.m    
%
%       DESCRIPTION: 
%       Main routine for n soil layers over homogenous halfspace
%       Load case: harmonic
%
%       REMARK:
%       Original Author: Julian Freisinger (main_hs_1L_fvar)
%       Modified:  Tom Hicks                 
%       Date:    23-05-2023 
%       Changed: 24-09-2024 - Hicks  
%__________________________________________________________________________


%% 1.) Initialize
disp('A - Initialize')

% Create relative paths
path.current = pwd; cd ..  
path.submain = pwd; cd ..  
path.main    = pwd; cd(path.current)

% Path Object
calc.folder                 = ['f_all_' loading.type '_' calc.kx];

path.functions_input        = fullfile(path.current, 'f_input');
path.functions_load         = fullfile(path.current, 'f_load');
path.functions_displ        = fullfile(path.current, 'f_displ');
path.functions_ifft         = fullfile(path.current, 'f_ifft');
path.functions_plot         = fullfile(path.current, 'f_plot');
path.functions_post         = fullfile(path.current, 'f_post');

path.results                = fullfile(path.current, ['\results_' calc.folder]);
path.latex                  = fullfile(path.results, 'latex');
path.GiD                    = fullfile(path.results, 'GiD');
path.input                  = fullfile(path.results, 'input');
path.figures                = fullfile(path.results, 'figures');
path.mat                    = fullfile(path.results, 'mat');
path.graphics               = ['results_' calc.folder '/figures/'];


addpath(path.functions_input);
addpath(path.functions_load);
addpath(path.functions_displ);
addpath(path.functions_ifft);
addpath(path.functions_plot);
addpath(path.functions_post);

mkdir(path.results); 
mkdir(path.GiD); 
mkdir(path.input);
mkdir(path.figures); 
mkdir(path.mat); 

% Start time measure
tic


%% 2.) Material
[mat] = f_material(mat,calc);


%% 3.) Variables
% Determine all necessary variables out of the input variables
[dis_itm,geo,mat] = f_variables(geo,dis_itm,mat,loading,calc);

%% 4.) Loading for TF
disp('B - Load for TF')
% Load with specified spatial distribution and constant frequency spectrum 
% for determination of the system response to a harmonic excitation 
% for each frequency specified in dis_itm.freq in the input 
% -> Result: System transfer function (TF)

% for rectangular loading
[loading] = f_load_itm_tf_rect(path,calc,plot_cont,loading,dis_itm,geo);

% for comparison with hs_polar
% [loading] = f_load_itm_tf_gaus(path,calc,plot_cont,loading,dis_itm,geo);

%% 5.) Plot geo, load
f_plot_geo(path,calc,plot_cont,loading,geo);
tic

%% 6.) Displ for hs, hs1L, hs1L_fixed
disp('C - Determine displ. TF')
[displ] = f_displ_itm(dis_itm,geo,mat,loading,calc,calc_cont);
toc
%% 7.) IFFT displ (kx,ky) → (x,y) 
% Inverse Transformation of the displacements 
% ux,uy,uz due to load Pz
disp('D - IFFT (x,y) - uiPz')
[displ] = f_ifft_xy_displ_TF('ux','Pz',displ,dis_itm,geo,loading,calc); 
[displ] = f_ifft_xy_displ_TF('uy','Pz',displ,dis_itm,geo,loading,calc); 
[displ] = f_ifft_xy_displ_TF('uz','Pz',displ,dis_itm,geo,loading,calc); 

%% 8.) Transient load and response: P(t), u(x,y,z,f) and u(x,y,z,t) 
disp('E - Transient soil response u(t)')
% Transient (time dependent) course of load with spatial distribution
% specified above. The transient system response is finally obtained by
% multiplication of the transfer function of the system (TF), due to the 
% harmonic loads with constant frequency spectrum, and the spectrum of 
% the transient load -> u(x,y,z,f) = TF(x,y,z,f)*P(f).
% The time response u(x,y,t,t) is obtained via the IFFT.
[loading,dis_itm] = f_load_hs1L_itm_trans(path,plot_cont,loading,dis_itm);

% Determine u_z(x,y,t,t) due to Pz
[displ] = f_ifft_t_displ_Pt('uz','Pz',displ,dis_itm,geo,loading,mat,calc,plot_cont,path);

% Plot displ at soil surface z=0 for one frequency and tmax
f_plot_displ_hs_surf('uz','Pz',dis_itm,geo,displ,loading,plot_cont,path,calc)

clc
disp('Calculation finished successfully');
disp('Outputfiles exported successfully');










