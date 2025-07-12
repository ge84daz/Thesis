%__________________________________________________________________________
%
%       main_routine_hs_cyl.m    
%
%       DESCRIPTION: 
%       - Main routine for homogeneous halfspace including a tunnel
%       - non-equidistant wavenumber sampling included to save computational effort
%       - Load case: harmonic
%
%       REMARK:
%       Original Author: Julian Freisinger (main_calc_hs_cyl_1L)
%       Modified:  Tom Hicks                 
%       Date:    21-09-2023 
%       Changed: 02-12-2024 - Hicks  
%__________________________________________________________________________


%% 1.) Initialize
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
num = 1;disp([Alphabet(num) ' - Initialize Parameter']);num = num +1;
warning('off','all')

% Create relative paths
path.current = pwd; cd ..  
path.submain = pwd; cd ..  
path.main    = pwd; cd(path.current)

% Path Object
path.functions_input        = fullfile(path.current, 'f_input');
path.functions_load         = fullfile(path.current, 'f_load');
path.functions_mex          = fullfile(path.current, 'f_mex');
path.functions_displ        = fullfile(path.current, 'f_displ');
path.functions_ifft         = fullfile(path.current, 'f_ifft');
path.functions_plot         = fullfile(path.current, 'f_plot');
path.functions_post         = fullfile(path.current, 'f_post');
switch calc.system  
    case 'tunnel'
path.aoad_T1                = fullfile(path.current, 'ansys_wd_T1');
    case 'trench'
path.aoad_T1                = fullfile(path.current, 'ansys_wd_trench/ansys_wd_T1');
end

path.results                = fullfile(path.current, ['/results_' calc.folder]);
path.latex                  = fullfile(path.results, 'latex');
path.GiD                    = fullfile(path.results, 'GiD');
path.input                  = fullfile(path.results, 'input');
path.figures                = fullfile(path.results, 'figures');
path.mat                    = fullfile(path.results, 'mat');
path.ansys                  = fullfile(path.results,    'ansys');
path.graphics               = ['results_' calc.folder '/figures/'];

addpath(path.functions_input);
addpath(path.functions_load);
addpath(path.functions_mex);
addpath(path.functions_displ);
addpath(path.functions_ifft);
addpath(path.functions_plot);
addpath(path.functions_post);
mkdir(path.results); 
mkdir(path.GiD); 
mkdir(path.input);
mkdir(path.figures); 
mkdir(path.mat); 
mkdir(path.aoad_T1);
mkdir(path.ansys);

% Create new mex files 
switch calc_cont.mex
    case 'Create_new'; f_stiff_fem_script % not working so far
end

% Export variables for Ansys Mesh generation for tunnel T1
writematrix(geo_T1.radius, [path.aoad_T1 '/radius_export.txt'])
writematrix(dis_fem_T1.elem_size, [path.aoad_T1 '/elem_size_export.txt'])
writematrix(dis_fem_T1.div, [path.aoad_T1 '/div_export.txt'])
warning('on','all')


%% 2.) Material
disp([Alphabet(num) ' - Initialize Material']);num = num +1;
% Define material properties and determine secondary material parameters
[mat] = f_material(mat,calc);

%% 3.) Mesh, Import and Sort
disp([Alphabet(num) ' - Mesh, Import and sort FEM subsystem']);num = num +1;

% Mesh the FEM Structure T1
f_mesh(path.aoad_T1,calc_cont,calc)

% Import and Sorting of Nodes and Elements T1
[dis_fem_T1] = f_sort_cyl(geo_T1,path.aoad_T1,dis_fem_T1,calc_cont,calc);

%% 4.) Variables
disp([Alphabet(num) ' - Variables']);num = num +1;
% Determine all necessary variables out of the input variables
[dis_itm,dis_fem_T1,geo_T1,mat,trunc_T1] = f_variables_cyl(geo_T1,dis_itm,dis_fem_T1,mat,loading,calc,calc_cont);

%% 5.) Apply Loading
disp([Alphabet(num) ' - Loading']);num = num +1;
% Load with specified spatial distribution and constant frequency spectrum
% for determination of the system response to a harmonic excitation for each
% frequency specified in dis_itm.freq in the input
% -> Result: System transfer function (TF) for given P(x,y,f)

switch lower(loading.distr_xy) % till now only load applied on halfspace surface implemented
    case 'rect'
        [loading] = f_load_hscyl_itm_tf_rect(path,plot_cont,loading,dis_itm,geo_T1);                            % hs_cyl_itm
        [loading] = f_load_hscyl_fem_tf_rect(path,plot_cont,loading,dis_itm,dis_fem_T1,geo_T1,calc_cont,calc);  % hs_cyl_fem
        
    case 'line'
        loading.bx_itm = geo_T1.Bx/2; 
        [loading] = f_load_hscyl_itm_tf_rect(path,plot_cont,loading,dis_itm,geo_T1);                            % hs_cyl_itm
        [loading] = f_load_hscyl_fem_tf_line(path,plot_cont,loading,dis_itm,dis_fem_T1,geo_T1,calc_cont,calc);  % hs_cyl_fem
        
    case 'point'
        [loading] = f_load_hscyl_itm_tf_point(path,plot_cont,calc_cont,loading,dis_itm,geo_T1);                 % hs_cyl_itm
        [loading] = f_load_hscyl_fem_tf_point(path,plot_cont,calc_cont,loading,dis_itm,dis_fem_T1,geo_T1);      % hs_cyl_fem
end

%% 6.) Plot Geometry 
disp([Alphabet(num) ' - Plot geometry']);num = num +1;
f_plot_geo_hs_cyl(geo_T1,dis_fem_T1,loading,path,plot_cont,calc); 

%% 7.) Determine displ TF ITM-FEM 
disp([Alphabet(num) ' - Determine displ. TF Hscyl']);num = num +1;
% Start time measure for developing K_GES and solving for u_GES
tic

[displ,dis_fem_T1] = f_displ_hscyl_tot(dis_itm,dis_fem_T1,geo_T1,loading,mat,displ,trunc_T1,calc_cont,calc);

% end time measure
dur = seconds(toc);
sec = seconds(dur);
min = minutes(dur);
calc.duration_min = [num2str(min) ' min'];
calc.duration_sec = [num2str(sec) ' sec'];

%% 8.) IFFT of displ ITM-FEM system
disp([Alphabet(num) ' - IFFT Hscyl']);num = num +1;

% IFFT of the displ. on the hs surface z=0 -> TF
[displ] = f_ifft_xy_displ_hscyl_itm_z0('ux','Pz',plot_cont,displ,dis_itm,geo_T1,calc_cont,calc);
[displ] = f_ifft_xy_displ_hscyl_itm_z0('uy','Pz',plot_cont,displ,dis_itm,geo_T1,calc_cont,calc);
[displ] = f_ifft_xy_displ_hscyl_itm_z0('uz','Pz',plot_cont,displ,dis_itm,geo_T1,calc_cont,calc);

% IFFT of the displ. on cyl. boundary and within FEM subsystem for ux,uy,uz -> TF
[displ] = f_ifft_xy_displ_hscyl_fem('Pz','T1',calc,calc_cont,plot_cont,displ,dis_itm,dis_fem_T1,geo_T1);

%% 9.) Combine ITM and FEM results (x,y,f,nnPiF)
% → for plots, flex and thus postprocessing
if strcmpi(calc.system,'trench')
    [displ,geo_T1] = f_comb_displ_tr_f_z0('ux','Pz',displ,geo_T1,dis_itm,dis_fem_T1,calc);
    [displ,geo_T1] = f_comb_displ_tr_f_z0('uy','Pz',displ,geo_T1,dis_itm,dis_fem_T1,calc);
    [displ,geo_T1] = f_comb_displ_tr_f_z0('uz','Pz',displ,geo_T1,dis_itm,dis_fem_T1,calc);
end

%% 10.) Evaluation for f=fev: Plots in frequency space
disp([Alphabet(num) ' - Plots']);num = num +1;
% Combination of Results and Plots
switch calc.system
    case 'tunnel'
        [geo_T1] = f_plot_displ_hscyl_f_tun('uz','Pz',displ,geo_T1,dis_itm,dis_fem_T1,mat,plot_cont,path,calc);
     
    case 'trench'
        [displ,geo_T1] = f_plot_displ_hscyl_f_tr_sl('uz','Pz',displ,geo_T1,dis_itm,dis_fem_T1,loading,mat,plot_cont,path,calc);
end


%% 11.) Transient response ITM/FEM - transient load P(t)
disp([Alphabet(num) ' - Transien loading and plots']);num = num +1;
% → u(x,y,z,f) and u(x,y,z,t) 
if strcmpi(calc.antype,'stationary') && any(strcmpi(calc.system,{'tunnel','trench'}))
    disp([Alphabet(num) ' - Displ hscyl SOIL - transient']);num = num +1;
    % Transient (time dependent) course of load with spatial distribution
    % specified above. The transient system response is finally obtained by
    % multiplication of the transfer function of the system (TF), due to the 
    % harmonic loads with constant frequency spectrum, and the spectrum of 
    % the transient load -> u(x,y,z,f) = TF(x,y,z,f)*P(f).
    % The time response u(x,y,t,t) is obtained via the IFFT.
    [loading,dis_itm] = f_load_hscyl_t_xy_harm( loading,dis_itm,calc_cont,plot_cont,path);
    
    % Determine u_i(x,y,t,t) due to Pi
    [displ] = f_ifft_t_displ_hscyl_harm('ux','Pz',displ,dis_itm,dis_fem_T1,geo_T1,loading,calc,calc_cont);
    [displ] = f_ifft_t_displ_hscyl_harm('uy','Pz',displ,dis_itm,dis_fem_T1,geo_T1,loading,calc,calc_cont);
    [displ] = f_ifft_t_displ_hscyl_harm('uz','Pz',displ,dis_itm,dis_fem_T1,geo_T1,loading,calc,calc_cont);
    
    % Combine ITM and FEM results
    if strcmpi(calc.system,'trench')
        [displ,geo_T1] = f_comb_displ_tr_t_z0('ux','Pz',displ,loading,geo_T1,dis_fem_T1,calc);
        [displ,geo_T1] = f_comb_displ_tr_t_z0('uy','Pz',displ,loading,geo_T1,dis_fem_T1,calc);
        [displ,geo_T1] = f_comb_displ_tr_t_z0('uz','Pz',displ,loading,geo_T1,dis_fem_T1,calc);
    end
    
    % Plots
    f_plot_displ_hscyl_t_harm( 'uz','Pz','ini',displ,dis_itm,geo_T1,loading,calc,plot_cont,path)
end



%% 12.) Save files and variables
disp([Alphabet(num) ' - Save variables']);num = num +1;
% savefile = [path.mat '\dis_itm.mat'];     save(savefile,'dis_itm');
% savefile = [path.mat '\dis_fem_T1.mat'];  save(savefile,'dis_fem_T1');
% savefile = [path.mat '\geo_T1.mat'];      save(savefile,'geo_T1');
% savefile = [path.mat '\trunc_T1.mat'];    save(savefile,'trunc_T1');
% savefile = [path.mat '\mat.mat'];         save(savefile,'mat');
% savefile = [path.mat '\path.mat'];        save(savefile,'path');
% savefile = [path.mat '\plot_cont.mat'];   save(savefile,'plot_cont');
% savefile = [path.mat '\loading.mat'];     save(savefile,'loading');
% savefile = [path.mat '\calc.mat'];        save(savefile,'calc');
% savefile = [path.mat '\calc_cont.mat'];   save(savefile,'calc_cont');
% clear savefile

%% 13:) GiD export 
disp([Alphabet(num) ' - Export GiD']);num = num +1;
if strcmpi(plot_cont.gid, 'yes')
% Export results hscyl 
 % 2D   -> cutting plane at x=0
 f_GiD_export_2D('Pz','hscyl_res_2D',dis_itm,dis_fem_T1,geo_T1,loading,displ,path,calc,calc_cont)    
 % 3D -> threedimensional model
 f_GiD_export_3D('Pz','hscyl_res_3D',dis_itm,dis_fem_T1,geo_T1,loading,displ,path,calc,calc_cont)    
end


%% 14.) Finalize
disp([Alphabet(num) ' - Finalize']);num = num +1;
disp('Calculation finished successfully');
disp('Outputfiles exported successfully');
disp(['duration = ' calc.duration_min]);
disp(['duration = ' calc.duration_sec]);











