%__________________________________________________________________________
%
%       main_input_hs_cyl.m    
%
%       DESCRIPTION: 
%       - Main input for n soil layers over homogenous halfspace including length-invariant tunnel
%       - Wavenumber space is nonequidistant sampled to safe computational effort
%       - Load case: harmonic
%
%       REMARK:
%       Original Author: Julian Freisinger (main_input_hs_cyl_1L)
%       Modified:  Tom Hicks              
%       Date:    21-09-2023 
%       Changed: 21-09-2023 - Hicks               
%__________________________________________________________________________

%% Explanations  
% calc.system           Hs_cyl      -> Homogeneous halfspace including a length-invariant tunnel      

% Bx, By                Width of ground surface in x- and y-direction [m]
%                       Equivalent to repetition length for Fourier Transformation

% Nx_hs                 Number of Fourier members = Number of Points x-direction    
% Ny_hs                 Number of Fourier members = Number of Points y-direction
% Nt                    Number of Frequencies and Time steps
% f                     Frequency in [Hz] 
% T                     Period length 

% radius                Tunnel radius
% H                     Depth of tunnel center on z-axis

% div                   Number of elements per half of radius
% Nphi_t                Number of Fourier members along circumference phi on Tunnel surface
%                       Calculated in d_variables (=Number of boundary Nodes of mesh)

% P0_itm_hs             Aplitude of harmonic cosine load on half space surface z=0 on ITM subsystem [N/m]
% bx_itm                Load width. Symmetric to block load midpoint. Half width of loading in x in [m] (must be a multiple of dx=Bx/Nx_hs) 
% by_itm                Load width. Symmetric to block load midpoint. Half width of loading in x in [m] (must be a multiple of dy=By/Ny_hs)                                  
% pos_load_x_itm        Load position (mid point of rectangular block load) in x-y-plane as distance from origin in [m] -> (must be a multiple of dx resp. dy) including sign
% pos_load_y_itm        Load position (mid point of rectangular block load) in x-y-plane as distance from origin in [m] -> (must be a multiple of dx resp. dy) including sign

% P0_FE_Point           Amplitude of harmonic cosine load within FEM subsystem 
% P0_FE_Line            Amplitude of harmonic cosine load within FEM subsystem
% P0_FE_Rect            Amplitude of harmonic cosine load within FEM subsystem
% bx_fem                Load width. Symmetric to block load midpoint. Half width of loading in x in [m] (must be a multiple of dx=Bx/Nx_hs)
% by_fem                Load width. Symmetric to block load midpoint. Half width of loading in y in [m] (must be a multiple of dy_FE=radius/(2*div)) 
% by_fem_line           Load width. Symmetric to Line load midpoint.  Half width of loading in y in [m] (must be a multiple of dy_FE=radius/(2*div)) 
%                       Minimum 3 loaded nodes -> by_fem_line = 1*dy_FE => dy_FE = geo.radius/(2*dis_fem.div)
%                       For load of full foundation (case:within soil) -> geo.b_found/2
% pos_load_x_fem        Load position (mid point of rectangular block load) in x-y-plane as distance from origin in [m] -> (must be a multiple of dx resp. dy_FE) including sign
% pos_load_y_fem        Load position (mid point of rectangular block load) in x-y-plane as distance from origin in [m] -> (must be a multiple of dx resp. dy_FE) including sign
%                       Must be positioned inside the FEM domain -> pos_load_fem < radius    
% pos_load_z_fem        Position of midpoint of FEM lineload in z direction.
%                       'z0'             -> Nodes on midline of cyl. cavity
% load_dir              Direction of load: x,y,z
%__________________________________________________________________________

%% 1.) Initialize: 

clc
if exist('wait', 'var'),delete(wait), end
clear all
close all
profile clear 
profile on


%% 2.) Input
clc
if exist('wait', 'var'),delete(wait), end
clear all
close all
profile clear 
profile on

%% 2.1.) Input Variables
% 2.1.1.) Calc parameters

calc.folder             = 'test_thesis';          %
calc.system             = 'tunnel';               % tunnel, trench (halftunnel at z=0)
calc.structure          = 'tunnel';               % see documentation Diss JF

calc.mat_hs             = 'stiff_v01';            % mat_hs
calc.mat_fem_soil_T1    = 'stiff_v01';            % fem soil tunnel/trench T1

calc.antype             = 'stationary';           % 'stationary' - moving not implemented so far

calc_cont.mesh          = 'load';                  % 'new'  'load'
calc_cont.mode          = 'serial';               % 'serial', 'parallel'
calc_cont.spline        = 'without';              % 'with' 'without'
calc_cont.kx            = 'kx=all';               % 'kx=all', 'kx=0'
calc_cont.mex           = 'Use_current';          % 'Create_new', 'Use_current';
calc_cont.ansys_lic     = 'ansys';                %  
calc_cont.postproc      = 'no';

% 2.1.2.) Geometry and Discretization                                                  
dis_itm.Nx_hs           = 2^7; 
dis_itm.Ny_hs           = 2^7; 
dis_itm.freq            = [0,30];                 % One or several frequencies 

dis_itm.fev             = dis_itm.freq(2);        % ??
dis_fem_T1.div          = 4;                      % 4 for 32 nodes, 8 for 64 nodes
dis_fem_T1.elem_size    = 0.15;                   % coarse: 0.5;  fine: 0.15


geo_T1.Bx               = 64;
geo_T1.By               = 64;
geo_T1.x_eval_xy        = geo_T1.Bx/2;
geo_T1.y_eval_xy        = geo_T1.By/2;

% Standard case: One Tunnel T1 -> y_Tc = 0
geo_T1.radius          = 2;                     % radius of tunnel
geo_T1.H               = 5;                     % z-coord of tunnel center within halfspace without add. layer
geo_T1.y_Tc            = 0;                     % y-coord of tunnel center    

% 2.1.3.) Loading distr xy, ampl, sys

loading.distr_xy       = 'rect';       % 'rect','point','line'
loading.type           = 'harmonic';   % harmonic

loading.P0_itm_hs       = 1.0;
loading.bx_itm          = 2.0;
loading.by_itm          = 2.0;
loading.pos_load_x_itm  = 0.0;
loading.pos_load_y_itm  = 0.0;

loading.P0_FE           = 0.0;
loading.bx_fem          = 2.0;
loading.by_fem          = 2.0;   
loading.bz_fem          = 0.0;   
loading.pos_load_x_fem  = 0.0;
loading.pos_load_y_fem  = 0.0;     
loading.pos_load_z_fem  = 'z0';  % 'z0'
loading.load_dir        = 'z';   % direction of load: 'x','y','z'
loading.P_FE_T1         = 'no';  % 'yes','no' - load 1st tunnel 

% 2.1.4.) Loading transient - till now only harmonic included
% Discretization for interpolation over f
loading.Nf_int         = 2^16;
loading.T_int          = 4;
% Harmonic
loading.harmonic.f  = [dis_itm.fev];    % [f1...fn]
loading.harmonic.P0 = [1];              % [p1...pn]
loading.harmonic.nT = 4;                % T=nT*T0 with T0=1/fmin
loading.harmonic.Nt = 2^8;              % Time sample rate
loading.harm  = {'harmonic'};


%% 2.2) Plot control

% Plot control
plot_cont.geo             = 'yes';
plot_cont.geo_found       = 'yes';
plot_cont.load            = 'yes';
plot_cont.load_video      = 'no';
plot_cont.stresses        = 'no';
plot_cont.displ_trans     = 'yes';
plot_cont.displ_mov       = 'no';
plot_cont.displ_comp      = 'no';
plot_cont.displ_video     = 'no';
plot_cont.export          = 'no';
plot_cont.gid             = 'no';
plot_cont.norm_lambda_R   = 'no';
plot_cont.mesh            = 'no';
plot_cont.u_anim          = 'no';
plot_cont.closefigures    = 'no';  % create figures but close directly after saving
plot_cont.postproc        = 'no';

plot_cont.Fontsize        = 10;
plot_cont.LineWidth       = 1.3;
plot_cont.MarkerSize      = 5;
plot_cont.color_green     = [0.4660, 0.6740, 0.1880];  % green
plot_cont.color_cyan      = [0, 0.4470, 0.7410];       % blue
plot_cont.color_blue      = [0 0 1];                   % Dark blue -> hs_1L
plot_cont.color_orange    = [1 0.5 0.1];               % Orange    -> ItmFem
plot_cont.color_orange2   = [0.8500, 0.3250, 0.0980];  % Orange 2
plot_cont.color_yellow    = [0.9290, 0.6940, 0.1250];  % yellow
plot_cont.color_darkgray  = [0.50,0.50,0.50];
plot_cont.color_blue      = [0.00,0.45,0.74];
plot_cont.exp_width       = 16;                        % small: 8;   middle: 16;  large: 20  alt: 14
plot_cont.exp_heigth      = 9.6;                       % small: 4.8; middle: 9.6; large: 16  alt: 8.4
plot_cont.fig_position    = [632,260,1403,1040];       % [632,260,1403,1040] ; [381,128,1034,802]
plot_cont.fig_position_1  = [632,799,702,501];         % [632,260,1403,1040] [381,529,1034,401]; one subplot632,799,702,501
plot_cont.fig_position_2  = [632,799,1403,501];        % [632,260,1403,1040] [381,529,1034,401]; two subplots
plot_cont.exp_width       = 16;                        % small: 8;   middle: 16;  large: 20  alt: 14
plot_cont.exp_heigth      = 9.6;                       % small: 4.8; middle: 9.6; large: 16  alt: 8.4
plot_cont.resolution      = 300;                       %[dpi]
plot_cont.limval_xy       = 6;                         % n*lambda_R
plot_cont.limval_AR       = 1.5;
plot_cont.limval_IL_pos   = 20;
plot_cont.limval_IL_neg   = -5;
plot_cont.leg             = 'yes';                     %'yes', 'no'
plot_cont.loc_ledg        = 'northeast';
plot_cont.FaceColor       = 'interp';                  % [0 0 0], 'interp'
plot_cont.leg_box         = 'off';                     % 'off', 'on'
plot_cont.map             = 'no';                      % 'Yes', 'No' -> for gray resp. black/white plots
plot_cont.colormap        = 'gray';                    %'autumn', 'summer', 'winter', 'jet', 'parula'
plot_cont.shading         = 'yes';
plot_cont.shading_typ     = 'interp';
plot_cont.interpreter     = 'latex';                   %'latex', 'tex'
plot_cont.exportfig       = {'Width',plot_cont.exp_width,'Height',plot_cont.exp_heigth,'Fontmode','fixed','Fontsize',plot_cont.Fontsize,'Color','cmyk','resolution',plot_cont.resolution};

% Plot names
plot_cont.name_soft       = 'soft barrier';
plot_cont.name_stiff      = 'stiff barrier';
plot_cont.name_open       = 'open trench';
plot_cont.name_homo       = 'homogeneous';

% Set figure properties
matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','auto')
set(0,'defaultfigurecolor',[1 1 1])

% struct predef
mat.def           = 'Material Properties';
displ.def         = 'Displacements';
displ_post.def    = 'Displacement from Postprocessing SSI';
loading_post.def     = 'Loading for postprocessing';
loading_post.P_FE_T1 = loading.P_FE_T1;


%% 2.3) Select FEM mesh
% Select FEM mesh - only tunnel and trench implemented so far (see diss_jf_code)
switch calc.system
    case {'tunnel'}
        switch calc.structure
            case {'tunnel'}
                calc.mesh_ansys     = ['FEmesh/' 'FE_full_cyl'];
                calc_cont.mesh_type = 'full_cyl';   % no DOF skip
        end
    case {'trench'}
        switch calc.structure
            case 'trench'
                calc.mesh_ansys     = ['FEmesh/' 'FE_half_cyl'];
                calc_cont.mesh_type = 'half_cyl';
        end
end


%% 2.4) Set various variables
% Set variables for parallelization
switch calc_cont.mode
    case 'serial';   calc_cont.M=0;
    case 'parallel'; calc.parPool=gcp;
                     calc_cont.M=calc.parPool.NumWorkers;
end

%% 2.5) Input Control
% Depth of tunnel center on z-axis
switch calc.system
    case {'tunnel'}
        if geo_T1.H == 0;                     error('ERROR: Depth of Tunnel midpoint can no be 0 for calc.system = tunnel');
        elseif abs(geo_T1.H) < geo_T1.radius; error('ERROR: Depth of Tunnel midpoint must be bigger than radius for calc.system = tunnel');   end
    case {'trench'}
        if geo_T1.H ~= 0;                     error('ERROR: H must be 0 in case of trench or slit'); end    
end

% Itm load in case of calc.system trench
% -> no ITM load within trench
if loading.P0_itm_hs > 0
    switch calc.system
        case {'trench'}
            if loading.pos_load_y_itm + loading.by_itm < geo_T1.radius + geo_T1.y_Tc && loading.pos_load_y_itm - loading.by_itm > -geo_T1.radius + geo_T1.y_Tc
                error('ERROR: ITM load on hs surf must be located outside of the FEM subsystem (trench)!!');
            end
        otherwise
    end
end
   
%% 3.) Run main routine
main_routine_hs_cyl
    


