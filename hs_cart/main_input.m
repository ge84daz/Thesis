%__________________________________________________________________________
%
%       main_input.m    
%
%       DESCRIPTION: 
%       Main input for homogenous halfspace in cartesian coordinates with rectangular load
%       Load case: harmonic
%
%       REMARK:
%       Original Author: Julian Freisinger (main_input)
%       Modified:  Tom Hicks              
%       Date:    23-05-2023 
%       Changed: 24-09-2024 - Hicks               
%__________________________________________________________________________

%% Explanations  

% Bx, By                Width of ground surface in x- and y-direction [m]
%                       Equivalent to repetition length for Fourier Transformation

% Nx_hs                 Number of Fourier members = Number of Points x-direction    
% Ny_hs                 Number of Fourier members = Number of Points y-direction
% Nt                    Number of Frequencies and Time steps

% f                     Frequency in [Hz] 
% freq                  Either one frequency or frequency interval [f1,f2,...,fn] -> Calculation for all of this frequencies
% T                     Period length 

% P0_itm_hs_1L          Amplitude of harmonic cosine load on halfspace surface

% bx_itm_hs_1L          Load width. Symmetric to block load midpoint. Half width of loading in x in [m] (must be a multiple of dx=Bx/Nx_hs)
% by_itm_hs_1L          Load width. Symmetric to block load midpoint. Half width of loading in x in [m] (must be a multiple of dy=By/Ny_hs) 

% pos_load_x_itm_hs_1L  Load position (mid point of rectangular block load) in x-y-plane as distance from origin in [m] -> (must be a multiple of dx resp. dy) including sign
% pos_load_y_itm_hs_1L  Load position (mid point of rectangular block load) in x-y-plane as distance from origin in [m] -> (must be a multiple of dx resp. dy) including sign
% pos_load_z_itm_hs_1L  Load position in z-direction -> on surface:    z0   -> on interlayer: zH

% calc.system           Hs          -> Only homog. hs solution for displ and flex     -> f_displ_hs
%__________________________________________________________________________



%% 1.) Initialize: 

clc
if exist('wait', 'var'),delete(wait), end
clear all
close all
profile clear 
profile on


%% 2.) Input

% Initialization
% Material parameters
calc.mat_hs         = 'stiff_v01';       % mat_hs

% Solving over kx
calc_cont.mode      = 'serial';      % 'serial', 'parallel'

% struct predef
mat.def  = 'Material Properties';

% for 2.5D or 3D approach
calc.kx        = 'kx=all';           % 'kx=all', 'kx=0'    

% Geometry and Discretization
dis_itm.Nx_hs       = 2^8;              % 2^10
dis_itm.Ny_hs       = 2^8;              % 2^10

% Frequency content:
dis_itm.freq     = [0,30];              
geo.Bx           = 128;                  % normally 128 - depends on frequency and material
geo.By           = 128;                  % normally 128

% for plot of wave amplitudes in kx-ky to avoid high memory when using Nx = 10^12
geo.x_eval_xy = 64;
geo.y_eval_xy = 64;



%% Loading
loading.type           = 'harmonic';    % 'harmonic' 
loading.Nf_int         = 2^16;
loading.T_int          = 4;

loading.P0_itm_hs1L    = 1.0;           % Load amplitude for TF
%loading.bx_itm_hs1L    = geo.Bx/2;     % 6*geo.Bx/dis_itm.Nx_hs; ; %geo.Bx/2; % half total load width
loading.bx_itm_hs1L    = 2.0;           % 6*geo.Bx/dis_itm.Nx_hs; ; %geo.Bx/2; % half total load width
loading.by_itm_hs1L    = 2.0;           % 6*geo.By/dis_itm.Ny_hs; ; %geo.By/2; % half total load width
loading.pos_x_itm_hs1L = 0;
loading.pos_y_itm_hs1L = 0;             % for example = 8.0m
loading.pos_z_itm_hs1L = 'z0';          %'z0'

% Specific load input (time history)
% Harmonic
loading.harmonic.f  = [30];  % [f1...fn]
loading.harmonic.P0 = [1]; % [p1...pn]
loading.harmonic.nT = 4;   % T=nT*T0 with T0=1/fmin
loading.harmonic.Nt = 2^8; % Time sample rate

            
%% 3.) Plot Control

%% 3.1.) Plot control
plot_cont.geo         = 'yes';
plot_cont.geo_found   = 'yes';
plot_cont.load        = 'yes';
plot_cont.displ_soil  = 'yes';
plot_cont.displ_found = 'yes';
plot_cont.export      = 'no';
plot_cont.GID         = 'no';
plot_cont.postproc    = 'no';
plot_cont.u_anim      = 'no';

plot_cont.lowMemory   = 'yes';

%% 3.2.) Plot parameters
plot_cont.Fontsize      = 12;
plot_cont.color_green   = [0.4660, 0.6740, 0.1880];  % green
plot_cont.color_cyan    = [0, 0.4470, 0.7410];       % blue
plot_cont.color_blue    = [0 0 1];                   % Dark blue -> hs_1L
plot_cont.color_orange  = [1 0.5 0.1];               % Orange    -> ItmFem
plot_cont.color_orange2 = [0.8500, 0.3250, 0.0980];  % Orange 2
plot_cont.color_yellow  = [0.9290, 0.6940, 0.1250];  % yellow
plot_cont.exp_width     = 16;                        % small: 8;   middle: 16;  large: 20 
plot_cont.exp_heigth    = 9.6;                       % small: 4.8; middle: 9.6; large: 16 
plot_cont.resolution    = 300;                       %[dpi]
plot_cont.fig_position  = [632,260,1403,1040];       % 4 plot
plot_cont.fig_position1 = [632,808,710,492];         % 1 plot
plot_cont.fig_position2 = [632,808,1403,492];        % 2 plots
plot_cont.leg           = 'yes';                     %'Yes', 'No'
plot_cont.loc_ledg      = 'northeast';
plot_cont.FaceColor     = 'interp';                  % [0 0 0], 'interp'
plot_cont.leg_box       = 'off';                     % 'off', 'on'
plot_cont.map           = 'No';                      % 'Yes', 'No' -> for gray resp. black/white plots
plot_cont.colormap      = 'gray';                    % 'autumn', 'summer', 'winter', 'jet', 'parula'
plot_cont.shading       = 'yes';
plot_cont.shading_typ   = 'flat';
plot_cont.closefigures  = 'No';                      % create figures but close directly after saving
plot_cont.interpreter   = 'latex';                   % 'latex', 'tex'

set(0,'defaultfigurecolor',[1 1 1])
plot_cont.exportfig     = ['Width',plot_cont.exp_width,'Height',plot_cont.exp_heigth,'Fontmode','fixed','Fontsize',plot_cont.Fontsize,'Color','cmyk','resolution',plot_cont.resolution];

% Set variables
switch calc_cont.mode
    case 'serial';   calc_cont.M  = 0;
    case 'parallel'; calc.parPool = gcp;
                     calc_cont.M  = calc.parPool.NumWorkers;
end


%% 5.) Run main routine
main_routine_hs


