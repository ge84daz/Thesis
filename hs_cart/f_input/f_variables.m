%__________________________________________________________________________
%
%       f_variables.m    
%
%       INPUT: 
%       dis_itm       = Discretisation
%       geo           = Geometry
%
%       OUTPUT: 
%       dis_itm       = Discretisation
%       geo           = Geometry
%
%       DESCRIPTION: 
%       Determine all secondary variables out of input variables.
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_variables_hs_1L)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________

function [dis_itm,geo,mat] = f_variables(geo,dis_itm,mat,loading,calc)
     
%% Initialize: Input

Bx          = geo.Bx;
By          = geo.By;

freq        = dis_itm.freq;
Nx_hs       = dis_itm.Nx_hs;
Ny_hs       = dis_itm.Ny_hs;


%% Discretization and geometry parameters
% Inkrements
dx   = Bx/Nx_hs;
dy   = By/Ny_hs;

% Space and time Variables
x    = -Bx/2:dx:Bx/2-dx;
y    = -By/2:dy:By/2-dy;
N_OF = Ny_hs*Nx_hs;          

f         = [-flip(freq)];                           % => [-fn ... -f1] with f1 <0  
f_tot     = [-flip(freq(freq >0)) 0 freq(freq > 0)]; % => [-fn ... -f1 0 f1 ... fn]
omega     = 2*pi*f;
omega_tot = 2*pi*f_tot;

Nf        = length(f);
Nf_tot    = length(f_tot);

% Variables in the wavenumber and frequency domain
kx0     = 2*pi/Bx;
kx      = -Nx_hs/2*kx0:kx0:(Nx_hs/2-1)*kx0;
ky0     = 2*pi/By;
ky      = -Ny_hs/2*ky0:ky0:(Ny_hs/2-1)*ky0;



%% Output
% Discretisation soil
dis_itm.dx  = dx;
dis_itm.dy  = dy;
dis_itm.kx  = kx;
dis_itm.kx0 = kx0;
dis_itm.ky  = ky;
dis_itm.ky0 = ky0;

dis_itm.f     = f;
dis_itm.omega = omega;
dis_itm.Nf    = Nf;

dis_itm.f_tot     = f_tot;
dis_itm.omega_tot = omega_tot; 
dis_itm.Nf_tot    = Nf_tot;

% Geometry
geo.x      = x;
geo.y      = y;
geo.N_OF   = N_OF;

