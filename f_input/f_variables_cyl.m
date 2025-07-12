%__________________________________________________________________________
%
%       f_variables_cyl.m    
%
%       INPUT: 
%       - dis_itm       = Discretization ITM
%       - dis_fem       = Discretization FEM
%       - geo           = Geometry
%       - calc_cont     = Control parameters
%       - radius        = Tunnel radius
%       - H             = Depth of tunnel center on z-axis
%       - Bx, By        = Width of ground surface in x- and y direction [m]
%       - T             = Time Period length
%       - Nx_hs, Ny_hs  = Number of Fourier members on ground surface 
%       - Nt            = Number of Frequencies or time steps 
%       - num_nod_bd    = Number of nodes on tunnel surface
%       - path_results  = Path for saving in Results folder
%       - nod_midline   = Nodes on midline in y- direction in fem sorted from -radius+(radius/2/div) ... radius-(radius/2/div)
%       - div           = subdivision of FE per half radius
%
%       OUTPUT: 
%       - dis_itm       = Discretization ITM
%       - dis_fem       = Discretization FEM
%       - geo           = Geometry
%       - trunc         = Truncation criteria 
%       - x,y,t                          = Space and time variables
%       - dx, dy, dphi                   = Inkrements (hs, tunnel)
%       - Nphi_t                         = Number of Points/Fouriermembers on tunnel surface
%       - kx,kx0,ky0,ky                  = Inkrements in wavenumber domain
%       - omega                          = Frequency
%       - delta_wh,delta                 = Truncation criteria (wh = repetition section)
%       - abstand_wh, abstand            = Distance to surface from tunnel center
% 
%       - N_OF, y_OF, r_OF, phi_OF,      = Geometry Parameter on hs surface
%       - y_OF_ges, r_OF_li, phi_OF_li, 
%       - r_OF_re, phi_OF_re
%
%       - N_LR, phi_LR, y_LR, z_LR,      = Geometry Parameter on tunnel surface
%       - phi_LR_ges, r_LR_li, phi_LR_li, 
%       - r_LR_re, phi_LR_re
%       - geo                            = Tunnel or Trench
%
%       DESCRIPTION: 
%       - Determine soil and tunnel geometry and discretization
%       - Determine trunction critera
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_variables_cyl)
%       Modified:   Tom Hicks                 
%       Date:       05-11-2019 
%       Changed:    21-09-2023 - Hicks   
%__________________________________________________________________________

function [dis_itm,dis_fem,geo,mat,trunc] = f_variables_cyl(geo,dis_itm,dis_fem,mat,loading,calc,calc_cont)
     
% Initialize: Input
Bx          = geo.Bx;
By          = geo.By;
H           = geo.H;
radius      = geo.radius;
y_Tc        = geo.y_Tc;     % y-coord. of tunnel center

Nx_hs       = dis_itm.Nx_hs;
Ny_hs       = dis_itm.Ny_hs;

freq        = dis_itm.freq;

num_nod_bd  = dis_fem.num_nod_bd;
div         = dis_fem.div;


%% Discretization and geometry parameters

% Number of Fourier members on cylindrical tunnel surface = Number of FE Nodes on boundary
Nphi_t = num_nod_bd;

% Inkrements
dx    = Bx/Nx_hs;
dy    = By/Ny_hs;
dphi  = 2*pi/Nphi_t;
dy_FE = radius/2/div;

% Space and time Variables
% phi already defined in f_sort
x = -Bx/2:dx:Bx/2-dx;
y = -By/2:dy:By/2-dy;

f         = -flip(freq);                             % => [-fn ... -f1] with f1 <0  
f_tot     = [-flip(freq(freq >0)) 0 freq(freq > 0)]; % => [-fn ... -f1 0 f1 ... fn]
omega     = 2*pi*f;
omega_tot = 2*pi*f_tot;

Nf        = length(f);
Nf_tot    = length(f_tot);

% Variables in the wavenumber and frequency domain
kx0    = 2*pi/Bx;
kx     = -Nx_hs/2*kx0:kx0:(Nx_hs/2-1)*kx0;
ky0    = 2*pi/By;
ky     = -Ny_hs/2*ky0:ky0:(Ny_hs/2-1)*ky0;

% Total depth of tunnel
H_cyl_tot = H;

% Define truncation criterion
delta_wh = 10^-6;
delta    = 10^-6;

abstand_wh = sqrt(H^2+(By/2-abs(y_Tc))^2);
abstand    = H-radius;

%% Check if tunnel or cyl. excavation (= trench)
if H >= radius
    geo.geo = 'tunnel'; 
else 
    geo.geo = 'trench'; 
end

%% Coordinates of points on ground surface
% y-, r-, phi-Coordinates on ground surface
N_OF   = Ny_hs;
y_OF   = ifftshift(-Ny_hs/2*dy:dy:Ny_hs/2*dy-dy);
r_OF   = sqrt(H^2 + (-y_Tc + y_OF).^2);
phi_OF = atan((-y_Tc + y_OF)/H)+(H>=0)*pi;

% All Points on ground surface
y_OF_ges   = y_OF;
r_OF_ges   = r_OF;
phi_OF_ges = phi_OF;


% Changed coordinates for trench
if strcmp(geo.geo,'trench')
    % y-,r-,phi-Coordinates fot Trench
    y_OF_ges = y_OF; 
    % Select points on and outside the cylinder
    % Points for radius and - radius 
    % radius... By/2-dy, -By/2... -radius
    y_OF   = y_OF_ges(  y_OF_ges<=+y_Tc-radius | y_OF_ges >=+y_Tc+radius);
    phi_OF = phi_OF_ges(y_OF_ges<=+y_Tc-radius | y_OF_ges >=+y_Tc+radius);
    r_OF   = r_OF_ges(  y_OF_ges<=+y_Tc-radius | y_OF_ges >=+y_Tc+radius);
    N_OF   = size(y_OF,2);
end

% y-,r-,phi-Coordinates on ground surface of adjacent repetion section
% Remark: Sign in y_LR_li & y_LR_re switched compared to COS definition
y_OF_li   = y_OF-By;                    % coordinates in left repetition section
r_OF_li   = sqrt(H^2+(-y_Tc+y_OF_li).^2);
phi_OF_li = atan((-y_Tc+y_OF_li)/H)+(H>=0)*pi;
y_OF_re   = y_OF+By;                    % coordinates in right repetition section
r_OF_re   = sqrt(H^2+(-y_Tc+y_OF_re).^2);
phi_OF_re = atan((-y_Tc+y_OF_re)/H)+(H>=0)*pi;

%% Coordinates of points on tunnel
% phi-,y-,z-Coordinates of cylindrical surface
N_LR   = Nphi_t;
phi_LR = fftshift(-Nphi_t/2*dphi:dphi:(Nphi_t/2-1)*dphi);

% Points of 0...-radius...0...radius....0
y_LR = -sin(phi_LR)*radius + y_Tc;
z_LR = cos(phi_LR)*radius+H;

% Phi defined from z anti clockwise
% y  <-----|.
%          | .
%          |  .
%          v   . phi
%          z

% All Points on cyl. surface
% Remark: Phi from 0... 2*pi
phi_LR_ges = 0:dphi:(Nphi_t-1)*dphi;

%  Changed coordinates for trench
if strcmp(geo.geo,'trench')
    % phi-,y-,z-Coordinates des Lochrands
    % phi_LR_ges=0:dphi:(Nphi_t-1)*dphi;
    phi_LR = 0:dphi:Nphi_t*dphi;        % for symmetric spline additional point at phi=+2pi
    y_LR   = -sin(phi_LR)*radius;       % Also here additional point -> 2 points at y=0
    z_LR   = cos(phi_LR)*radius+H;
    phi_LR = phi_LR(z_LR+10*eps>=0);    % Selection of phi_LR for z>=0
    y_LR   = y_LR(z_LR+10*eps>=0)+y_Tc; % Corresponding y_LR
    z_LR   = z_LR(z_LR+10*eps>=0);      % Corresponding z_LR
    N_LR   = size(phi_LR,2);            % Number of points on cyl. Trench
end

% phi-,y-,z-Coordinates of cylindrical surfaces of adjacent repetion sections
% Remark: Sign in y_LR_li & y_LR_re switched compared to COS definition
y_LR_li   = y_LR-By + y_Tc; % coordinates in right repetition section
z_LR_li   = z_LR;
r_LR_li   = sqrt((H-z_LR_li).^2+y_LR_li.^2);
phi_LR_li = atan(y_LR_li./(H-z_LR_li))+((H-z_LR_li+eps)>=0)*pi;
y_LR_re   = y_LR+By + y_Tc; % coordinates in left repetition section
z_LR_re   = z_LR;
r_LR_re   = sqrt((H-z_LR_re).^2+y_LR_re.^2);
phi_LR_re = atan(y_LR_re./(H-z_LR_re))+((H-z_LR_re+eps)>=0)*(-pi);


%% Discretization soil fine - coarse
dx_fine = dx;
dy_fine = dy;

Nx_fine = Nx_hs;
Ny_fine = Ny_hs;

x_fine  = x;
y_fine  = y;

[X_fine,Y_fine] = ndgrid(x_fine,y_fine);


%% Output
% Discretisation soil
dis_itm.dx      = dx;
dis_itm.dy      = dy;
dis_itm.dphi    = dphi;
dis_itm.Nphi_t  = Nphi_t;
dis_itm.kx      = kx;
dis_itm.kx0     = kx0;
dis_itm.ky      = ky;
dis_itm.ky0     = ky0;

dis_itm.f         = f;
dis_itm.omega     = omega;
dis_itm.Nf        = Nf;
dis_itm.f_tot     = f_tot;
dis_itm.omega_tot = omega_tot; 
dis_itm.Nf_tot    = Nf_tot;

% Discretisation fine/coarse and foundation
dis_itm.Nx_fine   = Nx_fine;
dis_itm.Ny_fine   = Ny_fine;
dis_itm.dx_fine   = dx_fine;
dis_itm.dy_fine   = dy_fine;
dis_fem.dy_FE     = dy_FE;

% Geometry
geo.x         = x;
geo.y         = y;
geo.H_cyl_tot = H_cyl_tot;

geo.y_OF      = y_OF;
geo.r_OF      = r_OF;
geo.phi_OF    = phi_OF;
geo.y_OF_ges  = y_OF_ges;
geo.r_OF_li   = r_OF_li;
geo.phi_OF_li = phi_OF_li;
geo.r_OF_re   = r_OF_re;
geo.phi_OF_re = phi_OF_re;
geo.N_OF      = N_OF;
geo.N_LR      = N_LR;

geo.phi_LR     = phi_LR;
geo.y_LR       = y_LR;
geo.z_LR       = z_LR;
geo.phi_LR_ges = phi_LR_ges;
geo.r_LR_li    = r_LR_li;
geo.phi_LR_li  = phi_LR_li;
geo.r_LR_re    = r_LR_re;
geo.phi_LR_re  = phi_LR_re;
                
% Discretisation fine/coarse and foundation                
geo.x_fine        = x_fine;
geo.y_fine        = y_fine;
geo.X_fine        = X_fine;
geo.Y_fine        = Y_fine;


% Truncation criteria
trunc.delta_wh   = delta_wh;
trunc.abstand_wh = abstand_wh;
trunc.delta      = delta;
trunc.abstand    = abstand;

end

