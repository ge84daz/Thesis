%__________________________________________________________________________
%
%       f_displ_sampling.m    
%
%       INPUT: 
%       kx,ky           = Inkrements in transformed domain
%       dk_coarse       = each dk_coarse-th value is evaluated for Kdyn within ]-kb kb[
%       dk_fine         = dk_fine amount of points within [-kb kb]
%       Bk              = increased wavenumber boundary kb = Bk*kr to get peak at kr   
%
%       OUTPUT: 
%       kx_ev, ky_ev         = points in x-, y-direction at which K_dyn is evaluated
%       kx_ev_tot, ky_ev_tot = original grid with increased sampling within [-kb kb]
%       kb                   = increased wavenumber boundary kb = Bk*kr
%
%       DESCRIPTION: 
%       - finer sampling within [-kb kb]    - take dk_fine*original points
%       - coarser sampling outside [-kb kb] - only each dk_coarse-th point
%       - interpolate to original grid for multiplication within f_displ_interp
%       - coarse: non-equidistant sampling
%       - no index: original sampling
%       - tot: original sampling with increased sampling within [-kb kb]
%
%       REMARK: 
%       Original Author: Tom Hicks 
%       Modified:   Tom Hicks                 
%       Date:       30-08-2023 
%       Changed:    26-09-2023
%__________________________________________________________________________

function [ky_ev, kx_ev, kb] = f_displ_sampling(mat, dis_itm, calc, calc_cont)
% function [ky_ev, ky_ev_tot, kx_ev, kx_ev_tot, kb] = f_displ_sampling(mat, dis_itm, calc, calc_cont)

dk_coarse = dis_itm.dk_coarse;          % only each dk_coarse-th value is evaluated for Kdyn outside [-kb kb]
dk_fine   = dis_itm.dk_fine;            % dk_fine amount of points within [-kb kb]
Bk = dis_itm.Bk;                        % determine boundary wavenumber kb = Bk*kr

freq = dis_itm.freq;
kx   = dis_itm.kx;
ky   = dis_itm.ky;
ky0  = dis_itm.ky0;
kx0  = dis_itm.kx0;
omega_calc = 2*pi*freq;                 % omega for which will be solved

% take lowest rayleigh wave velocity to set boundary of kb
switch calc.system    
    case 'tunnel'
        cr = [mat.hs.cr mat.fem_soil_T1.cr]; 
        kr = omega_calc(1,end)/min(abs(cr));
    case 'trench'
        cr = [mat.hs.cr mat.fem_soil_T1.cr]; 
        kr = omega_calc(1,end)/min(abs(cr));
end
kb = Bk * kr;                           % boundary to which Kdyn is interpolated with dk_coarse sample points

switch calc_cont.kx
    case 'kx=0'  % for 2.5D
                
        % Define ranges for different segments - only y-direction relevant
        kly = ky < -kb;                  % kl for k-low  - all values <-kb
        kmy = ky >= -kb & ky <= kb;      % km for k-mid  - all values within [-kb kb]
        khy = ky > kb;                   % kh for k-high - all values > kb
        
        % coarser sampling for kly and khy
        kly_ind_tot = find(kly);                         % indices of kly
        khy_ind_tot = find(khy);                         % indices of khy
        kmy_ind     = kmy;                               % indices of kmy

        kly_ind     = kly_ind_tot(1:dk_coarse:end);     % select every dk_coarse-th index
        khy_ind     = khy_ind_tot(1:dk_coarse:end);     % select every dk_coarse-th index

        kly_ev      = ky(kly_ind);                      % select every dk_coarse-th point of ky
        khy_ev      = ky(khy_ind);                      % select every dk_coarse-th point of ky

        %kly_ev_tot  = ky(kly_ind_tot);                  % points of ky <-kb
        %khy_ev_tot  = ky(khy_ind_tot);                  % points of ky > kb

        % finer sampling within [-kb kb]
        kmy_ev      = ky(kmy_ind);                       % points of ky within [-kb kb]
        % kmy_ev_tot  = kmy_ev(1):ky0/dk_fine:kmy_ev(end); % finer sampling of kmy
        % kmy_ev_tot(abs(kmy_ev_tot) < 1e-12) = 0;         % set kmy_ev_tot at kx=0 equal to zero
        
        % Create the total non-equidistant sampled vectors
        % ky_ev     = [kly_ev kmy_ev_tot khy_ev];         % grid over which K_dyn will be evaluated
        ky_ev     = [kly_ev kmy_ev khy_ev];         % grid over which K_dyn will be evaluated
        % ky_ev_tot = [kly_ev_tot kmy_ev_tot khy_ev_tot]; % original grid with finer sampling within [-kb kb]
        % kx_ev_tot = 0;                                   % for 2.5D case only kx = 0 is evaluated
        kx_ev     = zeros(1, length(kx));                
   
    case 'kx=all'  % for 3D
        
        % y-direction
        % Define ranges for different segments 
        kly = ky < -kb;                  % kl for k-low  - all values <-kb
        kmy = ky >= -kb & ky <= kb;      % km for k-mid  - all values within [-kb kb]
        khy = ky > kb;                   % kh for k-high - all values > kb

        % coarser sampling for kly and khy
        kly_ind_tot = find(kly);                         % indices of kly
        khy_ind_tot = find(khy);                         % indices of khy
        kmy_ind     = kmy;                               % indices of kmy
        
        kly_ind     = kly_ind_tot(1:dk_coarse:end);     % select every dk_coarse-th index
        khy_ind     = khy_ind_tot(1:dk_coarse:end);     % select every dk_coarse-th index

        kly_ev      = ky(kly_ind);                      % select every dk_coarse-th point of ky
        khy_ev      = ky(khy_ind);                      % select every dk_coarse-th point of ky

        % kly_ev_tot  = ky(kly_ind_tot);                  % points of ky <-kb
        % khy_ev_tot  = ky(khy_ind_tot);                  % points of ky > kb

        % finer sampling within [-kb kb]
        kmy_ev      = ky(kmy_ind);                       % points of ky within [-kb kb]
        % kmy_ev_tot = kmy_ev(1):ky0/dk_fine:kmy_ev(end); % finer sampling of kmy
        % kmy_ev_tot(abs(kmy_ev_tot) < 1e-12) = 0;         % set kmy_ev_tot at kx=0 equal to zero
        % Create the total non-equidistant sampled vectors
        ky_ev     = [kly_ev kmy_ev khy_ev];         % grid over which K_dyn will be evaluated
        % ky_ev     = [kly_ev kmy_ev_tot khy_ev];         % grid over which K_dyn will be evaluated
        % ky_ev_tot = [kly_ev_tot kmy_ev_tot khy_ev_tot]; % original grid with finer sampling within [-kb kb]

        % x-direction
        % Define ranges for different segments 
        klx = kx < -kb;                  % kl for k-low  - all values <-kb
        kmx = kx >= -kb & kx <= kb;      % km for k-mid  - all values within [-kb kb]
        khx = kx > kb;                   % kh for k-high - all values > kb

        % coarser sampling for klx and khx
        klx_ind_tot = find(klx);                         % indices of klx
        khx_ind_tot = find(khx);                         % indices of khx
        kmx_ind     = kmx;                               % indices of kmx
        
        klx_ind     = klx_ind_tot(1:dk_coarse:end);     % select every dk_coarse-th index
        khx_ind     = khx_ind_tot(1:dk_coarse:end);     % select every dk_coarse-th index

        klx_ev      = kx(klx_ind);                      % select every dk_coarse-th point of kx
        khx_ev      = kx(khx_ind);                      % select every dk_coarse-th point of kx

        % klx_ev_tot  = kx(klx_ind_tot);                  % points of kx <-kb
        % khx_ev_tot  = kx(khx_ind_tot);                  % points of kx > kb

        % finer sampling within [-kb kb]
        kmx_ev      = kx(kmx_ind);                      % points of kx within [-kb kb]
        % kmx_ev_tot = kmx_ev(1):kx0/dk_fine:kmx_ev(end); % finer sampling of kmx
        % kmx_ev_tot(abs(kmx_ev_tot) < 1e-12) = 0;        % set kmx_ev_tot at ky=0 equal to zero
        % Create the total non-equidistant sampled vectors
        kx_ev     = [klx_ev kmx_ev khx_ev];         % grid over which K_dxn will be evaluated
        % kx_ev     = [klx_ev kmx_ev_tot khx_ev];         % grid over which K_dxn will be evaluated
        % kx_ev_tot = [klx_ev_tot kmx_ev_tot khx_ev_tot]; % original grid with finer sampling within [-kb kb]

end
end



