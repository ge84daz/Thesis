%__________________________________________________________________________
%
%        f_load_hscyl_t_xy_harm.m    
%
%       INPUT: 
%       - loading, dis_itm, calc_cont. plot_cont, path
%
%
%       DESCRIPTION: 
%       - plot time history of load for harmonic case
%       - same function as for Hs1L
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_load_hscyl_t_xy_harm)
%       Modified:   Tom Hicks                 
%       Date:       25-09-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________

function [loading,dis_itm] = f_load_hscyl_t_xy_harm(loading,dis_itm,~,plot_cont,path)


% Initialize
f_calc = loading.harmonic.f;
P0     = loading.harmonic.P0;
nT     = loading.harmonic.nT;
Nt     = loading.harmonic.Nt;

f_tot  = dis_itm.f_tot;
Nf_tot = dis_itm.Nf_tot;

% Frequency spectrum of load

% location of load f_calc in f_tot
f_calc_tot     = [-flip(f_calc) f_calc];
[~,~,loc_ftot] = intersect(f_calc_tot,f_tot);

% Set vector of amplitude at load frequencies within f_tot
Pf_tot = zeros(1,1,Nf_tot);
P0_tot = [flip(P0) P0];

Pf_tot(loc_ftot) = P0_tot;

%% Time/frequency discretization for IFFT -> u(t)

fmin   = min(f_calc);
T_harm = nT*1/fmin;

dt_harm = T_harm/Nt;        % sample rate
df_harm = 1/Nt/dt_harm;     % Nyquist

% finer f discretization for finer t in u(t) after IFFT
f_harm  = -Nt/2*df_harm:df_harm:(Nt/2-1)*df_harm;
t_harm  = -Nt/2*dt_harm:dt_harm:(Nt/2-1)*dt_harm;   % t = [-T/2;T/2]

Nt_harm = length(t_harm);

% location of load frequencies f_calc in f_harm
[~,~,loc_f_harm] = intersect(f_calc_tot,f_harm);

%% Time dependent load function
Pt_harm = zeros(1,Nt);
for nnf = 1:length(f_calc)
    Pt_harm = Pt_harm + P0(nnf)*cos(2*pi*f_calc(nnf)*t_harm);
end
% Adaption of load amplitude as:
% Pt_harm = IFFT Pf_harm with amplitude 1 in freq domain
% -> Pt_harm actually has amplitude 2
Pt_harm = 2*Pt_harm;

% Excitation spectrum for plot
Pf_harm = 1/Nt_harm*fftshift(fft(ifftshift(Pt_harm)));%dt_harm*

%% Output
loading.Pf_tot  = Pf_tot;       % load vector for u(f)=TF(f)*Pf_tot with f=f_tot
loading.Pt_harm = Pt_harm;      % with t = t_harm
loading.Pf_harm = Pf_harm;      % with f = f_harm

loading.loc_f_tot  = loc_ftot;     % location of f_calc in f_tot
loading.loc_f_harm = loc_f_harm;   % loaction of f_calc in f_harm

dis_itm.t_harm   = t_harm;
dis_itm.dt_harm  = dt_harm;
dis_itm.df_harm  = df_harm;
dis_itm.T_harm   = T_harm;
dis_itm.f_harm   = f_harm;
dis_itm.Nt_harm  = Nt_harm;

%% Plots

%-----------------------------------------------------
mainfigurename = 'Hs_cyl - ITM P(t) and P(f)';
fig = figure('Name',mainfigurename);
fig.Position = plot_cont.fig_position; %[632,260,1403,1040];
%-----------------------------------------------------
sgtitle('Harmonic $P(t)$, $P(f$)','FontSize',10,'interpreter','latex')
%
subplot(1,2,1)
plot(t_harm,Pt_harm)
xlabel('$t$',    'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel('$P(t)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
%
subplot(1,2,2)
stem(f_harm,real(Pf_harm),'Marker','.')
xlabel('$f$',    'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel('$P(f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)

savename = '\load_Pz_t_f';
% saveas(fig,[path.figures savename],'fig')
%
if strcmpi(plot_cont.export,'yes')
    exportfig(gcf,[path.figures savename],'Width',plot_cont.exp_width,'Height',plot_cont.exp_heigth,'Fontmode','fixed','Fontsize',plot_cont.Fontsize,'Color','cmyk','resolution',plot_cont.resolution);
end
 switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end


end

