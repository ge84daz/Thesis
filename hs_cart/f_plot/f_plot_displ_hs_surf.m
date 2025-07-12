%__________________________________________________________________________
%
%       f_plot_displ_hs_surf.m    
%
%       INPUT: 
%       displ       = Displ. Hs_1L in (kx,ky,z,omega) in z or y-direction
%       dis_itm     = ITM discretization
%       geo         = Geometry
%
%       OUTPUT: 
%
%       DESCRIPTION: 
%       Plot of displ at soil surface z=0
%       - for wavenumber and space domain
%       - for frequency spectrum and time course
%       - displ course for one f and tmax
%
%
%       REMARK: 
%       Original Author:  Julian Freisinger  (f_plot_displ_hs_surf)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________



function [] = f_plot_displ_hs_surf(displ_dir,load_dir,dis_itm,geo,displ,loading,plot_cont,path,calc)

%% 1.) Initialize: Input

f_tot     = dis_itm.f_tot;
% x_eval_xy = geo.x_eval_xy;
% y_eval_xy = geo.y_eval_xy;

kx       = dis_itm.kx;
ky       = dis_itm.ky;
kx_ev    = dis_itm.kx;
ky_ev    = dis_itm.ky;

x           = geo.x;
y           = geo.y;
xev         = geo.x;
yev         = geo.y;

calc_sys = strrep('hs','_','');

displ_dir_sub = strrep(displ_dir,'u','');

if strcmpi(loading.type,'harmonic')
   f_harm     = dis_itm.f_harm;
   t_harm     = dis_itm.t_harm;
   loc_f      = loading.harmonic.loc_f;
   loc_f_harm = loading.harmonic.loc_f_harm; 
   Pt_harm    = loading.Pt_harm;
   Pf_harm    = loading.Pf_harm;
   by_itm_hs1L= loading.by_itm_hs1L;
else
   % not implemented for transient load as not meaningful for many excitation frequencies 
end

varname1 = [displ_dir load_dir '_xyf_z0'];
varname2 = [displ_dir load_dir '_kxkyf_z0'];
varname3 = [displ_dir load_dir '_xyt_z0'];

u_xyf_z0   = displ.Pt.(varname1);
u_kxkyf_z0 = displ.Pt.(varname2);
u_xyt_z0   = displ.Pt.(varname3);

%% 2.) Plots Displ Soil
  
if  (size(xev,2)>=2 && size(yev,2)>=2)
    %% 2D Abs,Re,Im kx,ky and x,y for f and t
    if strcmpi(plot_cont.lowMemory,'no')
    %-----------------------------------------------------
    mainfigurename = [calc_sys ' - Displ. SOIL 2D u' displ_dir_sub '(z=0) due to ' load_dir];
    fig = figure('Name',mainfigurename);
    fig.Position = [632,260,1403,1040];
    %-----------------------------------------------------
    %
    subplot(3,3,1)
    dkxky   = 2;
    [KX,KY] = ndgrid(kx_ev(1:dkxky:end),ky_ev(1:dkxky:end));
    mesh(KX,KY,abs(u_kxkyf_z0(1:dkxky:end,1:dkxky:end,loc_f(1))),'FaceColor',plot_cont.FaceColor);
    xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$|u_' displ_dir_sub '(k_x,k_y)|$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,2)
    mesh(KX,KY,real(u_kxkyf_z0(1:dkxky:end,1:dkxky:end,loc_f(1))),'FaceColor',plot_cont.FaceColor);
    xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Re} \;\; u_' displ_dir_sub '(k_x,k_y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,3)
    mesh(KX,KY,imag(u_kxkyf_z0(1:dkxky:end,1:dkxky:end,loc_f(1))),'FaceColor',plot_cont.FaceColor);
    xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Im} \;\; u_' displ_dir_sub '(k_x,k_y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,4)
    dxy   = 2;
    [X,Y] = ndgrid(xev(1:dkxky:end),yev(1:dkxky:end));
    mesh(X,Y,abs(u_xyf_z0(1:dxy:end,1:dxy:end,loc_f_harm(1))),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$|u_' displ_dir_sub '(x,y)|$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,5)
    mesh(X,Y,real(u_xyf_z0(1:dxy:end,1:dxy:end,loc_f_harm(1))),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Re} \;\;u_' displ_dir_sub '(x,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,6)
    mesh(X,Y,imag(u_xyf_z0(1:dxy:end,1:dxy:end,loc_f_harm(1))),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Im} \;\;u_' displ_dir_sub '(x,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,7)
    % find tmax and its index
    [max_uxyt,ind_max_uxyt] = max(u_xyt_z0(:));
    [max_row,max_col,ind_t_max] = ind2sub(size(u_xyt_z0),ind_max_uxyt);
    
    mesh(X,Y,abs(u_xyt_z0(1:dxy:end,1:dxy:end,ind_t_max)),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$|u_' displ_dir_sub '(x,y)|$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,8)
    mesh(X,Y,real(u_xyt_z0(1:dxy:end,1:dxy:end,ind_t_max)),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Re} \;\;u_' displ_dir_sub '(x,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,9)
    mesh(X,Y,imag(u_xyt_z0(1:dxy:end,1:dxy:end,ind_t_max)),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Im} \;\;u_' displ_dir_sub '(x,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    zlim([min(real(u_xyt_z0(1:dxy:end,1:dxy:end,ind_t_max)),[],'all') max(real(u_xyt_z0(1:dxy:end,1:dxy:end,ind_t_max)),[],'all')])
    
    savename = ['\displ_soil_kxky_xy_f_t'];
    saveas(fig,[path.figures savename],'fig')
    
    if  strcmpi(plot_cont.export,'yes');
        exportfig(gcf,[path.figures savename],plot_cont.exportfig);end
    if ~strcmpi(plot_cont.displ_soil,'yes');close(fig);end
    end
    
    %% 1D Abs,Re,Im seperated over ky/y
    %-----------------------------------------------------
    mainfigurename = [calc_sys ' - Displ. soil 1D u' displ_dir_sub '(z=0) due to ' load_dir '- 01'];
    fig = figure('Name',mainfigurename);
    fig.Position = [632,260,1403,1040];
    %-----------------------------------------------------
    %
    subplot(3,3,1)
    plot(ky_ev,abs(u_kxkyf_z0(kx_ev==0,:,loc_f(1))),'Color',plot_cont.color_cyan);
    xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$|u_' displ_dir_sub '(k_x=0,k_y)|$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,2)
    plot(ky_ev,real(u_kxkyf_z0(kx_ev==0,:,loc_f(1))),'Color',plot_cont.color_cyan);
    xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$\mathrm{Re}\;\; u_' displ_dir_sub '(k_x=0,k_y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,3)
    plot(ky_ev,imag(u_kxkyf_z0(kx_ev==0,:,loc_f(1))),'Color',plot_cont.color_cyan);
    xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$\mathrm{Im}\;\; u_' displ_dir_sub '(k_x=0,k_y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,4)
    plot(yev,abs(u_xyf_z0(xev==0,:,loc_f_harm(1))),'Color',plot_cont.color_cyan);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$|u_' displ_dir_sub '(x=0,y)|$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,5)
    plot(yev,real(u_xyf_z0(xev==0,:,loc_f_harm(1))),'Color',plot_cont.color_cyan);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$\mathrm{Re}\;\; u_' displ_dir_sub '(x=0,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,6)
    plot(yev,imag(u_xyf_z0(xev==0,:,loc_f_harm(1))),'Color',plot_cont.color_cyan);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$\mathrm{Im}\;\; u_' displ_dir_sub '(x=0,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,7)
    % find tmax and its index
    [max_uxyt,ind_max_uxyt] = max(u_xyt_z0(:));
    [max_row,max_col,ind_t_max] = ind2sub(size(u_xyt_z0),ind_max_uxyt);
    
    plot(yev,abs(u_xyt_z0(xev==0,:,ind_t_max)),'Color',plot_cont.color_cyan);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$|u_' displ_dir_sub '(x=0,y)|$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,8)
    plot(yev,real(u_xyt_z0(xev==0,:,ind_t_max)),'Color',plot_cont.color_cyan);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$\mathrm{Re}\;\; u_' displ_dir_sub '(x=0,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,3,9)
    plot(yev,imag(u_xyt_z0(xev==0,:,ind_t_max)),'Color',plot_cont.color_cyan);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$\mathrm{Im}\;\; u_' displ_dir_sub '(x=0,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    ylim([min(real(u_xyt_z0(xev==0,:,ind_t_max)),[],'all') max(real(u_xyt_z0(xev==0,:,ind_t_max)),[],'all')])
    
    savename = ['\displ_soil_kx0x0_kyy'];
    saveas(fig,[path.figures savename],'fig')
    
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig);end
    if ~strcmpi(plot_cont.displ_soil,'yes');close(fig);end
    
    %% 1D Re,Im u over ky/y and f/t for kx=0/x=0
    %-----------------------------------------------------
    mainfigurename = [calc_sys ' - Displ. soil 1D u' displ_dir_sub '(z=0) due to ' load_dir ' - 02'];
    fig = figure('Name',mainfigurename);
    fig.Position = [632,260,1403,1040];
    %-----------------------------------------------------
    %
    subplot(3,2,1)
    hold on
    plot(ky_ev,real(u_kxkyf_z0(kx_ev==0,:,loc_f(1))),'Color',plot_cont.color_cyan);
    plot(ky_ev,imag(u_kxkyf_z0(kx_ev==0,:,loc_f(1))),'Color',plot_cont.color_orange2);
    xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$u_' displ_dir_sub '(k_x=0,k_y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    %
    subplot(3,2,2)
    stem(f_tot,real(squeeze(u_kxkyf_z0(kx_ev==0,ky_ev==0,:))),'Color',plot_cont.color_cyan,'Marker','.');
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(k_x=0,k_y=0,z=0,f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(k_x=0,k_y=0,z=0,f)$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(3,2,3)
    hold on
    plot(yev,real(u_xyf_z0(xev==0,:,loc_f_harm(1))),'Color',plot_cont.color_cyan);
    plot(yev,imag(u_xyf_z0(xev==0,:,loc_f_harm(1))),'Color',plot_cont.color_orange2);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(x=0,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    %
    subplot(3,2,4)
    hold on
    stem(f_harm,squeeze(real(u_xyf_z0(xev==0,yev==0,:))),'Color',plot_cont.color_cyan,'Marker','.');
    stem(f_harm,squeeze(imag(u_xyf_z0(xev==0,yev==0,:))),'Color',plot_cont.color_orange2,'Marker','.');
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y=0,z=0,f)$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    %
    subplot(3,2,5)
    % find tmax and its index
    [max_uxyt,ind_max_uxyt] = max(u_xyt_z0(:));
    [max_row,max_col,ind_t_max] = ind2sub(size(u_xyt_z0),ind_max_uxyt);
    hold on
    plot(yev,real(u_xyt_z0(xev==0,:,ind_t_max)),'Color',plot_cont.color_cyan);
    plot(yev,imag(u_xyt_z0(xev==0,:,ind_t_max)),'Color',plot_cont.color_orange2);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(x=0,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    %
    subplot(3,2,6)
    % find tmax and its index
    [max_uxyt,ind_max_uxyt] = max(u_xyt_z0(:));
    [max_row,max_col,ind_t_max] = ind2sub(size(u_xyt_z0),ind_max_uxyt);
    hold on
    plot(t_harm,squeeze(real(u_xyt_z0(xev==0,yev==0,:))),'Color',plot_cont.color_cyan);
    plot(t_harm,squeeze(imag(u_xyt_z0(xev==0,yev==0,:))),'Color',plot_cont.color_orange2);
    xlabel('$t$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y=0,z=0,t)$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    
    savename = ['\displ_soil_kx0x0_kyy_f_t'];
    saveas(fig,[path.figures savename],'fig')
    
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig);end
    if ~strcmpi(plot_cont.displ_soil,'yes');close(fig);end
    
    %% Overview displacements
    %-----------------------------------------------------
    mainfigurename = [calc_sys ' - Displ. soil u' displ_dir_sub '(z=0) due to ' load_dir ' - Overview'];
    fig = figure('Name',mainfigurename);
    fig.Position = [632,260,1403,1040];
    %-----------------------------------------------------
    subplot(4,2,1)
    plot(t_harm,Pt_harm)
    xlabel('$t$',    'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P(t)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    title('Time History','Interpreter',plot_cont.interpreter)
    %
    subplot(4,2,2)
    hold on
    stem(f_harm,real(Pf_harm),'Marker','.')
    stem(f_harm,imag(Pf_harm),'Marker','.')
    xlabel('$f$',    'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    ylabel('$P(f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    title('Frequency course','Interpreter',plot_cont.interpreter)
    %
    subplot(4,2,3)
    % find tmax and its index
    [max_uxyt,ind_max_uxyt] = max(u_xyt_z0(:));
    [max_row,max_col,ind_t_max] = ind2sub(size(u_xyt_z0),ind_max_uxyt);
    hold on
    plot(t_harm,squeeze(real(u_xyt_z0(xev==0,yev==0,:))),'Color',plot_cont.color_cyan);
    plot(t_harm,squeeze(imag(u_xyt_z0(xev==0,yev==0,:))),'Color',plot_cont.color_orange2);
    xlabel('$t$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y=0,z=0,t)$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    %
    subplot(4,2,4)
    hold on
    stem(f_harm,squeeze(real(u_xyf_z0(xev==0,yev==0,:))),'Color',plot_cont.color_cyan,'Marker','.');
    stem(f_harm,squeeze(imag(u_xyf_z0(xev==0,yev==0,:))),'Color',plot_cont.color_orange2,'Marker','.');
    xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y=0,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    %
    subplot(4,2,5)
    % find tmax and its index
    [max_uxyt,ind_max_uxyt] = max(u_xyt_z0(:));
    [max_row,max_col,ind_t_max] = ind2sub(size(u_xyt_z0),ind_max_uxyt);
    hold on
    plot(yev,real(u_xyt_z0(xev==0,:,ind_t_max)),'Color',plot_cont.color_cyan);
    plot(yev,imag(u_xyt_z0(xev==0,:,ind_t_max)),'Color',plot_cont.color_orange2);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(x=0,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    %
    subplot(4,2,6)
    hold on
    plot(yev,real(u_xyf_z0(xev==0,:,loc_f_harm(1))),'Color',plot_cont.color_cyan);
    plot(yev,imag(u_xyf_z0(xev==0,:,loc_f_harm(1))),'Color',plot_cont.color_orange2);
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$ u_' displ_dir_sub '(x=0,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x=0,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    legend('Re','Im','Interpreter',plot_cont.interpreter)
    %
    subplot(4,2,7)
    %
    dxy   = 2;
    [X,Y] = ndgrid(xev(1:dxy:end),yev(1:dxy:end));
    % find tmax and its index
    [max_uxyt,ind_max_uxyt] = max(u_xyt_z0(:));
    [max_row,max_col,ind_t_max] = ind2sub(size(u_xyt_z0),ind_max_uxyt);
    %
    mesh(X,Y,real(u_xyt_z0(1:dxy:end,1:dxy:end,ind_t_max)),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Re} \;\;u_' displ_dir_sub '(x,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x,y,z=0,t_{max})$'],'Interpreter',plot_cont.interpreter)
    %
    subplot(4,2,8)
    mesh(X,Y,real(u_xyf_z0(1:dxy:end,1:dxy:end,loc_f_harm(1))),'FaceColor',plot_cont.FaceColor);
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Re} \;\;u_' displ_dir_sub '(x,y)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    title([calc_sys ' - $u_{' displ_dir_sub ',' load_dir '}(x,y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    
    savename = ['\displ_soil_x0y_xy_P_f_t_overview'];
    saveas(fig,[path.figures savename],'fig')
    
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig);end
    if ~strcmpi(plot_cont.displ_soil,'yes');close(fig);end
    
end

end