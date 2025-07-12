%__________________________________________________________________________
%
%         f_plot_displ_hscyl_t_harm.m    
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
%       - plot of time history of displacements
%
%       DESCRIPTION: 
%
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_plot_displ_hscyl_t_harm)
%       Modified:   Tom Hicks                 
%       Date:       25-09-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________



function [] = f_plot_displ_hscyl_t_harm(displ_dir,load_dir,calc_state,displ,dis_itm,geo,loading,calc,plot_cont,path)

% Initialize
f_tot      = dis_itm.f_tot;
t_harm     = dis_itm.t_harm;
f_harm     = dis_itm.f_harm;
Pt_harm    = loading.Pt_harm;
Pf_harm    = loading.Pf_harm;

x_eval_xy  = geo.x_eval_xy;
y_eval_xy  = geo.y_eval_xy;

xev       = geo.x(abs(geo.x)<=x_eval_xy);
yev       = geo.y(abs(geo.y)<=y_eval_xy);


calc_sys      = ['Hscyl - ' calc.system ];
displ_dir_sub = strrep(displ_dir,'u','');

u_xyf_itm_z0_TF      = getfield(displ.hscyl,calc.system,'itm_Pt',[displ_dir load_dir '_xyf_hscyl_itm_z0_TF_nint']);
u_xyf_itm_z0_Pf_harm = getfield(displ.hscyl,calc.system,'itm_Pt',[displ_dir load_dir '_xyf_hscyl_itm_z0']);
u_xyt_itm_z0_Pf_harm = getfield(displ.hscyl,calc.system,'itm_Pt',[displ_dir load_dir '_xyt_hscyl_itm_z0']);



%% Plot:  TF(f) u(f)
if ~strcmpi(calc.structure,'open_trench')
%-----------------------------------------------------
fig          = figure;
fig.Position = plot_cont.fig_position;
fig.Name     = [calc_sys ' - Resp. soil P(f) TF(f) u(f)'];
%-----------------------------------------------------
sgtitle(['$u_{' displ_dir_sub '}(x=y=z=0)$ of SOIL due to ' loading.distr_xy ' ' loading.type ' load $P_' load_dir(2) '(z=0)$'],'interpreter','latex','FontSize',12)
%
subplot(2,3,1)
plot(f_harm,real(Pf_harm))
xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel(['$ \mathrm{Re}\;\; P_' load_dir(2) ' (f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
title(['Load~spectrum $P_' displ_dir_sub '(f)$'],'Interpreter','latex','FontSize',10)
%
subplot(2,3,4)
plot(f_harm,imag(Pf_harm))
xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel(['$ \mathrm{Im}\;\; P_' load_dir(2) ' (f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
%
%
subplot(2,3,2)
hold on
plot(f_tot,real(squeeze(u_xyf_itm_z0_TF(xev==0,yev==0,:))),'Marker','.','LineStyle','--') % initial -> for interp 'no' and harmonic
xlabel('$f$','Interpreter','latex')
ylabel(['$\mathrm{Re}\;\; \mathrm{TF} \; u_' displ_dir_sub '(f)$'],'Interpreter','latex')
legend('initial','Interpreter','latex')
title(['Transfer Function $u_' displ_dir_sub '(f)$'],'Interpreter','latex')
%
subplot(2,3,5)
hold on
plot(f_tot,imag(squeeze(u_xyf_itm_z0_TF(xev==0,yev==0,:))),'Marker','.','LineStyle','--') % initial
xlabel('$f$','Interpreter','latex')
ylabel(['$\mathrm{Im} \;\; \mathrm{TF} \;  u_' displ_dir_sub '(f)$'],'Interpreter','latex')
legend('initial','Interpreter','latex')
%
subplot(2,3,3)
hold on
plot(f_harm,real(squeeze((u_xyf_itm_z0_Pf_harm(xev==0,yev==0,:))))) % non interpolated
xlabel('$f$','Interpreter','latex')
ylabel(['$\mathrm{Re}\; u_' displ_dir_sub '(f)$'],'Interpreter','latex')
title(['Displ.~spectrum $u_' displ_dir_sub '(f) = \mathrm{TF}(f)\cdot P_' load_dir(2) '(f)$'],'Interpreter','latex')
%
subplot(2,3,6)
hold on
plot(f_harm,imag(squeeze((u_xyf_itm_z0_Pf_harm(xev==0,yev==0,:))))) % non interpolated
xlabel('$f$','Interpreter','latex')
ylabel(['$\mathrm{Im}\; u_' displ_dir_sub '(f)$'],'Interpreter','latex')

if     strcmpi(calc_state,'ini'); savename = ['\displ_hscyl_harm_' displ_dir load_dir '_P_TF_u_x0y0z0f_orig'];
elseif strcmpi(calc_state,'post');savename = ['\displ_hscyl_harm_post_' displ_dir load_dir '_P_TF_u_x0y0z0f_orig'];   
end

saveas(fig,[path.figures savename],'fig')

if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
if ~strcmpi(displ_dir,'uz'); switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end;end
end
%% Plot: Response P(f),u(f),u(t)
if ~strcmpi(calc.structure,'open_trench')
%-----------------------------------------------------
mainfigurename = [calc_sys ' - u' displ_dir_sub ' SOIL due to ' load_dir ' - u(x=0,y=0,z=0,f),u(x=0,y=0,z=0,t)'];
fig = figure('Name',mainfigurename);
fig.Position = plot_cont.fig_position;
%-----------------------------------------------------
sgtitle(['$P_' load_dir(2) '$, TF and $u_{' displ_dir_sub '}$ of SOIL at $x=y=z=0$ for ' calc_sys],'interpreter','latex','FontSize',12)
%
subplot(3,2,1)
plot(t_harm,Pt_harm)
xlabel('$t$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel(['$P_' load_dir(2) '(t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
%
subplot(3,2,2)
stem(f_harm,real(Pf_harm),'Marker','.')
xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel(['$P_' load_dir(2) '(f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
%
subplot(3,2,4)
hold on
plot(f_tot,real(squeeze(u_xyf_itm_z0_TF(xev==0,yev==0,:))),'Color',plot_cont.color_cyan)
plot(f_tot,imag(squeeze(u_xyf_itm_z0_TF(xev==0,yev==0,:))),'Color',plot_cont.color_orange2)
xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel(['$\mathrm{TF} \;\; u_' displ_dir_sub '(x,y,z,f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
legend('Re','Im')
%
subplot(3,2,6)
hold on
stem(f_harm,real(squeeze(u_xyf_itm_z0_Pf_harm(xev==0,yev==0,:))),'Marker','.')
stem(f_harm,imag(squeeze(u_xyf_itm_z0_Pf_harm(xev==0,yev==0,:))),'Marker','.')
xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel(['$u_' displ_dir_sub '(f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
legend('Re','Im','Interpreter',plot_cont.interpreter)
title(['$u_' displ_dir_sub '(x=0,y=0,z=0,f)$'],'Interpreter',plot_cont.interpreter)
%
subplot(3,2,5)
hold on
plot(t_harm,real(squeeze(u_xyt_itm_z0_Pf_harm(xev==0,yev==0,:))))
plot(t_harm,imag(squeeze(u_xyt_itm_z0_Pf_harm(xev==0,yev==0,:))))
xlabel('$t$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
ylabel(['$u_' displ_dir_sub '(t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
legend('Re','Im','Interpreter',plot_cont.interpreter)
title(['$u_' displ_dir_sub '(x=0,y=0,z=0,t)$'],'Interpreter',plot_cont.interpreter)
%
if     strcmpi(calc_state,'ini'); savename = ['\displ_hscyl_harm_' displ_dir load_dir '_P_TF_u_x0y0z0_f_t'];
elseif strcmpi(calc_state,'post');savename = ['\displ_hscyl_harm_post_' displ_dir load_dir '_P_TF_u_x0y0z0_f_t'];   
end

saveas(fig,[path.figures savename],'fig')
%
if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
if ~strcmpi(displ_dir,'uz'); switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end;end
end
%% Plot: Anmiated plot u(t)

if strcmpi(plot_cont.u_anim,'yes') && (size(xev,2)>=2 && size(yev,2)>=2) && ~strcmpi(calc.structure,'open_trench')
    %-----------------------------------------------------
    fig          = figure;
    fig.Position = plot_cont.fig_position;
    fig.Name     = [calc_sys ' - u' displ_dir_sub ' SOIL due to ' load_dir];
    %-----------------------------------------------------
    
    %             pause_t = 1/Nt_harm*dt_harm;
    
    sgtitle(['$u_{' displ_dir_sub '}(x=y=z=0)$ SOIL due to ' loading.distr_xy ' ' loading.type ' load $P_' load_dir(2) '(z=0)$'],'interpreter','latex','FontSize',12)
    
    [X,Y] = ndgrid(xev(1:2:end),yev(1:2:end));
    Z     = real(u_xyt_itm_z0_Pf_harm(1:2:end,1:2:end,:));
    Z0    = Z(:,:,1);
    minz  = min(Z,[],'all');
    maxz  = max(Z,[],'all');
    %
    subplot(2,2,1)
    subpl_1 = mesh(X,Y,Z0,'FaceColor',plot_cont.FaceColor);
    view([-22.500000000000004,10.199999999999999])
    zlim([minz maxz])
    colormap(jet)
    clim([minz/4 maxz/4])
    xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    zlabel(['$\mathrm{Re} \;\; u_' displ_dir_sub '(x,y,t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ax = gca;
    ax.ZAxis.TickLabelFormat = '%.1f';
    
    % 2D u(x=0,y,z=0,t)
    subplot(2,2,2);
    u_x0yt = squeeze(real(u_xyt_itm_z0_Pf_harm(xev==0,:,:)));
    %     plot(t_harm,u_x0y0t);
    miny3  = min(u_x0yt,[],'all');
    maxy3  = max(u_x0yt,[],'all');
    subpl_2 = plot(yev,u_x0yt(:,1));
    ylim([miny3 maxy3])
    xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel(['$\mathrm{Re} \;\; u_' displ_dir_sub '(x=0,y,t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    
    % 2D: u(t)
    subplot(2,2,3);
    hold on
    %             h1 = animatedline;
    Pplot2 = real(squeeze(u_xyt_itm_z0_Pf_harm(xev==0,yev==0,:)));
    miny1 = min(Pplot2);
    maxy1 = max(Pplot2);
    plot(t_harm,Pplot2)
    s1 = scatter(t_harm(1),Pplot2(1),'MarkerEdgeColor',plot_cont.color_cyan,'MarkerFaceColor',plot_cont.color_cyan);
    axis([t_harm(1) t_harm(end) miny1 maxy1])
    xlabel('$t$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$P(t)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    % 2D: P(t)
    subplot(2,2,4);
    hold on
    %             h2 = animatedline;
    Pt = real(Pt_harm);
    plot(t_harm,Pt)
    s2 = scatter(t_harm(1),Pt(1),'MarkerEdgeColor',plot_cont.color_cyan,'MarkerFaceColor',plot_cont.color_cyan);
    miny2 = min(Pt_harm);
    maxy2 = max(Pt_harm);
    axis([t_harm(1) t_harm(end) miny2 maxy2])
    xlabel('$t$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    ylabel('$u(t)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    
    for nnt = 1:1:length(t_harm)
        
        % ISO displ
        subpl_1.ZData = real(u_xyt_itm_z0_Pf_harm(1:2:end,1:2:end,nnt)); % current Data
        
        % 2D u(x=0,y,z=0,t)
        subpl_2.YData = real(squeeze(u_xyt_itm_z0_Pf_harm(xev==0,:,nnt)));
        
        % 2D u(t)
        subplot(2,2,3)
        s1.XData = t_harm(nnt);
        s1.YData = real(squeeze(u_xyt_itm_z0_Pf_harm(xev==0,yev==0,nnt)));
        
        % 2D p(t)
        subplot(2,2,4)
        s2.XData = t_harm(nnt);
        s2.YData = real(Pt(nnt));
        
        del = 0.01;
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if     strcmpi(calc_state,'ini');  savename = ['\displ_hscyl_harm_' displ_dir load_dir '_u_xyz0_gif.gif'];
        elseif strcmpi(calc_state,'post'); savename = ['\displ_hscyl_harm_post_' displ_dir load_dir '_u_xyz0_gif.gif'];
        end
       
        filename = [path.figures savename];
        
        if nnt == 1
            imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
        end
    end
    
end

end