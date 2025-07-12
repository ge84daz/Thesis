%__________________________________________________________________________
%
%       f_ifft_t_displ_Pt.m    
%
%       INPUT: 
%       path      = pathes
%       plot_cont = plot control parameters
%       loading   = load amplitude, width, position
%       dis_itm   = discretization
%       geo       = geometry
%       mat       = material
%       displ     = displacements from ITM calc
%
%       OUTPUT: 
%
%       DESCRIPTION: 
%
%       REMARK: 
%       Original Author: Julian Freisinger  (f_ifft_t_displ_Pt)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________


function [displ] = f_ifft_t_displ_Pt(displ_dir,load_dir,displ,dis_itm,geo,loading,mat,calc,plot_cont,path)

%%  Initialize: Input

f_tot    = dis_itm.f_tot;

% x_eval_xy   = geo.x_eval_xy;
% y_eval_xy   = geo.y_eval_xy;
x           = geo.x;
y           = geo.y;
xev         = x;
yev         = y;

Pf = loading.Pf;

calc_sys = strrep('hs','_',' ');



varname1 = [displ_dir load_dir '_xyf_z0'];
varname2 = [displ_dir load_dir '_kxkyf_z0'];

u_xy_z0_omega_TF_nint   = displ.TF.(varname1)(:,:,:); % spatial domain over xev, yev
u_kxky_z0_omega_TF_nint = displ.TF.(varname2)(:,:,:); % wavenumber domain over xev, yev

% u_xy_z0_omega_TF_nint   = displ.TF.(varname1)(abs(x)<=x_eval_xy,abs(y)<=y_eval_xy,:); %spatial domain over xev, yev
% u_kxky_z0_omega_TF_nint = displ.TF.(varname2)(abs(x)<=x_eval_xy,abs(y)<=y_eval_xy,:); %wavenumber domain over xev, yev

%% Apply excitation spectrum P(f) of harmonic/transient load on TF

% Set Im(u_x_y_omega_z0(f=0))=0
u_xy_z0_omega_TF_nint(:,:,f_tot==0) = real(u_xy_z0_omega_TF_nint(:,:,f_tot==0));

%  -> u(x,y,z,f) = TF(x,y,z,f).*P(f) by
%  Multiplication of P(f) in 3rd dim on all displ for each (x,y)

% transient non interpolated and harmonic
    % u(x,y,z,f)
    u_xy_z0_omega_Pf_nint   = u_xy_z0_omega_TF_nint.*Pf;    % with f = f_tot
    u_kxky_z0_omega_Pf_nint = u_kxky_z0_omega_TF_nint.*Pf;  % with f = f_tot
    

%% IFFT omega -> t


        
        %% Initialize
        Nt_harm = dis_itm.Nt_harm;
        
        t_harm  = dis_itm.t_harm;
        f_harm  = dis_itm.f_harm;
        
        Pt_harm    = loading.Pt_harm;
        Pf_harm    = loading.Pf_harm;
        loc_f      = loading.harmonic.loc_f;        % location of f_calc in f_tot
        loc_f_harm = loading.harmonic.loc_f_harm;   % loaction of f_calc in f_harm
        
        %% Frequency/time response u(f), u(t)
        u_xyf_z0_Pf = zeros(length(x),length(y),Nt_harm);
        
        % sort results for f_calc in uf with higher resolution for IFFT
        u_xyf_z0_Pf(:,:,loc_f_harm) = u_xy_z0_omega_Pf_nint(:,:,loc_f);
        u_kxkyf_z0_Pf               = u_kxky_z0_omega_Pf_nint;
        
        % IFFT
        %  Multiply with Nt_harm as Pf not from forward FFT from
        %  Pt, but here directly introduced in f domain
        u_xyt_z0_Pf = Nt_harm*fftshift(ifft(ifftshift(u_xyf_z0_Pf),[],3));
        
        %% Plot: Response P(f),u(f),u(t)
        %-----------------------------------------------------
        mainfigurename = [calc_sys ' - u' displ_dir ' SOIL due to ' load_dir ' - u(x=0,y=0,z=0,f),u(x=0,y=0,z=0,t)'];
        fig = figure('Name',mainfigurename);
        fig.Position = [632,260,1403,1040];
        %-----------------------------------------------------
        sgtitle(['$P_' load_dir(2) '$, TF and $u_{' displ_dir '}$ of SOIL at $x=y=z=0$ for ' calc_sys],'interpreter','latex','FontSize',12)
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
        plot(f_tot,real(squeeze(u_xy_z0_omega_TF_nint(xev==0,yev==0,:))),'Color',plot_cont.color_cyan)
        plot(f_tot,imag(squeeze(u_xy_z0_omega_TF_nint(xev==0,yev==0,:))),'Color',plot_cont.color_orange2)
        plot(f_tot,abs(squeeze(u_xy_z0_omega_TF_nint(xev==0,yev==0,:))),'Color',plot_cont.color_green)
        xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel(['$\mathrm{TF} \;\; u_' displ_dir '(x,y,z,f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        legend('Re','Im','Abs')
        %
        %         subplot(3,2,4)
        %         hold on
        %         Pplot = squeeze(u_x_y_omega_z0(Nx_hs/2+1,Ny_hs/2+1,f_tot>=0));
        %         plot(f_tot(f_tot>=0),abs(atan2(imag(Pplot),real(Pplot))),'Color',plot_cont.color_cyan)
        %         xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        %         ylabel(['$|\mathrm{Phase \; TF}|$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        %
        subplot(3,2,6)
        hold on
        stem(f_harm,real(squeeze(u_xyf_z0_Pf(xev==0,yev==0,:))),'Marker','.')
        stem(f_harm,imag(squeeze(u_xyf_z0_Pf(xev==0,yev==0,:))),'Marker','.')
        xlabel('$f$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel(['$u_' displ_dir '(f)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        legend('Re','Im','Interpreter',plot_cont.interpreter)
        title(['$u_' displ_dir '(x=0,y=0,z=0,f)$'],'Interpreter',plot_cont.interpreter)
        %
        subplot(3,2,5)
        hold on
        plot(t_harm,real(squeeze(u_xyt_z0_Pf(xev==0,yev==0,:))))
        plot(t_harm,imag(squeeze(u_xyt_z0_Pf(xev==0,yev==0,:))))
        xlabel('$t$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel(['$u_' displ_dir '(t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        legend('Re','Im','Interpreter',plot_cont.interpreter)
        title(['$u_' displ_dir '(x=0,y=0,z=0,t)$'],'Interpreter',plot_cont.interpreter)
        %
        savename = ['\displ_soil_P_TF_u_x0y0z0_f_t'];
        saveas(fig,[path.figures savename],'fig')
        %
        if strcmpi(plot_cont.export,'yes')
            exportfig(gcf,[path.figures savename],'Width',plot_cont.exp_width,'Height',plot_cont.exp_heigth,'Fontmode','fixed','Fontsize',plot_cont.Fontsize,'Color','cmyk','resolution',plot_cont.resolution);
        end
        if ~strcmpi(plot_cont.displ_soil,'yes');close(fig);end
        
        %% Plot: Animated plot u(t)
        if strcmpi(plot_cont.u_anim,'yes') && (size(xev,2)>=2 && size(yev,2)>=2)
            %-----------------------------------------------------
            fig          = figure;
            fig.Position = [632,260,1403,1040];
            fig.Name     = [calc_sys ' - u' displ_dir ' SOIL due to ' load_dir];
            %-----------------------------------------------------
            
            %             pause_t = 1/Nt_harm*dt_harm;
            
            sgtitle(['$u_{' displ_dir '}(x=y=z=0)$ SOIL due to ' loading.distr_xy ' ' loading.type ' load $P_' load_dir(2) '(z=0)$'],'interpreter','latex','FontSize',12)
            
            [X,Y] = ndgrid(xev(1:2:end),yev(1:2:end));
            Z     = real(u_xyt_z0_Pf(1:2:end,1:2:end,:));
            Z0    = Z(:,:,1);
            minz  = min(Z,[],'all');
            maxz  = max(Z,[],'all');
            %
            subplot(2,2,1)
            subpl_1 = mesh(X,Y,Z0,'FaceColor',plot_cont.FaceColor);
            view([-22.500000000000004,10.199999999999999])
            zlim([minz maxz])
            colormap(jet)
            caxis([minz/4 maxz/4])
            xlabel('$x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
            ylabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
            zlabel(['$\mathrm{Re} \;\; u_' displ_dir '(x,y,t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
            
            % 2D u(x=0,y,z=0,t)
            subplot(2,2,2);
            u_x0yt = squeeze(real(u_xyt_z0_Pf(xev==0,:,:)));
            %     plot(t_harm,u_x0y0t);
            miny3  = min(u_x0yt,[],'all');
            maxy3  = max(u_x0yt,[],'all');
            subpl_2 = plot(yev,u_x0yt(:,1));
            ylim([miny3 maxy3])
            xlabel('$y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
            ylabel(['$\mathrm{Re} \;\; u_' displ_dir '(x=0,y,t)$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
            
            % 2D: u(t)
            subplot(2,2,3);
            hold on
            %             h1 = animatedline;
            Pplot2 = real(squeeze(u_xyt_z0_Pf(xev==0,yev==0,:)));
            miny1 = min(Pplot2);
            maxy1 = max(Pplot2);
            plot(t_harm,Pplot2)
            s1 = scatter(t_harm(1),Pplot2(1),'MarkerEdgeColor',plot_cont.color_cyan,'MarkerFaceColor',plot_cont.color_cyan);
            axis([t_harm(1) t_harm(end) miny1 maxy1])
            
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
            
            for nnt = 1:1:length(t_harm)
                
                % ISO displ
                subpl_1.ZData = real(u_xyt_z0_Pf(1:2:end,1:2:end,nnt)); % current Data
                
                % 2D u(x=0,y,z=0,t)
                subpl_2.YData = real(squeeze(u_xyt_z0_Pf(xev==0,:,nnt)));
                
                % 2D u(t)
                subplot(2,2,3)
                s1.XData = t_harm(nnt);
                s1.YData = real(squeeze(u_xyt_z0_Pf(xev==0,yev==0,nnt)));
                
                % 2D p(t)
                subplot(2,2,4)
                s2.XData = t_harm(nnt);
                s2.YData = real(Pt(nnt));
                
                del = 0.01;
                frame = getframe(gcf);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                
                savename = ['\displ_soil_u_xyz0_gif'];
                filename = [path.figures savename];
                
                if nnt == 1
                    imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
                end
            end
            
        end
        
   

%% Plot:  TF(f) u(f)
%-----------------------------------------------------
fig          = figure;
fig.Position = [632,260,1403,1040];
fig.Name     = ['hs' '- Resp. soil TF(f) u(f)'];
%-----------------------------------------------------
sgtitle(['$u_{' displ_dir '}(x=y=z=0)$ of SOIL due to ' 'rect' ' ' loading.type ' load $P_' load_dir(2) '(z=0)$'],'interpreter','latex','FontSize',12)
%
subplot(2,2,1)
hold on
plot(f_tot,real(squeeze(u_xy_z0_omega_TF_nint(xev==0,yev==0,:))),'Marker','.','LineStyle','--') % initial -> for interp 'no' and harmonic

xlabel('$f$','Interpreter','latex')
ylabel(['$\mathrm{Re}\;\; \mathrm{TF} \; u_' displ_dir '(f)$'],'Interpreter','latex')
%legend('initial','interpolated','Interpreter','latex')
legend('initial','Interpreter','latex')
title(['Transfer Function $u_' displ_dir '(f)$'],'Interpreter','latex')
%
subplot(2,2,3)
hold on
plot(f_tot,imag(squeeze(u_xy_z0_omega_TF_nint(xev==0,yev==0,:))),'Marker','.','LineStyle','--') % initial

xlabel('$f$','Interpreter','latex')
ylabel(['$\mathrm{Im} \;\; \mathrm{TF} \;  u_' displ_dir '(f)$'],'Interpreter','latex')
%legend('initial','interpolated','Interpreter','latex')
legend('initial','Interpreter','latex')
%
subplot(2,2,2)
hold on

plot(f_tot,real(squeeze((u_xy_z0_omega_Pf_nint(xev==0,yev==0,:))))) % non interpolated

xlabel('$f$','Interpreter','latex')
ylabel(['$\mathrm{Re}\; u_' displ_dir '(f)$'],'Interpreter','latex')
title(['Displ.~spectrum $u_' displ_dir '(f) = \mathrm{TF}(f)\cdot P_' load_dir(2) '(f)$'],'Interpreter','latex')
%
subplot(2,2,4)
hold on

plot(f_tot,imag(squeeze((u_xy_z0_omega_Pf_nint(xev==0,yev==0,:))))) % non interpolated

xlabel('$f$','Interpreter','latex')
ylabel(['$\mathrm{Im}\; u_' displ_dir '(f)$'],'Interpreter','latex')

savename = ['\displ_soil_TF_u_x0y0z0f_orig_int'];
saveas(fig,[path.figures savename],'fig')

if ~strcmpi(plot_cont.displ_soil,'yes');close(fig);end

%% Plot: u(x=x_dist,y=0,z=0,t)
%-----------------------------------------------------------
if ~strcmpi(loading.type,'harmonic')
fig          = figure;
fig.Position = [632,260,1403,1040];
fig.Name     = ['hs' ' - u' displ_dir ' SOIL due to ' load_dir '(x=xev,y=0,z=0,t)'];
%-----------------------------------------------------------
x_ev_pl = 2;
plot(t_trans*mat.hs.cs/x_ev_pl,mat.hs.G*x_ev_pl/loading.P0_itm_hs1L*real(squeeze(u_xyt_z0_Pf(xev==x_ev_pl,yev==0,:))))
set(gca,'Ydir','reverse')
xlabel('$t c_s/x_{ev}$','interpreter','latex')
ylabel(['$u_{' displ_dir '}(x=x_{ev},y=0,z=0,t)$'],'interpreter','latex')
xlim([-2 2])

savename = ['\displ_soil_u_xevy0z0_t_dimless'];
saveas(fig,[path.figures savename],'fig')
end 

%% Plot: u(x,y=0,z=0,t)
%-----------------------------------------------------------
if ~strcmpi(loading.type,'harmonic')
fig          = figure;
fig.Position = [632,260,1403,1040];
fig.Name     = ['hs' ' - u' displ_dir ' SOIL due to ' load_dir '(x,y=0,z=0,t)'];
%-----------------------------------------------------------
[Xev,Ttrans] = ndgrid(xev,t_trans);

mesh(Xev,Ttrans,real(squeeze(u_xyt_z0_Pf(:,yev==0,:))))
xlabel('$x$','interpreter','latex')
ylabel('$t$','interpreter','latex')
zlabel(['$u_{' displ_dir '}(x=x_{ev},y=0,z=0,t)$'],'interpreter','latex')
ylim([-1 2])
% 
savename = ['\displ_soil_u_xy0z0_t_dimless'];
saveas(fig,[path.figures savename],'fig')
end

%% Output

varname1 = [displ_dir load_dir '_xyf_z0'];
varname2 = [displ_dir load_dir '_kxkyf_z0'];
varname3 = [displ_dir load_dir '_xyt_z0'];

displ.Pt.(varname1) = u_xyf_z0_Pf;
displ.Pt.(varname2) = u_kxkyf_z0_Pf;
displ.Pt.(varname3) = u_xyt_z0_Pf;
 
end
