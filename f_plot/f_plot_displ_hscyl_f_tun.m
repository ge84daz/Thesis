%__________________________________________________________________________
%
%        f_plot_displ_hscyl_f_tun.m    
%
%       INPUT: 
%       -displ_dir     = Displ. direction
%       -load_dir      = Displ. due to load on load_dir
%       -displ         = Displacements
%       -dis_itm       = Discretization ITM
%       -dis_fem       = Discretization FEM
%       -geo           = Geometry
%       -loading       = Load (amplitudes, distribution,...)
%       -mat           = Material Parameters
%       -plot_cont     = plot control parameters
%       -path          = path parameters
%       -calc          = Calculation parameters
%       -fev           = Evaluation frequency
%
%
%       DESCRIPTION: 
%       - Select displ at z=0 hs_surf of system hs+tun+FEM with one layer on top
%       - Select displ at z=H_cyl_tot of tunnel midline inside FEM subsystem
%
%       REMARK: 
%
%       - only displ. in frequency domain plotted for negative frequency f
%       - if harmonic time dependent displ shall be plotted: 
%       Time of maximum displ t_max_itm_fem -> equal for trench
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_plot_displ_hscyl_f_tun)
%       Modified:   Tom Hicks                 
%       Date:       31-03-2020 
%       Changed:    02-12-2024 - Hicks   
%__________________________________________________________________________


function [geo_T1] = f_plot_displ_hscyl_f_tun(displ_dir,load_dir,displ,geo_T1,dis_itm,dis_fem_T1,mat,plot_cont,path,calc)

%%  Initialize: Input 
fev   = dis_itm.fev;
loc_f = find(fev == -dis_itm.f);


% Select displ hscyl
u_xyf_hscyl_z0          = getfield(displ.hscyl,calc.system,'itm',[displ_dir load_dir '_xyf_hscyl_itm_z0']);
% Tunnel T1
u_x_midline_f_fem_zH_T1 = getfield(displ.hscyl,calc.system,'fem','T1',[displ_dir load_dir '_xyf_hscyl_fem_midline']);

% Geo & Discretization
x       = geo_T1.x;
y       = geo_T1.y;
Bx      = geo_T1.Bx;
By      = geo_T1.By;
y_Tc_T1 = geo_T1.y_Tc;

Nx_hs          = dis_itm.Nx_hs;
nod_midline_T1 = dis_fem_T1.nod_midline;

% Plot, legend and title names
displ_dir_sub  = strrep(displ_dir,'u','');
leg_name_hscyl = ['ITM-FEM ',strrep(calc.system,'_',' ')];
    
% Normalize x,y  w.r.t. lambda_R of soil
if strcmpi(plot_cont.norm_lambda_R,'yes')
    
    % Rayleigh wavelength of soil for fev
    lambda_R_hs  = mat.hs.cr/fev; 
    
    x            = x./real(lambda_R_hs);
    y            = y./real(lambda_R_hs);  
    y_midline_T1 = nod_midline_T1(:,2)./real(lambda_R_hs);
    x_lim        = (Bx/2)./real(lambda_R_hs);
    y_lim        = (By/2)./real(lambda_R_hs);
    
else
    y_midline_T1 = nod_midline_T1(:,2) + y_Tc_T1;
    x_lim     = (Bx/2);
    y_lim     = (By/2);
end

[X,Y] = ndgrid(x,y);

%% Plots
if strcmpi(plot_cont.displ_comp,'yes')
         
    %% Re  u(x,y,z=0,omega) - ItmFem
    %--------------------------------------------------------------------------
    fig                     = figure;
    fig.Name                = ['ItmFem - Re u(x,y,z=0,f= ' num2str(fev) ') - u' displ_dir_sub ' ' load_dir];
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------------------------------------
    %
    mesh(X,Y,real(u_xyf_hscyl_z0(:,:,loc_f)),'FaceColor',plot_cont.FaceColor)
    
    switch lower(plot_cont.map),     case 'yes', colormap(plot_cont.colormap),   otherwise, end
    switch lower(plot_cont.shading), case 'yes', shading(plot_cont.shading_typ), otherwise, end
    
    xlim([-x_lim +x_lim])
    ylim([-y_lim +y_lim])
    
    switch lower(plot_cont.norm_lambda_R)
        case 'yes'; xlabel('$x/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter),...
                    ylabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$x$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter),...
                    ylabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    zlabel(['$\mathrm{Re}\;u_{' displ_dir_sub ' ' load_dir '}(x,y)$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(leg_name_hscyl);
        lgd.Location    = plot_cont.loc_ledg;
        lgd.Box         = plot_cont.leg_box;
        lgd.FontSize    = plot_cont.Fontsize;
        lgd.Interpreter = plot_cont.interpreter;
    end
    
    savename = ['\displ_hscyl_fev_' num2str(fev) '_' displ_dir load_dir '_xyf_z0_re'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
    
    %% Abs  u(x,y,z=0,omega) - ItmFem
    %--------------------------------------------------------------------------
    fig                     = figure;
    fig.Name                = ['ItmFem - Abs u(x,y,z=0,f= ' num2str(fev) ') - u' displ_dir_sub ' ' load_dir];
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------------------------------------
    %
    mesh(X,Y,abs(u_xyf_hscyl_z0(:,:,loc_f)),'FaceColor',plot_cont.FaceColor)
    
    switch lower(plot_cont.map),     case 'yes', colormap(plot_cont.colormap),   otherwise, end
    switch lower(plot_cont.shading), case 'yes', shading(plot_cont.shading_typ), otherwise, end
    
    xlim([-x_lim +x_lim])
    ylim([-y_lim +y_lim])
    
    switch lower(plot_cont.norm_lambda_R)
        case 'yes'; xlabel('$x/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter),...
                    ylabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$x$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter),...
                    ylabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    zlabel(['$|\;u_{' displ_dir_sub ' ' load_dir '}(x,y)$|'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(leg_name_hscyl);
        lgd.Location    = plot_cont.loc_ledg;
        lgd.Box         = plot_cont.leg_box;
        lgd.FontSize    = plot_cont.Fontsize;
        lgd.Interpreter = plot_cont.interpreter;
    end
    
    savename = ['\displ_hscyl_fev_' num2str(fev) '_' displ_dir load_dir '_xyf_z0_abs'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end, otherwise, end
      
    %% Re  u(x=0,y,z=0,omega)
    %--------------------------------------------------------------------------
    fig                     = figure;
    fig.Name                = ['Comp - u(x=0,y,z=0,f= ' num2str(fev) ') - u' displ_dir_sub ' ' load_dir];
    fig.Position            = [500,250,1800,800]; 
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------------------------------------
    hold on
    %
    subplot(1,3,1)
    %
    plot(y,real(u_xyf_hscyl_z0(Nx_hs/2+1,:,loc_f)),'-','Color',plot_cont.color_cyan,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth ) % Orange        -> ITM-FEM
    
    switch lower(plot_cont.norm_lambda_R)
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$\mathrm{Re}\;(u_{' displ_dir_sub ' ' load_dir '}(y))$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(leg_name_hscyl);
        lgd.Location    = plot_cont.loc_ledg;
        lgd.Box         = plot_cont.leg_box;
        lgd.FontSize    = plot_cont.Fontsize;
        lgd.Interpreter = plot_cont.interpreter;
    end
    
    % Set Figure and axis handles
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))                 % Set small figure surrounding
    set(gca,'FontSize',plot_cont.Fontsize)                                  % Fontsize axis
    set(findall(gca, 'Type', 'Line'),'LineWidth',plot_cont.LineWidth);      % Linewidth all lines
    set(gca,'LineWidth',plot_cont.LineWidth)                                % Linewidth axis
    
    savename = ['\displ_hscyl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_z0_re'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end, otherwise, end
    
    %% Im  u(x=0,y,z=0,omega)
    subplot(1,3,2)
    %
    plot(y,imag(u_xyf_hscyl_z0(Nx_hs/2+1,:,loc_f)),'-','Color',plot_cont.color_cyan,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth ) % Orange        -> ITM-FEM
    
    switch lower(plot_cont.norm_lambda_R)
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$\mathrm{Im}\;(u_{' displ_dir_sub ' ' load_dir '}(y))$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(leg_name_hscyl);
        lgd.Location    = plot_cont.loc_ledg;
        lgd.Box         = plot_cont.leg_box;
        lgd.FontSize    = plot_cont.Fontsize;
        lgd.Interpreter = plot_cont.interpreter;
    end
    
    % Set Figure and axis handles
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))                 % Set small figure surrounding
    set(gca,'FontSize',plot_cont.Fontsize)                                  % Fontsize axis
    set(findall(gca, 'Type', 'Line'),'LineWidth',plot_cont.LineWidth);      % Linewidth all lines
    set(gca,'LineWidth',plot_cont.LineWidth)                                % Linewidth axis
    
    savename = ['\displ_hscyl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_z0_im'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end, otherwise, end
    
    %% Abs u(x=0,y,z=0,omega)
    subplot(1,3,3)
    %
    plot(y,abs(u_xyf_hscyl_z0(Nx_hs/2+1,:,loc_f)),'-','Color',plot_cont.color_cyan,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth ) % Orange        -> ITM-FEM
    
    switch lower(plot_cont.norm_lambda_R)
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$|\;(u_{' displ_dir_sub ' ' load_dir '}(y))|$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(leg_name_hscyl);
        lgd.Location    = plot_cont.loc_ledg;
        lgd.Box         = plot_cont.leg_box;
        lgd.FontSize    = plot_cont.Fontsize;
        lgd.Interpreter = plot_cont.interpreter;
    end
    
    % Set Figure and axis handles
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))                 % Set small figure surrounding
    set(gca,'FontSize',plot_cont.Fontsize)                                  % Fontsize axis
    set(findall(gca, 'Type', 'Line'),'LineWidth',plot_cont.LineWidth);      % Linewidth all lines
    set(gca,'LineWidth',plot_cont.LineWidth)                                % Linewidth axis
    
    savename = ['\displ_hscyl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_z0_abs'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
    
    %% Re  u(x=0,y,z=H,omega) 
    % at tunnel midline over y
    %--------------------------------------------------------------------------
    fig                     = figure;
    fig.Name                = ['u(x=0,y,z=H,f= ' num2str(fev) ') - u' displ_dir_sub ' ' load_dir];
    fig.Position            = plot_cont.fig_position;
    fig.PaperPositionMode   = 'auto';
    %--------------------------------------------------------------------------
    hold on
    %
    subplot(1,3,1)
    %
    
    plot(y_midline_T1, real(u_x_midline_f_fem_zH_T1(Nx_hs/2+1,:,loc_f))','-','Color',plot_cont.color_cyan,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)

    switch lower(plot_cont.norm_lambda_R)
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$\mathrm{Re}\;(u_{' displ_dir_sub ' ' load_dir '_s}(y))$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(leg_name_hscyl);
        lgd.Location    = plot_cont.loc_ledg;
        lgd.Box         = plot_cont.leg_box;
        lgd.FontSize    = plot_cont.Fontsize;
        lgd.Interpreter = plot_cont.interpreter;
    end
    
    % Set Figure and axis handles
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))                 % Set small figure surrounding
    set(gca,'FontSize',plot_cont.Fontsize)                                  % Fontsize axis
    set(findall(gca, 'Type', 'Line'),'LineWidth',plot_cont.LineWidth);      % Linewidth all lines
    set(gca,'LineWidth',plot_cont.LineWidth)
    
    savename = ['\displ_hscyl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_zH_re'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end

    %% Im  u(x=0,y,z=H,omega) 
    % at tunnel midline over y
    subplot(1,3,2)
    
    plot(y_midline_T1, imag(u_x_midline_f_fem_zH_T1(Nx_hs/2+1,:,loc_f))','-','Color',plot_cont.color_cyan,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)
    switch lower(plot_cont.norm_lambda_R)
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$\mathrm{Im}\;(u_{' displ_dir_sub ' ' load_dir '_s}(y))$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(leg_name_hscyl);
        lgd.Location    = plot_cont.loc_ledg;
        lgd.Box         = plot_cont.leg_box;
        lgd.FontSize    = plot_cont.Fontsize;
        lgd.Interpreter = plot_cont.interpreter;
    end
    
    % Set Figure and axis handles
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))                 % Set small figure surrounding
    set(gca,'FontSize',plot_cont.Fontsize)                                  % Fontsize axis
    set(findall(gca, 'Type', 'Line'),'LineWidth',plot_cont.LineWidth);      % Linewidth all lines
    set(gca,'LineWidth',plot_cont.LineWidth)
    
    savename = ['\displ_hscyl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_zH_im'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
    
    %% Abs u(x=0,y,z=H,omega) 
    % at tunnel midline over y
    subplot(1,3,3)
    
    plot(y_midline_T1, abs(u_x_midline_f_fem_zH_T1(Nx_hs/2+1,:,loc_f))','-','Color',plot_cont.color_cyan,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)
    switch lower(plot_cont.norm_lambda_R)
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$|\;(u_{' displ_dir_sub ' ' load_dir '_s}(y))|$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(leg_name_hscyl);
        lgd.Location    = plot_cont.loc_ledg;
        lgd.Box         = plot_cont.leg_box;
        lgd.FontSize    = plot_cont.Fontsize;
        lgd.Interpreter = plot_cont.interpreter;
    end
    
    % Set Figure and axis handles
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))                 % Set small figure surrounding
    set(gca,'FontSize',plot_cont.Fontsize)                                  % Fontsize axis
    set(findall(gca, 'Type', 'Line'),'LineWidth',plot_cont.LineWidth);      % Linewidth all lines
    set(gca,'LineWidth',plot_cont.LineWidth)
    
    savename = ['\displ_hscyl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_zH_abs'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
    
 end

%%  Output

if strcmpi(plot_cont.norm_lambda_R,'yes')
    geo_T1.x_norm         = x;
    geo_T1.y_norm         = y;
    geo_T1.y_midline_norm = y_midline_T1;
    geo_T1.x_lim_norm     = x_lim;
    geo_T1.y_lim_norm     = y_lim;
end

end
