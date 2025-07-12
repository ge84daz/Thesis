%__________________________________________________________________________
%
%        f_plot_displ_hscyl_f_tr_sl.m    
%
%       INPUT: 
%       displ_dir     = Displ. direction
%       load_dir      = Displ. due to load on load_dir
%       displ         = Displacements
%       dis_itm       = Discretization ITM
%       dis_fem       = Discretization FEM
%       geo           = Geometry
%       loading       = Load (amplitudes, distribution,...)
%       mat           = Material Parameters
%       plot_cont     = plot control parameters
%       path          = path parameters
%       calc          = Calculation parameters
%       fev           = Evaluation frequency
%
%
%       DESCRIPTION: 
%       The function combines the solution parts on the halfspace surface of the
%       ITM and FEM System for a halfspace with a cylindrical trench/slit and
%       displayes them in different plots, also compared to a reference solution
%       calculated using a layered halfspace with one layer.
%       Furthermore the Amplitude Reduction Factors AR are calculated for
%       different setups.
%
%       REMARK: 
%       Original Author: Julian Freisinger (f_plot_displ_hscyl_f_tr_sl)
%       Modified:   Tom Hicks                 
%       Date:       31-03-2020 
%       Changed:    02-12-2024 - Hicks   
%__________________________________________________________________________

function [displ,geo_T1] = f_plot_displ_hscyl_f_tr_sl(displ_dir,load_dir,displ,geo_T1,dis_itm,~,~,mat,plot_cont,path,calc)

%% Initialize: Input

% index of fev within f
fev   = dis_itm.fev;
loc_f = find(fev == -dis_itm.f);

% Select dipl hscyl
u_xyf_comb_tr_z0 = getfield(displ.hscyl,'trench','itmfem',[displ_dir load_dir '_xyf_comb_tr_z0']);

% Geo
x       = geo_T1.x;
y       = geo_T1.y;
y_comb  = geo_T1.y_comb;

lambda_R_hs = mat.hs.cr/fev; % Rayleigh wavelength of soil


% Tunnel 1
y_Tc_T1        = geo_T1.y_Tc;

% Plot names
name_ItmFem_trench  = 'ITM-FEM - tr';    %'Coupled approach - trench'
displ_dir_sub       = strrep(displ_dir,'u','');

switch calc.system
    case 'trench';  lgd_names = {name_ItmFem_trench};
end

%% Combination of the solution parts ITM-FEM-ITM - Displ
% ITM(x<r) - FEM nod_middline_in - ITM(x>r)

% Select f=fev and nnPiF
u_x0yfev_comb_tr_z0  = squeeze(u_xyf_comb_tr_z0(x==0,:,loc_f));
u_xyfev_comb_tr_z0   = squeeze(u_xyf_comb_tr_z0(:,:,loc_f));
u_xyTcfev_comb_tr_z0 = squeeze(u_xyf_comb_tr_z0(:,y_comb==y_Tc_T1,loc_f));


%% Normalize x,y w.r.t. lambda_R of soil
if strcmpi(plot_cont.norm_lambda_R,'yes')
    x         = x./real(lambda_R_hs);
    y         = y./real(lambda_R_hs);
    y_comb    = y_comb./real(lambda_R_hs);
end

[X_comb,Y_comb] = ndgrid(x,y_comb);

%% Plots 

if strcmpi(plot_cont.displ_comp,'yes')
    %%  Trench - Re u(x,y,z=0,fev)
    if strcmpi(calc.system,'trench') 
        %----------------------------------------------------------------------
        fig          = figure;
        fig.Name     = ['Re u(x,y,z=0,f=fev=) Trench - u' displ_dir_sub ' ' load_dir];
        fig.Position = plot_cont.fig_position;
        %----------------------------------------------------------------------
        
        mesh(X_comb,Y_comb,real(u_xyfev_comb_tr_z0),'FaceColor',plot_cont.FaceColor) % coupled ITM-FEM -> trench
        
        switch plot_cont.map,     case 'yes', colormap(plot_cont.colormap),   otherwise, end
        switch plot_cont.shading, case 'yes', shading(plot_cont.shading_typ), otherwise, end
        
        switch plot_cont.norm_lambda_R
            case 'yes'; xlabel('$x/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter),...
                    ylabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
            case 'no';  xlabel('$x$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter),...
                    ylabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        end
        zlabel(['$\mathrm{Re}\;(u_' displ_dir_sub ' ' load_dir '(x,y))$'], 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        title(['ItmFem - Trench $f_{ev}$=' num2str(fev) ],'Interpreter',plot_cont.interpreter)
        
        savename = ['\displ_hscyl_tr_fev_' num2str(fev) '_' displ_dir load_dir '_xyf_z0_re'];
        saveas(fig,[path.figures savename],'fig')
        %
        if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
        switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
        
    end
    
     
    %%  Comp   - Abs u(x=0,y,z=0,omega)
    %----------------------------------------------------------------------
    fig          = figure;
    fig.Name     = ['u(x=0,y,z=0,-fev=' num2str(fev) ') - u' displ_dir_sub ' ' load_dir];
    fig.Position = [500,250,1800,800]; 
    %----------------------------------------------------------------------
    hold on
    subplot(1,3,1)
    
    plot(y_comb, abs(u_x0yfev_comb_tr_z0), '-','Color',plot_cont.color_blue,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)
   
    % Label
    switch plot_cont.norm_lambda_R
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$|u_' displ_dir_sub ' ' load_dir '(y)|$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    % Legend
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(lgd_names);
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
    
    savename = ['\displ_hscyl_trsl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_z0_abs'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
    
    %%  Comp   - Re u(x=0,y,z=0,omega)
    subplot(1,3,2)   
    
    plot(y_comb, real(u_x0yfev_comb_tr_z0), '-','Color',plot_cont.color_blue,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)

    % Label
    switch plot_cont.norm_lambda_R
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$\mathrm{Re} \; u_' displ_dir_sub ' ' load_dir '(y)$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    % Legend
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(lgd_names);
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
    
    savename = ['\displ_hscyl_trsl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_z0_re'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
    
    %%  Comp   - Im u(x=0,y,z=0,omega)
    subplot(1,3,3)  
    
    plot(y_comb, imag(u_x0yfev_comb_tr_z0), '-','Color',plot_cont.color_blue,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)

    
    % Label
    switch plot_cont.norm_lambda_R
        case 'yes'; xlabel('$y/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$y$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$\mathrm{Im} \; u_' displ_dir_sub ' ' load_dir '(y)$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    % Legend
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(lgd_names);
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
    
    savename = ['\displ_hscyl_trsl_fev_' num2str(fev) '_' displ_dir load_dir '_x0yf_z0_im'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
     
    %%  Comp   - Abs u(x,y=yTc,z=0,omega)
    %----------------------------------------------------------------------
    fig          = figure;
    fig.Name     = ['Abs u(x,y=' num2str(y_Tc_T1) ',z=0,-fev=' num2str(fev) ') - u' displ_dir_sub ' ' load_dir];
    fig.Position = [500,250,1800,800]; 
    %----------------------------------------------------------------------
    hold on
    subplot(1,3,1)
    
    plot(x, abs(u_xyTcfev_comb_tr_z0), '-','Color',plot_cont.color_blue,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)
    
    % Label
    switch plot_cont.norm_lambda_R
        case 'yes'; xlabel('$x/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$x$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$|u_' displ_dir_sub ' ' load_dir '(x)|$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    % Legend
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(lgd_names);
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
    
    savename = ['\displ_hscyl_trsl_fev_' num2str(fev) '_' displ_dir load_dir '_xyTcf_z0_abs'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
       
    %%  Comp   - Re u(x,y=yTc,z=0,omega)
    subplot(1,3,2)
    
    plot(x, real(u_xyTcfev_comb_tr_z0), '-','Color',plot_cont.color_blue,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)
    
    % Label
    switch plot_cont.norm_lambda_R
        case 'yes'; xlabel('$x/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$x$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$\mathrm{Re} \; u_' displ_dir_sub ' ' load_dir '(x)$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    % Legend
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(lgd_names);
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
    
    savename = ['\displ_hscyl_trsl_fev_' num2str(fev) '_' displ_dir load_dir '_xyTcf_z0_re'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
    
    %%  Comp   - Im u(x,y=yTc,z=0,omega)
    subplot(1,3,3)
    
    plot(x, imag(u_xyTcfev_comb_tr_z0), '-','Color',plot_cont.color_blue,'MarkerSize',plot_cont.MarkerSize,'LineWidth',plot_cont.LineWidth)
    
    % Label
    switch plot_cont.norm_lambda_R
        case 'yes'; xlabel('$x/\lambda_R$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
        case 'no';  xlabel('$x$','FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
    end
    ylabel(['$\mathrm{Im} \; u_' displ_dir_sub ' ' load_dir '(x)$'],'FontSize',plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
    
    % Legend
    if strcmpi(plot_cont.leg,'yes')
        lgd             = legend(lgd_names);
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
    
    savename = ['\displ_hscyl_trsl_fev_' num2str(fev) '_' displ_dir load_dir '_xyTcf_z0_im'];
    saveas(fig,[path.figures savename],'fig')
    %
    if  strcmpi(plot_cont.export,'yes');exportfig(gcf,[path.figures savename],plot_cont.exportfig{:});end
    switch plot_cont.closefigures, case 'yes', close(gcf), otherwise, end
      
end











