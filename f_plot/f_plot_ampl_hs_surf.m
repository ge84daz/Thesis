%__________________________________________________________________________
%
%       f_plot_ampl_hs_surf.m    
%
%       INPUT: 
%       displ       = Displ. hs_cyl in (kx,ky,z,omega) in z-direction
%       dis_itm     = ITM discretization
%
%       OUTPUT: 
%
%       DESCRIPTION: 
%       - Plot of one impedance at z0 original kx-ky grid
%       - Plot of displ at z0
%       - in image space over different sampling (kx-ky)
%       - no index: original sampling
%       - tot: original sampling with increased sampling within [-kb kb]
%
%       REMARK: 
%       Original Author:  Tom Hicks      
%       Date:       30-08-2023 
%       Changed:    25-09-2023 - Hicks   
%__________________________________________________________________________



function [] = f_plot_ampl_hs_surf(dis_itm,dis_fem_T1,geo_T1,displ,plot_cont,path,calc)

% 1.) Initialize: Input

% get frequency
f_tot     = dis_itm.f_tot;
loc_f     = 1;

% original sampling
kx       = dis_itm.kx;
ky       = dis_itm.ky;

Nx_hs            = dis_itm.Nx_hs; 
Ny_hs            = dis_itm.Ny_hs; 


varname = {'uPx_kxkyf_hscyl_tot','uPy_kxkyf_hscyl_tot','uPz_kxkyf_hscyl_tot'};
% get all displacements of all nodes due to given load case
uPz_kxkyf        = displ.hscyl.(calc.system).tot.(varname{3});
% get total stiffness matrix for all nodes hs surf + cyl bound + cyl interior
K_GES_hs_cyl_tun = displ.hscyl.(calc.system).tot.K_GES_kx0;

% get displacement uz for hs surf (z=0) due to load Pz
uz_kxkyf         =  uPz_kxkyf(Nx_hs/2+1,3:3:3*Ny_hs);

% get impedance in z-direction at hs surf due to unit load in z-direction at hs surf
K_zz_z0    = K_GES_hs_cyl_tun(3:3:3*Ny_hs, 3:3:3*Ny_hs);
Z_zz_z0    = inv(K_zz_z0); 
Z_zz       = diag(Z_zz_z0);

% get impedance in z-direction at cyl bound due to unit load in z-direction at hs surf
K_zz_z0    = K_GES_hs_cyl_tun(3:3:3*Ny_hs, 3:3:3*Ny_hs);
Z_zz_z0    = inv(K_zz_z0); 
Z_zz       = diag(Z_zz_z0);





% user input until what wavenumber the plots are evaluated
% dkxky     = dis_itm.dkxky;           % for low memory plots of total domain
% kb        = dis_itm.kb;              % boundary wavenumber for fine and coarse sampling
% focus on range [-kb kb]
% kx_eval           = dis_itm.kx(abs(dis_itm.kx) <= kb);
% ky_eval           = dis_itm.ky(abs(dis_itm.ky) <= kb);
% focus on range [-kb kb]
% Z_zz_eval = Z_zz(abs(kx)<=kb,abs(ky)<=kb);


% displacement at z0
% uz_kxkyf       =  displ.TF.uzPz_kxkyf_z0;
% uz_kxkyf       =  displ.TF.uzPz_kxkyf_zH_L1;
% uz_kxkyf       =  displ.TF.uzPz_kxkyf_zH_L2;


%% Complex stiffness plots - original grid
% 2D plot over ky
mainfigurename = 'Impedance and displacement in image space';
fig = figure('Name',mainfigurename);
fig.Position = [632,260,1403,1040];
%
subplot(3,2,1)
plot(ky, abs(Z_zz(:)));
xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$|Z_{dyn,zz}(k_x=0,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
%
subplot(3,2,3)
plot(ky, real(Z_zz(:)));
xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$\mathrm{Re}\;\; Z_{dyn,zz}(k_x=0,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
%
subplot(3,2,5)
plot(ky, imag(Z_zz(:)));
xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$\mathrm{Im}\;\; Z_{dyn,zz}(k_x=0,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)

%% displacement plots in image space - original grid
subplot(3,2,2)
plot(ky,abs(uz_kxkyf),'Color',plot_cont.color_cyan);
xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$|u_z(k_x=0,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
title(['$u_{z,Pz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
%
subplot(3,2,4)
plot(ky,real(uz_kxkyf),'Color',plot_cont.color_cyan);
xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$\mathrm{Re}\;\; u_z(k_x=0,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
title(['$u_{z,Pz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
    %
subplot(3,2,6)
plot(ky,imag(uz_kxkyf),'Color',plot_cont.color_cyan);
xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$\mathrm{Im}\;\; u_z(k_x=0,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
title(['$u_{z,Pz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% 
% % export
% if strcmp(plot_cont.export, 'yes')                  
%        savename = '\impedance and displacement - original grid';
%        saveas(fig, fullfile(path.figures, savename));
% end
% if strcmp(plot_cont.closefigures, 'yes')                     % close figure after export      
%           close(fig);
% end

% %% Complex stiffness plots - fine grid (tot)
% if strcmpi(calc.tot,'yes')
% % 2D plot over ky
% mainfigurename = 'Impedance in image space - fine grid';
% fig = figure('Name',mainfigurename);
% fig.Position = [632,260,1403,1040];
% %
% subplot(3,2,1)
% plot(ky_tot, abs(Z_zz_tot(kx_tot_zero,:)));
% xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$|Z_{dyn,zz}(k_x=0,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,3)
% plot(ky_tot, real(Z_zz_tot(kx_tot_zero,:,loc_f(1))));
% xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$\mathrm{Re}\;\; Z_{dyn,zz}(k_x=0,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,5)
% plot(ky_tot, imag(Z_zz_tot(kx_tot_zero,:,loc_f(1))));
% xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$\mathrm{Im}\;\; Z_{dyn,zz}(k_x=0,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,2)
% plot(ky_tot_eval, abs(Z_zz_tot_eval(kx_tot_eval_zero,:,loc_f(1))));
% xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$|Z_{dyn,zz}(k_x=0,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,4)
% plot(ky_tot_eval, real(Z_zz_tot_eval(kx_tot_eval_zero,:,loc_f(1))));
% xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$\mathrm{Re}\;\; Z_{dyn,zz}(k_x=0,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,6)
% plot(ky_tot_eval, imag(Z_zz_tot_eval(kx_tot_eval_zero,:,loc_f(1))));
% xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$\mathrm{Im}\;\; Z_{dyn,zz}(k_x=0,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x=0,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% 
% % export
% if strcmp(plot_cont.export, 'yes')                  
%        savename = '\impedance - fine grid';
%        saveas(fig, fullfile(path.figures, savename));
% end
% if strcmp(plot_cont.closefigures, 'yes')                     % close figure after export      
%           close(fig);
% end
% 
% end
% %% 3D plot Kdyn over kx-ky plane - original grid
% if strcmpi(plot_cont.lowMemory,'no')
% mainfigurename = 'Dynamic stiffness over kx-ky - original';
% fig = figure('Name',mainfigurename);
% fig.Position = [632,260,1403,1040];
% %
% subplot(3,2,1)
% [KX_eval,KY_eval] = ndgrid(kx_eval,ky_eval);
% mesh(KX_eval, KY_eval, abs(Z_zz_eval(:,:)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$|Z_{dyn,zz}(k_x,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,2)
% [KX,KY] = ndgrid(kx(1:dkxky:end),ky(1:dkxky:end));
% mesh(KX,KY,abs(Z_zz(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$|Z_{dyn,zz}(k_x,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,3)
% mesh(KX_eval, KY_eval, real(Z_zz_eval(:,:)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$\mathrm{Re}\;\; Z_{dyn,zz}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,4)
% mesh(KX,KY,real(Z_zz(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$\mathrm{Re}\;\; Z_{dyn,zz}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,5)
% mesh(KX_eval, KY_eval, imag(Z_zz_eval(:,:)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$\mathrm{Im}\;\; Z_{dyn,zz}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,6)
% mesh(KX,KY,imag(Z_zz(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$\mathrm{Im}\;\; Z_{dyn,zz}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% 
% % export
% if strcmp(plot_cont.export, 'yes')                  
%        savename = '\3D impedance - original grid';
%        saveas(fig, fullfile(path.figures, savename));
% end
% if strcmp(plot_cont.closefigures, 'yes')                     % close figure after export      
%           close(fig);
% end
% 
% end
% 
% %% 3D plot Kdyn over kx-ky plane - fine grid (tot)
% if strcmpi(plot_cont.lowMemory,'no')
% if strcmpi(calc.tot,'yes')
% mainfigurename = 'Dynamic stiffness over kx-ky - fine grid';
% fig = figure('Name',mainfigurename);
% fig.Position = [632,260,1403,1040];
% %
% subplot(3,2,1)
% [KX_tot_eval,KY_tot_eval] = ndgrid(kx_tot_eval,ky_tot_eval);
% mesh(KX_tot_eval, KY_tot_eval, abs(Z_zz_tot_eval(:,:)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$|Z_{dyn,zz}(k_x,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,2)
% [KX_tot,KY_tot] = ndgrid(kx_tot(1:dkxky:end),ky_tot(1:dkxky:end));
% mesh(KX_tot,KY_tot,abs(Z_zz_tot(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$|Z_{dyn,zz}(k_x,k_y)|$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,3)
% mesh(KX_tot_eval, KY_tot_eval, real(Z_zz_tot_eval(:,:)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$\mathrm{Re}\;\; Z_{dyn,zz}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,4)
% mesh(KX_tot,KY_tot,real(Z_zz_tot(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$\mathrm{Re}\;\; Z_{dyn,zz}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,5)
% mesh(KX_tot_eval, KY_tot_eval, imag(Z_zz_tot_eval(:,:)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$\mathrm{Im}\;\; Z_{dyn,zz}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% %
% subplot(3,2,6)
% mesh(KX_tot,KY_tot,imag(Z_zz_tot(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('$\mathrm{Im}\;\; Z_{dyn,zz}(k_x,k_y)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% title(['$Z_{dyn,zz}(k_x,k_y,z=0,-f=' num2str(f_tot(loc_f(1))) ')$'],'Interpreter',plot_cont.interpreter)
% 
% % export
% if strcmp(plot_cont.export, 'yes')                  
%        savename = '\3D impedance - fine grid';
%        saveas(fig, fullfile(path.figures, savename));
% end
% if strcmp(plot_cont.closefigures, 'yes')                     % close figure after export      
%           close(fig);
% end
% end
% end

end % end function