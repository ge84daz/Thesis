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



function [] = f_plot_ampl_hs_surf(dis_itm,geo,displ,plot_cont,path)

%% 1.) Initialize: Input

% f_tot     = dis_itm.f_tot;
geo.x_eval_xy = 10;
geo.y_eval_xy = 10;

x_eval_xy = geo.x_eval_xy;
y_eval_xy = geo.y_eval_xy;

kx       = dis_itm.kx;
ky       = dis_itm.ky;

Nx = dis_itm.Nx_hs;
Ny = dis_itm.Ny_hs;

x           = geo.x;
y           = geo.y;

xev         = geo.x(abs(geo.x)<=x_eval_xy);
yev         = geo.y(abs(geo.y)<=y_eval_xy);

kx_ev    = dis_itm.kx(abs(geo.x)<=x_eval_xy);
ky_ev    = dis_itm.ky(abs(geo.y)<=x_eval_xy);

% coefficients only for first negative omega
A2_neg =  displ.TF.C_A2_z0(:,:,1);
Bx2_neg = displ.TF.C_Bx2_z0(:,:,1);
By2_neg = displ.TF.C_By2_z0(:,:,1);

dkxky   = 6; % for low memory plots of total domain
A2_neg_ev =  displ.TF.C_A2_z0(abs(geo.x)<=x_eval_xy,abs(geo.y)<=y_eval_xy,1);
Bx2_neg_ev =  displ.TF.C_Bx2_z0(abs(geo.x)<=x_eval_xy,abs(geo.y)<=y_eval_xy,1);
By2_neg_ev =  displ.TF.C_By2_z0(abs(geo.x)<=x_eval_xy,abs(geo.y)<=y_eval_xy,1);

%% 2D plot over kx or ky
mainfigurename = 'Wave amplitude A2 over ky';
fig = figure('Name',mainfigurename);
fig.Position = [632,260,1403,1040];
plot(ky, abs(A2_neg(Nx/2+1,:)));
xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% plot(ky, abs(A2_neg(Nx/2,:)));
% xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('Wave amplitudes $A_2$');
% xlim([-12, 12]);
hold on
sgtitle('Wave amplitude A2 over ky','Interpreter',plot_cont.interpreter)
savename = '\amplitude_A2_ky';
saveas(fig,[path.figures savename],'fig')

mainfigurename = 'Wave amplitude Bx2 over ky';
fig = figure('Name',mainfigurename);
fig.Position = [632,260,1403,1040];
plot(ky, abs(Bx2_neg(Nx/2+1,:)));
xlabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('Wave amplitudes $B_{x2}$');
hold on
%xlim([-50, 50]);
sgtitle('Wave amplitude Bx2 over ky','Interpreter',plot_cont.interpreter)
savename = '\amplitude_Bx2_ky';
saveas(fig,[path.figures savename],'fig')

mainfigurename = 'Wave amplitude By2 over kx';
fig = figure('Name',mainfigurename);
fig.Position = [632,260,1403,1040];
plot(kx, abs(By2_neg(:,Ny/2+1)));
xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('Wave amplitudes $B_{y2}$');
hold on
sgtitle('Wave amplitude By2 over kx','Interpreter',plot_cont.interpreter)
savename = '\amplitude_By2_kx';
saveas(fig,[path.figures savename],'fig')

%% 3D plot over kx or ky
mainfigurename = 'Wave amplitudes A2, Bx2, By2 over kx-ky';
fig = figure('Name',mainfigurename);
fig.Position = [632,260,1403,1040];
% for A2
subplot(3,2,1)
% mesh(kx, ky, abs(A2_neg)); %too high memory
mesh(kx_ev, ky_ev, abs(A2_neg_ev),'FaceColor',plot_cont.FaceColor);
xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
zlabel('A2');
hold on
% low memory plot - only each dkxky point is plotted
subplot(3,2,2)
[KX,KY] = ndgrid(kx(1:dkxky:end),ky(1:dkxky:end));
mesh(KX,KY,abs(A2_neg(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
zlabel('A2');
hold on

% for Bx2
subplot(3,2,3)
mesh(kx_ev, ky_ev, abs(Bx2_neg_ev),'FaceColor',plot_cont.FaceColor);
xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
zlabel('Bx2');
hold on
% low memory plot - only each dkxky point is plotted
subplot(3,2,4)
[KX,KY] = ndgrid(kx(1:dkxky:end),ky(1:dkxky:end));
mesh(KX,KY,abs(Bx2_neg(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
zlabel('Bx2');
hold on

% for By2
subplot(3,2,5)
mesh(kx_ev, ky_ev, abs(By2_neg_ev),'FaceColor',plot_cont.FaceColor);
xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
zlabel('By2');
hold on
% low memory plot - only each dkxky point is plotted
subplot(3,2,6)
[KX,KY] = ndgrid(kx(1:dkxky:end),ky(1:dkxky:end));
mesh(KX,KY,abs(By2_neg(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
zlabel('By2');
hold on

%save plots
sgtitle('Wave amplitudes over kx-ky','Interpreter',plot_cont.interpreter)
savename = '\amplitude_A2_Bx2By2_kxky';
saveas(fig,[path.figures savename],'fig')


%% old plots
% subplot(3,1,2)
% %mesh(kx, ky, abs(Bx2_neg));
% % low memory
% mesh(KX,KY,abs(Bx2_neg(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('Wave amplitudes Bx2');
% hold on
% % xlim([-2, 2]);
% % ylim([-2, 2]);
% 
% subplot(3,1,3)
% %mesh(kx, ky, abs(By2_neg));
% % low memory
% mesh(KX,KY,abs(By2_neg(1:dkxky:end,1:dkxky:end)),'FaceColor',plot_cont.FaceColor);
% xlabel('$k_x$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% ylabel('$k_y$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter);
% zlabel('Wave amplitudes By2');
% hold on
% % xlim([-2, 2]);
% % ylim([-2, 2]);
% sgtitle('Wave amplitudes over kx-ky','Interpreter',plot_cont.interpreter)
% %savename = '\amplitude_A2_Bx2By2_kxky';
% %saveas(fig,[path.figures savename],'fig')














end