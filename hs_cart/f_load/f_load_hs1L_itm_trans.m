%__________________________________________________________________________
%
%       f_load_hs1L_itm_transF.m    
%
%       INPUT: 
%
%       OUTPUT: 
%
%       DESCRIPTION: 
%
%       REMARK: 
%       Original Author: Julian Freisinger  (f_load_hs1L_itm_trans)
%       Modified:   Tom Hicks                 
%       Date:       23-05-2023 
%       Changed:    23-05-2023 - Hicks   
%__________________________________________________________________________



function [loading,dis_itm] = f_load_hs1L_itm_trans(path,plot_cont,loading,dis_itm)
        
        % Initialize
        f_calc = loading.harmonic.f;
        P0     = loading.harmonic.P0;
        nT     = loading.harmonic.nT;
        Nt     = loading.harmonic.Nt;
        
        f_tot  = dis_itm.f_tot;
        Nf_tot = dis_itm.Nf_tot;
             
        %% Frequency spectrum of load
        
        % location of load f_calc in f_tot
        if f_calc == 0; f_calc_tot  = f_calc;
        else;           f_calc_tot  = [-flip(f_calc) f_calc];
        end
        [~,~,loc_f] = intersect(f_calc_tot,f_tot);
        
        % Set vector of amplitude at load frequencies within f_tot
        Pf        = zeros(1,1,Nf_tot);
         if f_calc == 0; P0_tot = P0;
         else;           P0_tot    = [flip(P0) P0];
         end
        Pf(loc_f) = P0_tot;
        
        %% Time/frequency discretization for IFFT -> u(t)
        
        fmin   = min(f_calc);
        T_harm = nT*1/fmin;         % observation time
        
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
        % (documentation Julian - due to fft definition in matlab)
        Pt_harm = 2*Pt_harm;
        
        % Excitation spectrum for plot
        Pf_harm = 1/Nt_harm*fftshift(fft(ifftshift(Pt_harm)));%dt_harm*
        
        %% Output
        loading.Pf      = Pf;      % load vector for multiplication with TF(f) with f=f_tot
%         loading.Pt      = Pt_harm; % with t=t_harm 
        loading.Pt_harm = Pt_harm;      % with t = t_harm
        loading.Pf_harm = Pf_harm;      % with f = f_harm
        
        loading.harmonic.loc_f      = loc_f;        % location of f_calc in f_tot
        loading.harmonic.loc_f_harm = loc_f_harm;   % loaction of f_calc in f_harm
        
        dis_itm.t_harm   = t_harm;
        dis_itm.dt_harm  = dt_harm;
        dis_itm.df_harm  = df_harm;
        dis_itm.T_harm   = T_harm;
        dis_itm.f_harm   = f_harm;
        dis_itm.Nt_harm  = Nt_harm;
        
        %% Plots
        
        %-----------------------------------------------------
        mainfigurename = 'Hs1L - ITM P(t) and P(f)';
        fig = figure('Name',mainfigurename);
        fig.Position = [632,260,1403,1040];
        %-----------------------------------------------------
        sgtitle('Harmonic P(t), P(f)')
        subplot(1,2,1)
        plot(t_harm,Pt_harm)
        xlabel('$t$',    'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel('$P(t)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        subplot(1,2,2)
        stem(f_harm,real(Pf_harm),'Marker','.')
        xlabel('$f$',    'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        ylabel('$P(f)$', 'FontSize', plot_cont.Fontsize,'Interpreter',plot_cont.interpreter)
        
        savename = ['\load_Pz_t_f'];
        saveas(fig,[path.figures savename],'fig')
        %
        if strcmpi(plot_cont.export,'yes')
            exportfig(gcf,[path.figures savename],'Width',plot_cont.exp_width,'Height',plot_cont.exp_heigth,'Fontmode','fixed','Fontsize',plot_cont.Fontsize,'Color','cmyk','resolution',plot_cont.resolution);
        end
        if strcmpi(plot_cont.load,'yes');close(fig);end
        

end

