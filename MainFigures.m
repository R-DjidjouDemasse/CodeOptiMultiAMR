clear
close all
clc

fpath = 'your directory here';

LETTERS = 'A':'Z';

%% Colors
color = [0.4660 0.6740 0.1880; ... % green: S
         0.55   0.18   0.89;   ... % purple: I
         0.00   0.45   0.74;   ... % blue: Re
         1.00   0.00   1.00;   ... % magenta
         0.6350 0.0780 0.1840; ... % brown
         1.00   0.41   0.16;   ... % orange: treatment failure
         0.94   0.74   0.58;   ... % light orange
         0.50   0.50   0.50];      % grey: R


%% ================================================================
% Parameters
%% ================================================================
nbyr=15;
TmaxAB = nbyr*365;

T_Recovery = 4;
T_treatmentStart = 4;

qs = 0.7;

Vect_PropTreated = [0.5; 0.75; 0.95];
l_PropTreated = length(Vect_PropTreated);

Vect_PropTreatB = [0.1; 0.5; 0.95];
l_PropTreatB = length(Vect_PropTreatB);

VectDImmunity = 365;% [180,365];

Vect_RelEfficacy = 0.85;%[0.5,0.85, 0.95];
l_RelEfficacy = length(Vect_RelEfficacy);

% Duration displayed before introduction of B
Tpre = 365;   % one year


%% ================================================================
% Loop over immunity duration
%% ================================================================

for DImmunity = VectDImmunity

    %% ============================================================
    % Loop over relative efficacy of B
    %% ============================================================

    for l1 = 1:l_RelEfficacy

        RelEfficacy = Vect_RelEfficacy(l1);

        figure
        set(gcf,'position',[100,100,1300,2000])

        tiledlayout(1+l_PropTreatB,l_PropTreated,...
                    'TileSpacing','compact',...
                    'Padding','compact');

        FigName = ['FigDynaWithAB_qs_' num2str(qs) ...
                   '_reff_' num2str(RelEfficacy) ...
                   '_DImmunity_' num2str(DImmunity) '.pdf'];

        ColNum = 0;


        %% ========================================================
        % Loop over treatment coverage
        %% ========================================================

        for l2 = 1:l_PropTreated

            RowNum = 0;
            ColNum = ColNum + 1;

            PropTreated = Vect_PropTreated(l2);


            %% ====================================================
            % Loop over proportion receiving treatment B
            %% ====================================================

            for l3 = 1:l_PropTreatB

                RowNum = RowNum + 1;

                PropTreatB = Vect_PropTreatB(l3);


                %% =================================================
                % Model
                %% =================================================

                [R0_AB,Re,PropITF_equiAB,...
                 R0_A,PropITF_equiA,...
                 x,T_AB,...
                 PropCompartment_AB,PropTreatFail_AB,...
                 N0TFA,N0TFB,...
                 EtaStar,xEtaStar,...
                 SR_AB,Infected_AB,...
                 TotPopDynamics_AB,bart_year,tTF_year] = ...
                    MainFunctions(...
                    qs,...
                    T_treatmentStart,...
                    T_Recovery,...
                    PropTreated,...
                    RelEfficacy,...
                    PropTreatB,...
                    TmaxAB,...
                    DImmunity);
                

                %% =================================================
                % Time after introduction of B
                %% =================================================

                TimeYear = T_AB(:)/365;


                %% =================================================
                % AB compartment proportions
                %% =================================================

                PropS  = PropCompartment_AB(:,1);
                PropR  = PropCompartment_AB(:,2);
                PropTS = PropCompartment_AB(:,3);
                PropTF = PropCompartment_AB(:,4);
                PropU  = PropCompartment_AB(:,5);

                PropITF = PropTreatFail_AB(:,1);

                Itf_AB = PropITF_equiAB;


                %% =================================================
                % Pre-period:
                % equilibrium under A only
                %% =================================================

                dt = T_AB(2)-T_AB(1);

                T_Apre = (-Tpre:dt:-dt)';
                TimePre = T_Apre/365;

                % Initial state of AB = equilibrium under A
                PropS_A  = PropS(1);
                PropR_A  = PropR(1);
                PropTS_A = PropTS(1);
                PropTF_A = PropTF(1);
                PropU_A  = PropU(1);

                % Constant equilibrium trajectories
                PropS_pre  = PropS_A  * ones(length(TimePre),1);
                PropR_pre  = PropR_A  * ones(length(TimePre),1);
                PropTS_pre = PropTS_A * ones(length(TimePre),1);
                PropTF_pre = PropTF_A * ones(length(TimePre),1);
                PropU_pre  = PropU_A  * ones(length(TimePre),1);

                % Treatment-failure proportion at equilibrium under A
                PropITF_pre = ...
                    PropITF_equiA * ones(length(TimePre),1);

                % At the positive endemic equilibrium under A:
                % effective reproduction number = 1
                Re_pre = ones(length(TimePre),1);


                %% =================================================
                % Complete trajectories:
                % equilibrium A + AB dynamics
                %% =================================================

                TimePlot = [TimePre; TimeYear];

                PropS_plot  = [PropS_pre;  PropS];
                PropR_plot  = [PropR_pre;  PropR];
                PropTS_plot = [PropTS_pre; PropTS];
                PropTF_plot = [PropTF_pre; PropTF];
                PropU_plot  = [PropU_pre;  PropU];

                PropITF_plot = ...
                    [PropITF_pre; PropITF];

                Re_plot = ...
                    [Re_pre; Re(:)];

                % Total infected population
                PropI_plot = ...
                    PropU_plot + ...
                    PropTF_plot + ...
                    PropTS_plot;


                %% =================================================
                % FIRST ROW:
                % within-host reproduction numbers
                %% =================================================

                if RowNum == 1

                    IdGraph = ColNum;

                    ax = nexttile(IdGraph);
                    ax.FontSize = 14;

                    hold on

                    plot(x,N0TFA,'k','LineWidth',2,'LineStyle','-');

                    plot(x,N0TFB,'r','LineWidth',2,'LineStyle','--');

                    xline(EtaStar,'b','LineWidth',2,'LineStyle','-.');

                    %% eta_A^* label

                    yl = ylim;

                    text(EtaStar,...
                         yl(1)+0.08*(yl(2)-yl(1)),...
                         '$\eta_{\rm A}^*$',...
                         'Interpreter','latex',...
                         'FontSize',14,...
                         'Color',[0 0 1]);


                    %% Relative efficacy indication

                    if ~isempty(xEtaStar)

                        idx = xEtaStar(1);

                        plot([x(idx) x(idx)],...
                             [N0TFB(idx) N0TFA(idx)],...
                             'MarkerFaceColor',[0 0 1],...
                             'MarkerEdgeColor',[0 0 1],...
                             'MarkerSize',3,...
                             'Marker','diamond',...
                             'LineWidth',2,...
                             'Color',[0 0 1]);

                        idxText = min(idx+3,length(x));

                        text(x(idxText),...
                             (N0TFB(idxText)+N0TFA(idxText))/2,...
                             '$r_{\rm eff}$',...
                             'Interpreter','latex',...
                             'FontSize',18,...
                             'Color',[0.49,0.18,0.56]);

                    end


                    xlabel(...
                        'Resistance level $(x)$',...
                        'Interpreter','latex',...
                        'FontSize',14);


                    if IdGraph == 1

                        ylabel(...
                            '\fontsize{9}{0}\selectfont {\rm Within Repro. number}',...
                            'Interpreter','latex');

                    end


                    text(0.5,1.01,...
                        ['\fontsize{20}{0}\selectfont $q^{\rm T}=$' ...
                         num2str(PropTreated,2),...
                         newline,newline,newline ],'Units','normalized',...
                        'Interpreter','latex');
                    
                    
                    ht = title(...
                        ['\fontsize{10}{0}\selectfont\textbf{(' ...
                         LETTERS(IdGraph) ')}' ...
                         '$\quad\mathcal R_{0}^{\rm A}=$' ...
                         num2str(R0_A,2)],...
                        'Interpreter','latex');

                    ht.HorizontalAlignment = 'left';
                    ht.Position(1) = min(x);
                    

                    if IdGraph == 1

                        legend(...
                            '$\mathcal{N}_{\rm A}^{\rm N}$',...
                            '$\mathcal{N}_{\rm B}^{\rm N}$',...
                            '$\eta_{\rm A}^{*}$',...
                            'Interpreter','latex',...
                            'NumColumns',1,...
                            'Location','northeast',...
                            'FontSize',10);

                        legend boxoff

                    end

                    hold off

                end


                %% =================================================
                % BETWEEN-HOST DYNAMICS
                %% =================================================

                IdGraph = RowNum*l_PropTreated + ColNum;

                ax = nexttile(IdGraph);
                ax.FontSize = 14;

                hold on


                %% =================================================
                % LEFT AXIS
                %% =================================================

                yyaxis left


                %% -------------------------------------------------
                % Shaded region:
                % equilibrium under A before introduction of B
                %% -------------------------------------------------

                fill(...
                    [-Tpre/365, 0, 0, -Tpre/365],...
                    [0, 0, 1, 1],...
                    [0.88 0.88 0.88],...
                    'EdgeColor','none',...
                    'FaceAlpha',0.35,...
                    'HandleVisibility','off');

                hold on


                %% -------------------------------------------------
                % Dynamics
                %% -------------------------------------------------

                hS = plot(...
                    TimePlot,...
                    PropS_plot,...
                    'Color',color(1,:),...
                    'LineWidth',2,...
                    'LineStyle','-','Marker','none');


                hI = plot(...
                    TimePlot,...
                    PropI_plot,...
                    'Color',color(2,:),...
                    'LineWidth',2,...
                    'LineStyle','-','Marker','none');


                hR = plot(...
                    TimePlot,...
                    PropR_plot,...
                    'Color',color(8,:),...
                    'LineWidth',2,...
                    'LineStyle','-','Marker','none');


                hITF = plot(...
                    TimePlot,...
                    PropITF_plot,...
                    'Color',color(6,:),...
                    'LineWidth',2,...
                    'LineStyle','-','Marker','none');
                
                xline(bart_year,...
                        'k',...
                        'LineWidth',1.5,...
                        'LineStyle','--',...
                        'HandleVisibility','off');
                    
                xline(tTF_year,...
                    'r',...
                    'LineWidth',1.5,...
                    'LineStyle','-.',...
                    'HandleVisibility','off');


                ylim([0 1])

                xlim([-Tpre/365 TmaxAB/365])


                %% Introduction of B

                xline(0,...
                    'k',...
                    'LineWidth',1.5,...
                    'LineStyle',':',...
                    'HandleVisibility','off');


                %% Left ylabel

                if ismember(IdGraph,[4,7,10])
                    
                    ylabel(...
                        ['\fontsize{14}{0}\selectfont $q_{\rm B}^{\rm T}=$' ...
                         num2str(PropTreatB,2),...
                         newline,newline,...
                         '\fontsize{10}{0}\selectfont {\rm Prop. host pop.}'],...
                        'Interpreter','latex');

                end


                %% =================================================
                % RIGHT AXIS
                %% =================================================

                yyaxis right


                hRe = plot(...
                    TimePlot,...
                    Re_plot,...
                    'Color',color(3,:),...
                    'LineWidth',2,...
                    'LineStyle','--','Marker','none');


                % Re = 1 threshold
                plot(...
                    [TimePlot(1) TimePlot(end)],...
                    [1 1],...
                    'Color',color(3,:),...
                    'LineWidth',1,...
                    'LineStyle',':','Marker','none',...
                    'HandleVisibility','off');


                %% Robust y-limits for Re

                ReMin = min(Re_plot);
                ReMax = max(Re_plot);

                if abs(ReMax-ReMin) < 1e-10

                    ylim([ReMin-0.1 ReMax+0.1])

                else

                    marginRe = 0.05*(ReMax-ReMin);

                    ylim([ReMin-marginRe ...
                          ReMax+marginRe])

                end


                %% Right ylabel

                if ismember(IdGraph,[6,9,12])

                    ylabel(...
                        '\fontsize{10}{0}\selectfont {\rm Effec. repro.} $(\mathcal R_t)$',...
                        'Interpreter','latex');

                end


                %% X axis

                if ismember(IdGraph,10:12)

                    xlabel(...
                        '\fontsize{12}{0}\selectfont {\rm Time since introduction of B (year)}',...
                        'Interpreter','latex');

                end

                xticks(0:1:TmaxAB/365)


                %% Axis colors

                ax.YAxis(1).Color = 'k';
                ax.YAxis(2).Color = color(3,:);


                %% =================================================
                % PANEL TITLE
                %% =================================================


%                 ht = title(...
%                     ['\fontsize{10}{0}\selectfont\textbf{(' ...
%                      LETTERS(IdGraph) ')}' ...
%                      '$\quad\mathcal R_e(0^+)=$' ...
%                      num2str(Re(1),2) ...
%                      '$,\quad {\rm i}_{\rm AB}^{\rm N}=$' ...
%                      num2str(Itf_AB,2)],...
%                     'Interpreter','latex');
                DE=1-bart_year/nbyr; DV=1-tTF_year/nbyr;
                ht = title(...
                    ['\fontsize{10}{0}\selectfont\textbf{(' ...
                     LETTERS(IdGraph) ')}'...
                     ' $\bar t_{\rm E}=$' num2str(bart_year,2) ...
                     ' ($D_E=$' num2str(DE,2) ')' ...
                     '$,\quad \bar t_{\rm V}=$' num2str(tTF_year,2) ...
                     ' (i.e. $D_V=$' num2str(DV,2) ')'],...
                    'Interpreter','latex');
                
                

                ht.HorizontalAlignment = 'left';
                ht.Position(1) = -Tpre/365;



                %% =================================================
                % LEGEND
                %% =================================================

                if IdGraph == 12

                    legend(...
                        [hS,hI,hR,hITF,hRe],...
                        {'$S$',...
                         '$I$',...
                         '$R$',...
                         '${\rm i}^{\rm N}$',...
                         '$\mathcal R_t$'},...
                        'Interpreter','latex',...
                        'NumColumns',1,...
                        'Orientation','horizontal',...
                        'Location','northeast',...
                        'FontSize',10);

                    legend boxoff

                end


                hold off

            end
        end


        %% =========================================================
        % Export current figure
        %% =========================================================

        if ~exist(fpath,'dir')
            mkdir(fpath)
        end

        exportgraphics(...
            gcf,...
            fullfile(fpath,FigName),...
            'ContentType','vector');

    end

end