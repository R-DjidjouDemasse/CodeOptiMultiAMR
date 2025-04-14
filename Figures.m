clear
close all
clc

fpath ='D:\PROJECT IN PROGRESS AND FINISHED';

LETTERS='A':'Z';


%% Effet (taus,qs,qT) sur R0 et IT_Equi (avec A uniquement)
TabR0Avs_taus_qs_qT=0;
if TabR0Avs_taus_qs_qT
    Vect_treatmentStart=[4,7,10]; l_treatmentStart=length(Vect_treatmentStart);
    Vect_PropTreated=0.1:0.05:1; l_PropTreated=length(Vect_PropTreated);
    Vect_qs=0.1:0.05:1; l_qs=length(Vect_qs);
    n_rows=l_treatmentStart*l_PropTreated*l_qs;

    T_Recovery=3; 
    RelEfficacy=0;
    PropTreatB=0;
    TmaxAB=0;
    TmaxA=365*0.6;

    %STRORAGE
    Qs=zeros(n_rows,1);
    Taus=zeros(n_rows,1); 
    Qt=zeros(n_rows,1);
    R0=zeros(n_rows,1);
    ItfEqui=zeros(n_rows,1);

    %BEGIN EXPERIENCE_PLAN
    row_number=0;
    plan_exp=-111*ones(n_rows,3);
    for l0=1:l_treatmentStart
       for l1=1:l_qs
        for l2=1:l_PropTreated 
            row_number=row_number+1;
            plan_exp(row_number,1)=Vect_treatmentStart(l0);
            plan_exp(row_number,2)=Vect_qs(l1);
            plan_exp(row_number,3)=Vect_PropTreated(l2);
         end
       end
    end

    %BEGIN computation on experience plan
    fprintf('Progress:\n');
    fprintf(['\n' repmat('.',1,n_rows ) '\n\n']);

    parfor nl=1:n_rows 
        fprintf('\b|\n');

        param=plan_exp(nl,:);
        T_treatmentStart=param(1);
        qs=param(2);
        PropTreated=param(3);

        [R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
            PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
            EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
            (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB);

        Qs(nl)=qs;
        Taus(nl)=T_treatmentStart; 
        Qt(nl)=PropTreated;
        R0(nl)=R0_A;
        ItfEqui(nl)=PropITF_equiA;
    end

    %nom du tableau des sorties
    tabName='TabR0ItfvsTausQsQt.xlsx';

    %Ici on gènère le contenu du tableau des simus
    Tab=table(Qs,Taus,Qt,R0,ItfEqui);  
    writetable(Tab,fullfile(fpath,tabName))
end

FigTabR0Avs_taus_qs_qT=0;
if FigTabR0Avs_taus_qs_qT
    data=readmatrix([fpath '/TabR0ItfvsTausQsQt.xlsx']);
    %Ordre des colonnes: qs-Taus-qT-R0-ItfEqui

    color= [0 0 0;...
           0.4660 0.6740 0.1880;...  % vert (S)
           0.55 0.18 0.89;...        % violet (U)
           0.00 0.45 0.74;...        %bleu clair(I-TSA)
           1.00 0.50 0.75;...        % rose (I-TFA)
           0.6350 0.0780 0.1840;...  %marron I-TSB
           1.00,0.41,0.16;...  %orange I-TFB
           0.94 0.74 0.58;...  %orange clair (RA)
            0.50,0.50,0.50];        % gris (RB)

    qs=unique(data(:,1));
    taus=unique(data(:,2));
    qT=unique(data(:,3));
    IdGraph=0;
     FigName='FigR0Avs_taus_qs_qT.jpeg';
    figure
    set(gcf,'position',[100,100,1300,300])%3eCoor=lageur, 4e=hauteur
    for IdTaus=1:length(taus)
        IdGraph=IdGraph+1;
        ax=subplot(1,length(taus),IdGraph);
        ax.FontSize = 12;

        nl=length(qs);nc=length(qT);
        R0z=-10*ones(nl,nc);
        for l=1:nl
            for c=1:nc
                row=nl*nc*(IdTaus-1)+(l-1)*nc+c;
                if data(row,2)==taus(IdTaus) && data(row,1)==qs(l) && data(row,3)==qT(c)
                    R0z(l,c)=data(row,4);
                end
            end
        end

        hold on
        j=1;
        plot(qs,R0z(:,1),'color',color(j,:),'linewidth',2,'linestyle','-')
        if IdGraph==1
        annotation('textarrow',[0.2794 0.2554],...
    [0.7589 0.7189],'linewidth',0.05,'linestyle',':');
        text(0.7157,5.019,sprintf('$q^{T}$=0.1'), ...
        'Interpreter', 'latex', 'fontsize', 11);
        annotation('textarrow',[0.1934 0.2174],...
    [0.3478 0.3878],'linewidth',0.05,'linestyle',':');
        text(0.1005,1.4641,sprintf('$q^{T}$=0.95'), ...
        'Interpreter', 'latex', 'fontsize', 11);
        end
        plot(qs,R0z(:,nc-1),'color',color(j,:),'linewidth',2,'linestyle','-.')
        for c=2:2:(nc-2)
            j=j+1;
            plot(qs,R0z(:,c),'color',color(j,:),'linewidth',2)
        end
        yline(1,'k','linewidth',1,'linestyle','--');
        hold off

        title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')}' ...
            '\fontsize{16}{0}\selectfont$\:\tau_s=$' num2str(taus(IdTaus),2)], ...
                           'interpreter','latex');
        xlabel(('\fontsize{15}{0} $q_{\rm C}$'), ...
                                    'interpreter','latex');  
        if IdGraph==1
            ylabel(('\fontsize{12}{0}\selectfont {\rm Reproduction number } ($\mathcal R_0^{\rm A}$)'), ...
                                        'interpreter','latex');  
        end
    end
    saveas(gca, fullfile(fpath,FigName)); 
end

 %% Couleurs
   color= [0.4660 0.6740 0.1880;...  % vert (S)
       0.55 0.18 0.89;...        % violet (U)
       0.00 0.45 0.74;...        %bleu clair(I-TSA)
       1.00 0.00 1.00;...        % rose (I-TFA)
       0.6350 0.0780 0.1840;...  %marron I-TSB
       1.00,0.41,0.16;...  %orange I-TFB
       0.94 0.74 0.58;...  %orange clair (RA)
        0.50,0.50,0.50];        % gris (RB)

%% Figure: Dynamique avec A uniquement
DynaA_uniquement=0;
if DynaA_uniquement
    T_Recovery=3;       % delay before recovery
    T_treatmentStart=4; % delay before treatment
    RelEfficacy=0; % pas utile
    PropTreatB=0;
    TmaxAB=0;
    TmaxA=365;
    Vect_PropTreated=[0.5,0.75,0.95]; l_PropTreated=length(Vect_PropTreated);
    Vect_qs=[0.5,0.75,0.95]; l_Vect_qs=length(Vect_qs);
     FigName='FigDynaWithA.jpeg';
     figure
     set(gcf,'position',[100,100,1300,1300])%3eCoor=lageur, 4e=hauteur
     ColNum=0;

    for l=1:l_PropTreated 
        ColNum=ColNum+1
        RowNum=0;
        PropTreated=Vect_PropTreated(l);
        for k=1:l_Vect_qs
            RowNum=RowNum+1
            qs=Vect_qs(k);


            [R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
                        PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
                        EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
                        (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB);

            PropS=PropCompartment_A(:,1); % S
            PropR=PropCompartment_A(:,2); % R
            PropTF=PropCompartment_A(:,4);
            PropU=PropCompartment_A(:,5);
            PropITF=PropTreatFail_A(:,1);

            % proportions of treatment failure at equilibrium
            Itf_A=PropITF_equiA;

            % Figure basic reproduction number
            if RowNum==1
                 IdGraph=ColNum;
                  ax=subplot(1+l_Vect_qs,l_PropTreated,IdGraph);
                  ax.FontSize = 14;
                  hold on,
                  plot(x,N0TFA,'k','linewidth',2,'linestyle','-');
                  xline(EtaStar,'b','linewidth',2,'linestyle','-.');
                  text(EtaStar,-0.3,'$\eta_A^*$','Interpreter','latex','fontsize',14,'color',[0 0 1])
            %       xticks([min(x) max(x)])
                  plot([x(xEtaStar) x(xEtaStar)], [N0TFB(xEtaStar) N0TFA(xEtaStar)],...
                      'MarkerFaceColor',[0 0 1],...
                      'MarkerEdgeColor',[0 0 1],...
                      'MarkerSize',3,...
                'Marker','diamond','LineWidth',2,'Color',[0 0 1]);

                 xlabel(('\fontsize{10}{0}\selectfont {\rm Resistance level } ($x$)'), ...
                            'interpreter','latex');  

                if IdGraph==1
                    ylabel(('\fontsize{9}{0}\selectfont {\rm Within Repro. number} ($\mathcal N_{\rm A}^{\rm N}$)'), ...
                        'FontWeight','bold',...
                    'interpreter','latex');
                end

                title(['\fontsize{14}{0}\selectfont$q^{\rm T}=$' num2str(PropTreated,2), newline, newline, ...
                   '\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')}'], ...
                   'interpreter','latex');

                hold off
            end


            % Figure dynamics of the between-host model
            time=T_A/365;
            IdGraph=RowNum*l_PropTreated+ColNum;
            ax=subplot(1+l_Vect_qs,l_PropTreated,IdGraph);
            ax.FontSize = 14;
            hold on
            yyaxis left
            plot(time,PropS,'color',color(1,:),'linewidth',2,'linestyle','-');
            plot(time,PropU,'color',color(2,:),'linewidth',2,'linestyle','-');
            plot(time,PropTF,'color',color(6,:),'linewidth',2,'linestyle','-');
            plot(time,PropR,'color',color(8,:),'linewidth',2,'linestyle','-');

            xlim([0 max(time)])
            ylim([0 1])

            if ismember(IdGraph,[4,7,10])
                ylabel(['\fontsize{14}{0}\selectfont $q_{\rm C}=$',num2str(qs,2),newline, newline, ...
                    '\fontsize{10}{0}\selectfont {\rm Prop. host pop.}'], ...
                'interpreter','latex');
            end

            yyaxis right
            plot(time,PropITF,'color',color(3,:),'linewidth',2,'linestyle','--');

            if ismember(IdGraph,[6,9,12])
                ylabel(('\fontsize{10}{0}\selectfont {\rm Prop. failure} $({\rm i}^{\rm N})$'), ...
                        'interpreter','latex');
            end



            if ismember(IdGraph,10:12)
                xlabel(('\fontsize{12}{0}\selectfont {\rm Time (year)}'), ...
                        'interpreter','latex');   
            end


            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = color(3,:);

             title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')} $\:\mathcal R^{\rm A}_0=$' num2str(R0_A,2) ','...
                        '$\:{\rm i}_{\rm A}^{\rm N}=$' num2str(Itf_A,2)], ...
               'interpreter','latex');


            if IdGraph==12
                legend('$S$','$I^{\rm U}$','$I^{\rm N}$','$R$', '${\rm i}^{\rm N}$',...
                         'Interpreter','latex','location','northeast','FontSize', 14,'Orientation','horizontal','NumColumns',2)
                legend boxoff  
            end


            hold off
        end
    end
    saveas(gca, fullfile(fpath,FigName)); 
end

%% Typical Dynamics with A alone
Typic_DynaA_uniquement=0;
if Typic_DynaA_uniquement
    T_Recovery=3;       % delay before recovery
    T_treatmentStart=4; % delay before treatment
    RelEfficacy=0; % pas utile
    PropTreatB=0;
    TmaxAB=0;
    TmaxA=365;
    FigName='FigtypicalDynaWithA.jpeg';
    figure
    set(gcf,'position',[100,150,1734,700])%3eCoor=lageur, 4e=hauteur

    qs=0.5; PropTreated=0.5;
    [R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
                PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
                EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
                (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB);

    PropS=PropCompartment_A(:,1); % S
    PropR=PropCompartment_A(:,2); % R
    PropTS=PropCompartment_A(:,3);
    PropTF=PropCompartment_A(:,4);
    PropU=PropCompartment_A(:,5);
    PropITF=PropTreatFail_A(:,1);

    % proportions of treatment failure at equilibrium
    Itf_A=PropITF_equiA;

    % Figure within-host basic reproduction number
    IdGraph=1;
    ax=subplot(1,3,IdGraph);
    ax.FontSize = 14;
    hold on,
    plot(x,N0TFA,'k','linewidth',3,'linestyle','-');
    xline(EtaStar,'b','linewidth',3,'linestyle','-.');
    xticks([-0.5 0 1 1.5 2 2.5]) %       xticks([min(x) max(x)]) '
%     plot([x(xEtaStar) x(xEtaStar)], [N0TFB(xEtaStar) N0TFA(xEtaStar)],...
%         'MarkerFaceColor',[0 0 1],...
%         'MarkerEdgeColor',[0 0 1],...
%         'MarkerSize',3,...
%         'Marker','diamond','LineWidth',2,'Color',[0 0 1]);
    
    xlabel(('\fontsize{14}{0}\selectfont {\rm Resistance level } ($x$)'), ...
        'interpreter','latex');  
    
    
    ylabel(('\fontsize{15}{0}\selectfont {\rm Within-host Reproduction number} ($\mathcal N_{\rm A}^{\rm N}$)'), ...
        'FontWeight','bold', 'interpreter','latex');
    
    
    title(['\fontsize{14}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')}'], ...
        'interpreter','latex');
    box
    hold off

    % Figure dynamics of the between-host model
    time=T_A/365;
    IdGraph=IdGraph+1;
    ax=subplot(1,3,IdGraph);
    ax.FontSize = 14;
    hold on
    yyaxis left
    plot(time,PropS,'color',color(1,:),'linewidth',2.5,'linestyle','-');
    plot(time,PropU,'color',color(2,:),'linewidth',2.5,'linestyle','-');
    plot(time,PropTS,'color',color(4,:),'linewidth',2.5,'linestyle','-');
    plot(time,PropTF,'color',color(6,:),'linewidth',2.5,'linestyle','-');
    plot(time,PropR,'color',color(8,:),'Marker','none','linewidth',2.5,'linestyle','-');

    xlim([0 max(time)])
    ylim([0 1])
    
    ylabel(('\fontsize{18}{0}\selectfont {\rm Prop. host population}'), ...
                    'interpreter','latex');
    
    yyaxis right
    plot(time,PropITF,'color',color(3,:),'linewidth',2.5,'linestyle','--');
    
    xlabel(('\fontsize{14}{0}\selectfont {\rm Time (year)}'), ...
                        'interpreter','latex');   
    ylabel(('\fontsize{16}{0}\selectfont {\rm Prop. of therapeutic failure} $({\rm i}^{\rm N})$'), ...
                        'interpreter','latex');
    ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = color(3,:);
            
    title(['\fontsize{14}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')} $\:\mathcal R^{\rm A}_0=$' num2str(R0_A,2) ',\quad'...
        '$\:{\rm i}_{\rm A}^{\rm N}=$' num2str(Itf_A,2)], ...
        'interpreter','latex');
           
    legend('$S$','$I^{\rm U}$','$I^{\rm C}$','$I^{\rm N}$','$R$','${\rm i}^{\rm N}$',...
             'Interpreter','latex','location','northeast','FontSize', 18,'Orientation','horizontal','NumColumns',3)
    legend boxoff  
    box
    hold off

     % Figure dynamics of the between-host model
    IdGraph=3;
    ax=subplot(1,3,IdGraph);
    ax.FontSize = 14;
    hold on
     plot(time,EtaTf_tA,'b','linewidth',2.5,'linestyle','-');
     xlabel(('\fontsize{14}{0}\selectfont {\rm Time (year)}'), ...
                        'interpreter','latex'); 
     ylabel(('\fontsize{14}{0}\selectfont {\rm Between-host Aver. resistance level} $(\eta)$'), ...
                    'interpreter','latex');
     ylim([0 0.5])
     title(['\fontsize{14}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')}'], ...
        'interpreter','latex');

     annotation('line',[0.115033967428795 0.114457266160052],...
    [0.78281512605042 0.109957983193277],'Color',[0 0 1],'LineWidth',3,...
    'LineStyle','-.');

     annotation('line',[0.135263580992957 0.137570386067928],[0.87 0.11],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',1.5,...
    'LineStyle','--');

     annotation('textbox',...
    [0.140715109573241 0.29781502955426 0.0224913494809688 0.0554706847314546],...
    'String',{'$\bar{x}_{\rm A}^{*\rm N}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off');

     annotation('textbox',...
    [0.106102709976438 0.0544476518080581 0.0205913410770858 0.0610062882585346],...
    'Color',[0 0 1],...
    'String',{'$\eta_{\rm A}^*$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off');

     annotation('textarrow',[0.955357142857143 0.99632931305716],...
    [0.683823529411765 0.779503105590058],'Color',[0 0 1],'TextLineWidth',1,...
    'TextColor',[0 0 1],...
    'String',{'$\eta_{\rm A}^*=0.41$'},...
    'LineWidth',1,...
    'LineStyle','-.',...
    'Interpreter','latex',...
    'HeadWidth',7,...
    'HeadLength',8,...
    'FontWeight','bold',...
    'FontSize',18);

    box
    hold off
    saveas(gca, fullfile(fpath,FigName)); 
end
%% Figure: Dynamique avec A et B 
DynaAB=1;
if DynaAB
    TmaxAB=2*365+2;
    TmaxA=365*0.6;
    T_Recovery=3;       % delay before recovery
    T_treatmentStart=4; % delay before treatment
    qs=0.7; 
    Vect_RelEfficacy=0.85;           l_RelEfficacy=length(Vect_RelEfficacy);
    FigName=['FigDynaWithAB_qs_' num2str(qs) '_reff_' num2str(Vect_RelEfficacy),'.jpeg'];

    Vect_PropTreated=[0.5;0.75;0.95];    l_PropTreated=length(Vect_PropTreated);
    Vect_PropTreatB=[0.1;0.5;0.95];      l_PropTreatB=length(Vect_PropTreatB);



    figure
    set(gcf,'position',[100,100,1300,2000])%3eCoor=lageur, 4e=hauteur

    for l1=1:l_RelEfficacy
        RelEfficacy=Vect_RelEfficacy(l1); 
        
        ColNum=0;
        for l2=1:l_PropTreated
            RowNum=0;
            ColNum=ColNum+1,
            PropTreated=Vect_PropTreated(l2);
            for l3=1:l_PropTreatB
                RowNum=RowNum+1,
                PropTreatB=Vect_PropTreatB(l3);

                [R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
                    PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
                    EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
                    (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB);

                n1=length(T_AB); 
                % Combinaison of two dynamics
                PropS=[PropCompartment_A(:,1);PropCompartment_AB(:,1)]; % S
                PropR=[PropCompartment_A(:,2);PropCompartment_AB(:,2)]; % R
                PropTF=[PropCompartment_A(:,4);PropCompartment_AB(:,4)];
                PropU=[PropCompartment_A(:,5);PropCompartment_AB(:,5)];
                PropITF=[PropTreatFail_A(:,1);PropTreatFail_AB(:,1)];

                % proportions of treatment failure at equilibrium
                Itf_A=PropITF_equiA;
                Itf_AB=PropITF_equiAB;

                 % Figure basic reproduction number
                if RowNum==1
                   IdGraph=ColNum;
                   ax=subplot(1+l_PropTreatB,l_PropTreated,IdGraph);
                   ax.FontSize = 14;
                   hold on,
                   
                   plot(x,N0TFA,'k','linewidth',2,'linestyle','-');
                   plot(x,N0TFB,'r','linewidth',2,'linestyle','--');
                   xline(EtaStar,'b','linewidth',2,'linestyle','-.');
                   text(EtaStar,-0.3,'$\eta_{\rm A}^*$','Interpreter','latex','fontsize',14,'color',[0 0 1])
                   xticks([min(x) max(x)])
                   plot([x(xEtaStar) x(xEtaStar)], [N0TFB(xEtaStar) N0TFA(xEtaStar)],...
                      'MarkerFaceColor',[0 0 1],...
                      'MarkerEdgeColor',[0 0 1],...
                      'MarkerSize',3,...
                  'Marker','diamond','LineWidth',2,'Color',[0 0 1]);
                  text(x(xEtaStar+3),(N0TFB(xEtaStar+3)+N0TFA(xEtaStar+3))/2,...
                 '${\rm r_{eff}}$','Interpreter','latex','fontsize',18,'color',[0.49,0.18,0.56])

                  xlabel('Resistance level $(x)$','Interpreter','latex','fontsize',14);

                  if IdGraph==1
                      ylabel('\fontsize{9}{0}\selectfont {\rm Within Repro. number}', ...
                            'interpreter','latex');
                  end
                  
                  title(['\fontsize{14}{0}\selectfont$q^{\rm T}=$' num2str(PropTreated,2), newline, newline, ...
                   '\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')}'], ...
                   'interpreter','latex');

                  if IdGraph==1
                     legend('$\mathcal{N}_{\rm A}^{\rm N}$','$\mathcal{N}_{\rm B}^{\rm N}$',...
                                'Interpreter','latex','NumColumns',1,'location','northeast','FontSize', 10)
                     legend boxoff  
                  end

                  hold off
               end

                % Figure dynamics of the between-host model
                IdGraph=RowNum*l_PropTreated+ColNum;
                ax=subplot(1+l_PropTreatB,l_PropTreated,IdGraph);
                ax.FontSize = 14;
                hold on
                yyaxis left

                plot(PropS,'color',color(1,:),'linewidth',2,'linestyle','-');
                plot(PropU,'color',color(2,:),'linewidth',2,'linestyle','-');
                plot(PropTF,'color',color(6,:),'linewidth',2,'linestyle','-');
                plot(PropR,'color',color(8,:),'linewidth',2,'linestyle','-');
                xlim([0 n0+n1])
                ylim([0 1])
        
                if ismember(IdGraph,[4,7,10])
                   ylabel(['\fontsize{14}{0}\selectfont $q_{\rm B}^{\rm T}=$',num2str(PropTreatB,2),newline, newline, ...
                        '\fontsize{10}{0}\selectfont {\rm Prop. host pop.}'], ...
                        'interpreter','latex');
                end

                yyaxis right
                plot(PropITF,'color',color(3,:),'linewidth',2,'linestyle','--');
                fill([0,0,n0,n0], [0,1,1,0], 'c', 'EdgeColor', 'none','FaceAlpha', 0.1);
                ylim([0 0.6])
                xline(n0,'k','linewidth',1,'linestyle',':');
                xticks([0 n0 n0+365 n0+730 n0+n1])
                xticklabels({'','0','1','2',''})

                if ismember(IdGraph,[6,9,12])
                   ylabel(('\fontsize{10}{0}\selectfont {\rm Prop. failure} $({\rm i}^{\rm N})$'), ...
                           'interpreter','latex');
                end

                if ismember(IdGraph,10:12)
                   xlabel(('\fontsize{12}{0}\selectfont {\rm Time (year)}'), ...
                           'interpreter','latex');   
                end
                  
                ax.YAxis(1).Color = 'k';
                ax.YAxis(2).Color = color(3,:);

                title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')} $\:\mathcal R^{\rm AB}_0=$' num2str(R0_AB,2) ','...
                            '$\:{\rm i}_{\rm AB}^{\rm N}=$' num2str(Itf_AB,2)],'interpreter','latex');

                if IdGraph==12
                    legend('$S$','$I^{\rm U}$','$I^{\rm N}$','$R$','${\rm i}^{\rm N}$','$\mbox{Only}\;{\rm A}$',...
                             'Interpreter','latex','NumColumns', 2,'Orientation','horizontal','location','northeast','FontSize', 10)
                    legend boxoff  
                end
              
                hold off

            end
        end
    end

    saveas(gca, fullfile(fpath,FigName));
    
end

%% Figure: le success avec A uniquement demande un qs fort: qs peut être reduit avec l'introduction de B
DynaA_B_discuss=0;
if DynaA_B_discuss
    T_Recovery=3;       % delay before recovery
    T_treatmentStart=4; % delay before treatment
    FigName='FigDynaA_and_B_reduce_qs.jpeg';
    figure
    set(gcf,'position',[100,100,1400,450])%3eCoor=lageur, 4e=hauteur
     
    %modèle avec A uniquqement
    IdGraph=1;
    RelEfficacy=0; % pas utile
    PropTreatB=0;
    TmaxAB=0;
    TmaxA=365;
    PropTreated=0.95;
    
    qs=0.8;
    
    [R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
        PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
        EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
        (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB);
                    
    PropS=PropCompartment_A(:,1); % S
    PropR=PropCompartment_A(:,2); % R
    PropTF=PropCompartment_A(:,4);
    PropU=PropCompartment_A(:,5);
    PropITF=PropTreatFail_A(:,1);
    
    % proportions of treatment failure at equilibrium
    Itf_A=PropITF_equiA;
    
    
    time=T_A/365;
    ax=subplot(1,4,IdGraph);
    ax.FontSize = 14;
    hold on
    yyaxis left
    plot(time,PropS,'color',color(1,:),'linewidth',2,'linestyle','-');
    plot(time,PropU,'color',color(2,:),'linewidth',2,'linestyle','-');
    plot(time,PropTF,'color',color(6,:),'linewidth',2,'linestyle','-');
    plot(time,PropR,'color',color(8,:),'linewidth',2,'linestyle','-');

    xlim([0 max(time)])
    ylim([0 1])
    yticks([0 0.5 1])
    ylabel(('\fontsize{12}{0}\selectfont {\rm Prop. host population}'), ...
                    'interpreter','latex');
    
    yyaxis right
    plot(time,PropITF,'color',color(3,:),'linewidth',2,'linestyle','--');
    
    xlabel(('\fontsize{12}{0}\selectfont {\rm Time (year)}'), ...
                        'interpreter','latex');   
    ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = color(3,:);
            
    title(['\fontsize{12}{0}\selectfont$q_{\rm C}=$' num2str(qs,2), newline, newline, ...
        '\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')} $\:\mathcal R^{\rm A}_0=$' num2str(R0_A,2) ','...
        '$\:q^{\rm T}=$' num2str(PropTreated,2)], ...
        'interpreter','latex');
            
    if IdGraph==1
        legend('$S$','$I^{\rm U}$','$I^{\rm N}$','$R$','${\rm i}^{\rm N}$',...
                 'Interpreter','latex','location','northeast','FontSize', 12,'Orientation','horizontal','NumColumns',2)
        legend boxoff  
    end
    
    hold off
    
    
    qs=0.6;
    IdGraph=IdGraph+1;
    [R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
        PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
        EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
        (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB);
                    
    PropS=PropCompartment_A(:,1); % S
    PropR=PropCompartment_A(:,2); % R
    PropTF=PropCompartment_A(:,4);
    PropU=PropCompartment_A(:,5);
    PropITF=PropTreatFail_A(:,1);
    
    % proportions of treatment failure at equilibrium
    Itf_A=PropITF_equiA;
    
    
    time=T_A/365;
    ax=subplot(1,4,IdGraph);
    ax.FontSize = 14;
    hold on
    yyaxis left
    plot(time,PropS,'color',color(1,:),'linewidth',2,'linestyle','-');
    plot(time,PropU,'color',color(2,:),'linewidth',2,'linestyle','-');
    plot(time,PropTF,'color',color(6,:),'linewidth',2,'linestyle','-');
    plot(time,PropR,'color',color(8,:),'linewidth',2,'linestyle','-');

    xlim([0 max(time)])
    ylim([0 1])
    yticks([0 0.5 1])
    yyaxis right
    plot(time,PropITF,'color',color(3,:),'linewidth',2,'linestyle','--');
    
    xlabel(('\fontsize{12}{0}\selectfont {\rm Time (year)}'), ...
                        'interpreter','latex');   
    ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = color(3,:);
            
    title(['\fontsize{12}{0}\selectfont$q_{\rm C}=$' num2str(qs,2), newline, newline, ...
        '\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')} $\:\mathcal R^{\rm A}_0=$' num2str(R0_A,2) ','...
        '$\:q^{\rm T}=$' num2str(PropTreated,2)], ...
        'interpreter','latex');
    
    hold off
            
   
    %modèle avec A et B 
    RelEfficacy=0.85;
    TmaxAB=4*365+2;
    TmaxA=365*0.6;
    qs=0.6; %on réduit la valeur de qs par rappoer au modèle avec A uniquement
    Vect_PropTreatB=[0.65,0.9]; l_PropTreatB=length(Vect_PropTreatB);
    
    for l=1:l_PropTreatB
        IdGraph=IdGraph+1;
        PropTreatB=Vect_PropTreatB(l);
        
        [R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
            PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
            EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
            (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB);
        
        n1=length(T_AB); 
        % Combinaison of two dynamics
        PropS=[PropCompartment_A(:,1);PropCompartment_AB(:,1)]; % S
        PropR=[PropCompartment_A(:,2);PropCompartment_AB(:,2)]; % R
        PropTF=[PropCompartment_A(:,4);PropCompartment_AB(:,4)];
        PropU=[PropCompartment_A(:,5);PropCompartment_AB(:,5)];
        PropITF=[PropTreatFail_A(:,1);PropTreatFail_AB(:,1)];

        % proportions of treatment failure at equilibrium
        Itf_A=PropITF_equiA;
        Itf_AB=PropITF_equiAB;
        
        ax=subplot(1,4,IdGraph);
        ax.FontSize = 14;
        hold on
        yyaxis left
        plot(PropS,'color',color(1,:),'linewidth',2,'linestyle','-');
        plot(PropU,'color',color(2,:),'linewidth',2,'linestyle','-');
        plot(PropTF,'color',color(6,:),'linewidth',2,'linestyle','-');
        plot(PropR,'color',color(8,:),'linewidth',2,'linestyle','-');
        xlim([0 n0+n1])
        ylim([0 1])
        yticks([0 0.5 1])
        yyaxis right
        plot(PropITF,'color',color(3,:),'linewidth',2,'linestyle','--');
        xline(n0,'k','linewidth',1,'linestyle',':');
        xticks([0 n0 n0+365 n0+730 n0+1095 n0+1460 n0+n1])
        xticklabels({'','0','1','2','3', '4',''})
        fill([0,0,n0,n0], [0,1,1,0], 'c', 'EdgeColor', 'none','FaceAlpha', 0.1);
        if IdGraph==2
            ylabel(('\fontsize{10}{0}\selectfont {\rm Prop. of therapeutic failure} $({\rm i}^{\rm N})$'), ...
                    'interpreter','latex');
        end
        if IdGraph==4
            ylabel(('\fontsize{10}{0}\selectfont {\rm Prop. of therapeutic failure} $({\rm i}^{\rm N})$'), ...
                    'interpreter','latex');
        end
        
        xlabel('Time (year)','Interpreter','latex','fontsize',20) 
         
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = color(3,:);
        
        title(['\fontsize{12}{0}\selectfont$q_{\rm C}=$' num2str(qs,2), newline, newline, ...
           '\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(IdGraph) ')} $\:\mathcal R^{\rm AB}_0=$' num2str(R0_AB,2) ','...
                    '$\:q^{\rm T}_{\rm B}=$' num2str(PropTreatB,2)], ...
           'interpreter','latex');
       
        hold off
        
    end
    
    saveas(gca, fullfile(fpath,FigName)); 
end

%% Pour (reff,taus,qs,qT) donné, on cherche le plus petit qTB (s'il existe) à partir duquel on a R0AB<1
Tab_Min_qTB=0;
if Tab_Min_qTB
    Vect_reff=[0.5,0.7,0.85]; l_reff=length(Vect_reff);
    Vect_taus=[4,7,10]; l_taus=length(Vect_taus);
    Vect_qs=0.05:0.05:0.95; l_qs=length(Vect_qs);
    Vect_qt=0.05:0.05:0.95; l_qt=length(Vect_qt);
    n_rows=l_reff*l_taus*l_qs*l_qt;

    T_Recovery=3; 
    TmaxAB=1;
    TmaxA=365*0.6;

    %STRORAGE
    reff=zeros(n_rows,1);
    Taus=zeros(n_rows,1); 
    Qs=zeros(n_rows,1);
    Qt=zeros(n_rows,1);
    R0A=zeros(n_rows,1);
    R0AB=zeros(n_rows,1);
    qTB_opti=ones(n_rows,1);

    %BEGIN EXPERIENCE_PLAN
    row_number=0;
    plan_exp=-111*ones(n_rows,3);
    for l0=1:l_reff
       for l1=1:l_taus
        for l2=1:l_qs
            for l3=1:l_qt
                row_number=row_number+1;
                plan_exp(row_number,1)=Vect_reff(l0);
                plan_exp(row_number,2)=Vect_taus(l1);
                plan_exp(row_number,3)=Vect_qs(l2);
                plan_exp(row_number,4)=Vect_qt(l3);
            end
         end
       end
    end

    %BEGIN computation on experience plan
    fprintf('Progress:\n');
    fprintf(['\n' repmat('.',1,n_rows ) '\n\n']);

    parfor nl=1:n_rows 
        fprintf('\b|\n');

        param=plan_exp(nl,:);
        RelEfficacy=param(1);
        T_treatmentStart=param(2);
        qs=param(3);
        PropTreated=param(4);

        [R0_A,R0_AB,Min_qTB]=Main_Min_qTB...
            (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,TmaxA,TmaxAB)

        reff(nl)=RelEfficacy;
        Taus(nl)=T_treatmentStart; 
        Qs(nl)=qs;
        Qt(nl)=PropTreated;
        R0A(nl)=R0_A;
        R0AB(nl)=R0_AB;
        qTB_opti(nl)=Min_qTB;
    end

    %nom du tableau des sorties
    tabName='Tab_Min_qTB.xlsx';

    %Ici on gènère le contenu du tableau des simus
    Tab=table(reff,Taus,Qs,Qt,R0A,R0AB,qTB_opti);  
    writetable(Tab,fullfile(fpath,tabName))
end


%% QUELQUES fonctions

function[R0_A,R0_AB,Min_qTB]=Main_Min_qTB(qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,TmaxA,TmaxAB)
    qTB=0:0.05:1;
    cp=1;Id=0;Min_qTB=2;
    while cp
       Id=Id+1;
       if Id>length(qTB)
           cp=0;
       else
           PropTreatB=qTB(Id);
           [R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
               PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
               EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
               (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB);
           if R0_AB<1
               Min_qTB=PropTreatB;
               cp=0;
           end
       end
    end     
end
    
