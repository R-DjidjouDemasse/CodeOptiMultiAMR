function [within_out]=Test00_WithinHostQUASAROPTI() 
clear all
close all
clc

% Fixed parameters
T_treatmentStart=4;  % Percentage increase of the initial bacterial load
T_Recovery=3;
    %% Parameters coming from the within-host scale
   VectDrugGroup={'A','B'}; NbTreat=2; 
    bTF=cell(1,NbTreat); bTS=cell(1,NbTreat);
    totpopbTF=cell(1,NbTreat); totpopbTS=cell(1,NbTreat);
   
    betaTF=cell(1,NbTreat); betaTS=cell(1,NbTreat);
    alphaTF=cell(1,NbTreat); alphaTS=cell(1,NbTreat);
    
    gammaTF=cell(1,NbTreat); gammaTS=cell(1,NbTreat);
    omegaTF_U=cell(1,NbTreat); omegaTS_U=cell(1,NbTreat); 
    
    PgammaTF=cell(1,NbTreat); PgammaTS=cell(1,NbTreat);
    PomegaTF_U=cell(1,NbTreat); PomegaTS_U=cell(1,NbTreat);
    
    xBarTF=cell(1,NbTreat); xBarTS=cell(1,NbTreat);
    N0TF=cell(1,NbTreat); N0TS=cell(1,NbTreat);
   
    for l=1:NbTreat
        Drug=VectDrugGroup{l}; EtaStar=0; N0AEtaStar=1; N0BEtaStar=1;
       [Tau,dtau,Ntau,x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,PomegaTU,...
           PomegaUT,Pgamma,beta,alpha,xBar,ds,IndexT_treatmentStart,IndexT_Recovery]=...
        Withinhost('U',Drug,T_treatmentStart,T_Recovery,1,0,0);
        omegaU_T=omegaUT;  PomegaU_T=PomegaUT; PgammaU=Pgamma;
        betaU=beta; alphaU=alpha; xBarU=xBar; N0U=N0; 
        bU=b; totpopbU=totpopb; dsU=ds;
        gammaU=gamma
     
     [Tau,dtau,Ntau,x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,PomegaTU,...
         PomegaUT,Pgamma,beta,alpha,xBar,ds,IndexT_treatmentStart,IndexT_Recovery]=...
     Withinhost('TF',Drug,T_treatmentStart,T_Recovery,EtaStar,N0AEtaStar,N0BEtaStar);
     omegaTF_U{l}=omegaTU;  PomegaTF_U{l}=PomegaTU; PgammaTF{l}=Pgamma; 
     betaTF{l}=beta; alphaTF{l}=alpha;  xBarTF{l}=xBar; N0TF{l}=N0; 
     bTF{l}=b; totpopbTF{l}=totpopb;
     gammaTF{l}=gamma
       
     [Tau,dtau,Ntau,x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,PomegaTU,...
         PomegaUT,Pgamma,beta,alpha,xBar,ds,IndexT_treatmentStart,IndexT_Recovery]=...
     Withinhost('TS',Drug,T_treatmentStart,T_Recovery,EtaStar,N0AEtaStar,N0BEtaStar);
       omegaTS_U{l}=omegaTU; gammaTS{l}=gamma; PomegaTS_U{l}=PomegaTU; PgammaTS{l}=Pgamma; 
       betaTS{l}=beta; alphaTS{l}=alpha;  xBarTS{l}=xBar; N0TS{l}=N0;
       bTS{l}=b; totpopbTS{l}=totpopb;
    end

%% N0^U, N0^N & N0^C : Within-host basic reproduction number
    Figure2=0;
    if Figure2
        %Figure 2A
        figure()
        N0UA = tiledlayout(1,1);
        nexttile
        hold on,
        plot(x,N0U,'r','linewidth',2,'linestyle','-');
        xlabel('Level of resistance $(x)$','Interpreter','latex','fontsize',15);
        ylabel('Basic reproduction number $(\mathcal{N}^{\rm U})$','Interpreter','latex','fontsize',14)
        box
        hold off
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig2A','.jpeg'];
        N0UA.Units = 'centimeters';
        exportgraphics(N0UA,FigName)

        %Figure 2D
        figure()
        N0TFA = tiledlayout(1,1);
        nexttile
        hold on,
        plot(x,N0TF{1},'k','linewidth',2,'linestyle','-');
        annotation('line',[0.550595238095238 0.550595238095238],...
            [0.873809523809524 0.137301587301587],...
            'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
            'LineWidth',1.5,...
            'LineStyle','--');
        annotation('textbox',...
            [0.550595238095236 0.318453619396879 0.0809523809523809 0.0924511425078867],...
            'String',{'$\bar{x}^{*{\rm N}}$'},...
            'LineStyle','none',...
            'Interpreter','latex',...
            'FontWeight','bold',...
            'FontSize',16,...
            'FitBoxToText','off');
        xlabel('Level of resistance $(x)$','Interpreter','latex','fontsize',15);
        ylabel('Basic reproduction number $(\mathcal{N}^{\rm N})$','Interpreter','latex','fontsize',14) 
        box
        hold off
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig2D','.jpeg'];
        N0TFA.Units = 'centimeters';
        exportgraphics(N0TFA,FigName)

        % Figure 2G
        figure()
        N0TSA = tiledlayout(1,1);
        nexttile
        hold on,
        plot(x,N0TS{1},'b','linewidth',2,'linestyle','-');
        annotation('line',[0.589285714285714 0.588095238095238],...
            [0.887301587301587 0.136507936507937],'Color',[0 0 1],'LineWidth',1.5,...
            'LineStyle','--');
        annotation('textbox',...
            [0.591666666666666 0.306937665969608 0.0845238095238103 0.0952369372049971],...
            'Color',[0 0 1],...
            'String',{'$\bar{x}^{*{\rm C}}$'},...
            'LineStyle','none',...
            'Interpreter','latex',...
            'FontSize',16,...
            'FitBoxToText','off');
        xlabel('Level of resistance $(x)$','Interpreter','latex','fontsize',15);
        ylabel('Basic reproduction number $(\mathcal{N}^{\rm C})$','Interpreter','latex','fontsize',14)  
        box
        hold off
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig1A','.jpeg'];
        N0TSA.Units = 'centimeters';
        exportgraphics(N0TSA,FigName)

        %% bU(tau,x), bN(tau,x) & bC(tau,x)  
        %Figure 2B
        figure()
        bu = tiledlayout(1,1);
        nexttile
        h=surf(x(1:1:Nx),Tau(1:1:Ntau),bU(1:1:Ntau,1:1:Nx),'FaceColor','interp',...
           'EdgeColor','none',...
           'FaceLighting','gouraud');
        axis tight
        view(-37,83)
        camlight left
        set(h,'edgecolor','none')
        xlim([-1 1])
        xlabel('Level of resistance (x)','fontsize',13,'position',[-0.011754178861501,-4.781539005774364,-0.008358898580141],'rotation',28);
        ylabel('Time since infection (days)','fontsize',13,'position',[-1.302961670952575,15.015508814912636,0.024062364396574],'rotation',-43);
        zlabel('b^U(\tau,x)','fontsize',14,'FontWeight','bold');
        title('','Interpreter','latex','FontSize', 25);
        annotation('textarrow',[0.384523809523809 0.354166666666667],...
            [0.82875605815832 0.746825396825397],'String','$\bar{x}^{*{\rm U}}$',...
            'LineWidth',0.6,...
            'LineStyle','--',...
            'Interpreter','latex',...
            'HeadWidth',8,...
            'HeadStyle','cback3',...
            'HeadLength',7,...
            'FontSize',14,...
            'FontName','Arial');
        FigName=['fig2B','.png'];
        bu.Units = 'centimeters';
        exportgraphics(bu,FigName) 

        % Figure 2E
        figure()
        btf = tiledlayout(1,1);
        nexttile
        h=surf(x(1:1:Nx),Tau(1:1:Ntau),bTF{1}(1:1:Ntau,1:1:Nx),'FaceColor','interp',...
           'EdgeColor','none',...
           'FaceLighting','gouraud');
        axis tight
        view(-37,83)
        camlight left
        set(h,'edgecolor','none')
        xlabel('Level of resistance (x)','fontsize',13,'position',[0.606588140593461,-4.205845873736592,0.055221325802734],'rotation',27);
        ylabel('Time since infection (days)','fontsize',13,'position',[-1.443431457162228,15.033076820497499,0.08490065183554],'rotation',-43);
        zlabel('b^{N}(\tau,x)','fontsize',14,'FontWeight','bold');
        title('','Interpreter','latex','FontSize', 25);
        annotation('textarrow',[0.399404761904762 0.376785714285714],...
            [0.824603174603175 0.737301587301587],'String','$\bar{x}^{*{\rm N}}$',...
            'LineWidth',0.6,...
            'LineStyle','--',...
            'Interpreter','latex',...
            'HeadWidth',8,...
            'HeadStyle','cback3',...
            'HeadLength',7,...
            'FontSize',14,...
            'FontName','Arial');
        FigName=['fig2E','.png'];
        btf.Units = 'centimeters';
        exportgraphics(btf,FigName)

        % Figure 2H
        figure()
        bts = tiledlayout(1,1);
        nexttile
        h=surf(x(1:1:Nx),Tau(1:1:Ntau),bTS{1}(1:1:Ntau,1:1:Nx),'FaceColor','interp',...
           'EdgeColor','none',...
           'FaceLighting','gouraud');
        axis tight
        view(-37,83)
        camlight left
        set(h,'edgecolor','none')
        xlim([-1 1])
        xlabel('Level of resistance (x)','fontsize',13,'position',[-0.011754178861501,-4.781539005774364,-0.008358898580141],'rotation',28);
        ylabel('Time since infection (days)','fontsize',13,'position',[-1.302961670952575,15.015508814912636,0.024062364396574],'rotation',-43);
        zlabel('b^{C}(\tau,x)','fontsize',14,'FontWeight','bold');
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig2H','.png'];
        bts.Units = 'centimeters';
        exportgraphics(bts,FigName)
        %% B^U, B^N & B^C Total bacterial load
        %Figure 2C
        figure()
        BU = tiledlayout(1,1);
        nexttile
        hold on,
        plot(Tau,(totpopbU),'r','linewidth',2,'linestyle','-');
        xlabel('Time since infection (days)','Interpreter','latex','fontsize',15);
        ylabel('$B^{\rm U}(\tau)$','Interpreter','latex','fontsize',15);
        box
        hold off
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig2C','.png'];
        BU.Units = 'centimeters';
        exportgraphics(BU,FigName)

        % Figure 2F
        figure()
        BTF = tiledlayout(1,1);
        nexttile
        hold on,
        plot(Tau,(totpopbTF{1}),'k','linewidth',2,'linestyle','-');
        yline((totpopbTF{1}(1)),'k','linewidth',1.5,'linestyle',':');
        xlabel('Time since infection (days)','Interpreter','latex','fontsize',15);
        ylabel('$B^{\rm N}(\tau)$','Interpreter','latex','fontsize',15);
        box
        hold off
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig2F','.png'];
        BTF.Units = 'centimeters';
        exportgraphics(BTF,FigName)

        % Figure 2I
        figure()
        BTS = tiledlayout(1,1);
        nexttile
        hold on,
        plot(Tau,(totpopbTS{1}),'b','linewidth',2,'linestyle','-');
        xlabel('Time since infection (days)','Interpreter','latex','fontsize',15);
        ylabel('$B^{\rm C}(\tau)$','Interpreter','latex','fontsize',15);
        box
        hold off
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig2I','.jpeg'];
        BTS.Units = 'centimeters';
        exportgraphics(BTS,FigName)
    end
%% Figure 3 
    Figure3=0;
    if Figure3
        % Figure 3A Prob. of Remaining Untreated
        figure()
        PU_T = tiledlayout(1,1);
        nexttile
        hold on,
        plot(Tau,PomegaU_T,'b','linewidth',2,'linestyle','-');
        xlabel('Time since infection (days)','Interpreter','latex','fontsize',15);
        ylabel(['Prob. of Remaining Untreated                             ';'$\Big(e^{-\int_0^{\tau}\omega^{\rm U}_{\rm T}(s)ds}\Big)$'],...
            'FontWeight','bold',...
            'FontSize',16,...
            'Interpreter','latex');
        xticks([0 6 10 15 20 25 30])
        annotation('textbox',...
            [0.255952380952381 0.0561428571428574 0.0803571428571427 0.1],...
            'String',{'\tau_s'},...
            'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',16,...
            'FontName','Serif',...
            'FitBoxToText','off');
        box
        hold off
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig3A','.jpeg'];
        PU_T.Units = 'centimeters';
        exportgraphics(PU_T,FigName)

        %Figure 3B Prob. of Remaining Non-compliant
        figure()
        PTF_U = tiledlayout(1,1);
        nexttile
        hold on,
        plot(Tau,PomegaTF_U{1},'b','linewidth',2,'linestyle','-');
        xlabel('Time since infection (days)','Interpreter','latex','fontsize',15);
        ylabel(['Prob. of Remaining Non-compliant                         ';'$\Big(e^{-\int_0^{\tau}\omega^{\rm N}_{\rm U}(s)ds}\Big)$'],...
            'FontWeight','bold',...
            'FontSize',14,...
            'Interpreter','latex');
        box
        hold off
        FigName=['fig3B','.jpeg'];
        PTF_U.Units = 'centimeters';
        exportgraphics(PTF_U,FigName)

        % Figure 3C Prob. of Remaining Compliant
        figure()
        PgamT = tiledlayout(1,1);
        nexttile
        hold on,
        plot(Tau,PgammaTS{1},'b','linewidth',2,'linestyle','-');
        xlabel('Time since infection (days)','Interpreter','latex','fontsize',15);
        ylabel(['Prob. of Remaining Compliant                     ';'$\Big(e^{-\int_0^{\tau}\gamma^{\rm C}(s)ds}\Big)$'],...
            'FontWeight','bold',...
            'FontSize',16,...
            'Interpreter','latex');
        annotation('textbox',...
            [0.233333333333333 0.0696349206349208 0.0714285714285713 0.1],...
            'String','\tau_r',...
            'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',16,...
            'FontName','Serif',...
            'FitBoxToText','off');
        box
        hold off
        title('','Interpreter','latex','FontSize', 25);
        FigName=['fig3C','.jpeg'];
        PgamT.Units = 'centimeters';
        exportgraphics(PgamT,FigName)
    end
end
%% Fonction within-host Model
function[Tau,dtau,Ntau,x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,PomegaTU,...
    PomegaUT,Pgamma,beta,alpha,xBar,ds,IndexT_treatmentStart,IndexT_Recovery]=...
     Withinhost(TreatStatus,Drug,T_treatmentStart,T_Recovery,EtaStar,N0AEtaStar,N0BEtaStar)
% Discretization
dx=0.01; % discretization step
%xmin=-0.5; xmax=1.5; % Resistance space 
xmin=-1; xmax=2; % Resistance space
x=xmin:dx:xmax; % discretization of resistance level x
Nx=length(x); % vector length of x

dtau=0.1/20; 
Taumax=30; % Maximal time of infection age
Tau=0:dtau:Taumax;
Ntau=length(Tau); % vector length of infection age

% Fix parameters 
pm=10; p0=0.95*pm; kappa=1; m0=0.05; sigma0=0.05; p1_0=0.5;
p1=p0*p1_0;
r1=9e3; beta0=1.848e-3/1.5; alpha0=7.5e-2; % For alpha and beta

% immune effect
  N00=21;
  mu=1.81*p0/N00;
  
%functional parameters p
  pNum=1+(pm/p0-1)*(p0*(pm-p1)/(p1*(pm-p0))).^x;
  pp=(pm./pNum);
 
%functional parameter k
if strcmp(TreatStatus,'U')
     eps=0.05/5;
     kk=zeros(1,Nx); 
elseif strcmp(TreatStatus,'TS')
    eps=0.05;
    if strcmp(Drug,'A')
        k0=15; k1_0=0.3;  k1=k0*k1_0; % TS
        kk=(k0*(k1/k0).^(x));
    elseif strcmp(Drug,'B')
        k0=15; k1_0=0.3;  k1=k0*k1_0;
        kk=((mu*(1-(N0BEtaStar/N0AEtaStar))+(k0*(k1/k0)^EtaStar))/(N0BEtaStar/N0AEtaStar))*...
            (k1/k0).^(x-EtaStar);
    end
elseif strcmp(TreatStatus,'TF')
    eps=0.05;
    if strcmp(Drug,'A')
        k0=3; k1_0=0.01;  k1=k0*k1_0;  
        kk=(k0*(k1/k0).^(x));
    elseif strcmp(Drug,'B')
        k0=3; k1_0=0.01;  k1=k0*k1_0;  
        kk=((mu*(1-(N0BEtaStar/N0AEtaStar))+(k0*(k1/k0)^EtaStar))/(N0BEtaStar/N0AEtaStar))*...
            (k1/k0).^(x-EtaStar);
    end
end

    % Mutation kernel
         J0=normpdf(dx*(1-Nx:Nx-1),0,eps)';
         
  % Within-host basic reproduction
         N0=pp./(mu+kk);

% Initial conditions
         b=zeros(Ntau,Nx); % Initialization of bacterial state
         totpopb=zeros(1,Ntau); % Initialization of total bacterial state
         b(1,:)=m0*normpdf(x,0,sigma0); % initial condition of b
         bb=b(1,:); totpopb(1)=sum(bb); 
         
        PomegaTU=ones(1,Ntau);
        PomegaUT=ones(1,Ntau);
        Pgamma=ones(1,Ntau);
%% Solve of the within-host model
   for tau=1:Ntau-1
      B=(1+totpopb(tau))^(-kappa);
      LT=B*dx*conv(pp.*b(tau,:),J0,'same');    
      b(tau+1,:)=(b(tau,:)+dtau*LT)./(1+dtau*mu+dtau*kk); 
      bb=b(tau+1,:);
      totpopb(tau+1)=sum(bb);
   end 
   
      IndexT_Recovery=find(Tau>0 & abs(Tau-T_Recovery)<0.00001);
      IndexT_treatmentStart=find(Tau>0 & abs(Tau-T_treatmentStart)<0.00001);
      ds=totpopb(IndexT_treatmentStart)/totpopb(1);
      SeuilGuerison=totpopb(IndexT_Recovery)/totpopb(1);
   for tau=1:Ntau-1
      Id=1:(tau+1);
      kz=100;
%       SeuilGuerison=10^(-3);
%       ds=2;
      PomegaTU(tau+1)= exp(-trapz(Tau(Id),  kz*(totpopb(Id)>totpopb(1))   )    );
      PomegaUT(tau+1)= exp( -trapz( Tau(Id), kz*(totpopb(Id)>ds*totpopb(1))    )   ); 
      Pgamma(tau+1)= exp(-trapz(Tau(Id),kz*(totpopb(Id)<=SeuilGuerison*totpopb(1))));
   end
   
   % Individual average level of resistance
    xBar=zeros(1,Ntau);
  for n=1:Ntau
    xBar(n)=sum(x.*b(n,:))/totpopb(n);
  end
  
    % Transition rates of between-host model
    omegaTU=kz*(totpopb>totpopb(1));
    omegaUT=kz*(totpopb>ds*totpopb(1));
    gamma=kz*(totpopb<=SeuilGuerison*totpopb(1));
    
    % Transmission and disease-induced mortality rates
    beta=(beta0*totpopb)./(r1+totpopb);
    alpha=(alpha0*totpopb)./(r1+totpopb);
end