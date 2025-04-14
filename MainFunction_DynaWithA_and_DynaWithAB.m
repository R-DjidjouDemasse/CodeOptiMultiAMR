function[R0_AB,PropITF_equiAB,R0_A,PropITF_equiA,x,n0,T_A,T_AB,PropCompartment_A,...
    PropTreatFail_A,PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
    EtaStar,xEtaStar,EtaTf_tA]=MainFunction_DynaWithA_and_DynaWithAB...
    (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxA,TmaxAB)
%% Fixed parameters
    muh=5.2675e-2;      % death rate of humans
    Lambdah=5e4;        % recruitment rate
    prev=0.2;           % Initial prevalence
    DFE=Lambdah/muh;

%% Dynamics with only one treatment qTA=1

    % Age vector to define initial conditions for infectious states
    dtau=0.1/20; 
    Taumax=30; % Maximal time of infection age
    Tau=0:dtau:Taumax;
    Ntau=length(Tau); % vector length of infection age

    % Initial conditions for S, R_A, R_B, ITFA, ITSA, ITFA, ITFB, IU
    InitSR_UniqA=[DFE*(1-prev) 0 0]; 

    gamma0=10*log(10); 
    InitIU=(0.99*prev*DFE*gamma0)*exp(-gamma0.*Tau); 
    InitTSA=(0.009*prev*DFE*gamma0)*exp(-gamma0.*Tau); 
    InitTFA=zeros(1,Ntau);
    InitTFB=zeros(1,Ntau);
    InitTSB=(0.001*prev*DFE*gamma0)*exp(-gamma0.*Tau); 

    InitInfected_UniqA=[InitTFA InitTSA InitTFB InitTSB InitIU];

    % Model run with only one treatment
    %TmaxA=200;
    [x,T_A,N0_UniqA,R0_A,SR_UniqA,Infected_UniqA,PropCompartment_A,PropTreatFail_A,EtaTf_UniqA,PropITF_equiA,EtaTf_tA,Eta_tA]...
    =MainBetweenhost(TmaxA,Tau,dtau,Ntau,PropTreated,qs,0,...
     T_treatmentStart,T_Recovery,muh,Lambdah,InitSR_UniqA,InitInfected_UniqA,1,1,0);

    n0=length(T_A); nx=length(x); 
    N0TFA=N0_UniqA(1:nx); 

%% Dynamics with two treatments AB
    % Initial conditions 
    InitSR_AB=[SR_UniqA(n0,1) SR_UniqA(n0,2) SR_UniqA(n0,3)]; %S RA RB

    InitTFA=Infected_UniqA(:,n0);   %ITFA 
    InitTSA=Infected_UniqA(:,2*n0); %ITSA
    InitTFB=Infected_UniqA(:,3*n0); %ITFB
    InitTSB=Infected_UniqA(:,4*n0); %ITSB
    InitIU=Infected_UniqA(:,5*n0);  %IU
    InitInfected_AB=[InitTFA InitTSA InitTFB InitTSB InitIU];

    EtaStar=EtaTf_UniqA;
    xEtaStar=find(abs(x-round(EtaStar,2))<0.0001);
    N0AEtaStar=N0TFA(xEtaStar);
    N0BEtaStar=(1-RelEfficacy)*N0AEtaStar;

    %TmaxAB=2000;
    [x,T_AB,N0_AB,R0_AB,SR_AB,Infected_AB,PropCompartment_AB,PropTreatFail_AB,EtaTf_AB,PropITF_equiAB,EtaTf_tAB,Eta_tAB]...
    =MainBetweenhost(TmaxAB,Tau,dtau,Ntau,PropTreated,qs,PropTreatB,...
    T_treatmentStart,T_Recovery,muh,Lambdah,InitSR_AB,InitInfected_AB,N0BEtaStar,N0AEtaStar,EtaStar);
    
    N0TFB=N0_AB(2*nx+1:3*nx);
    
    
    
    %% Function Between-host model
function[x,T,N0_WH,R0,SR,Infected,PropCompartment,PropTreatFail,EtaTf,PropITF_equi,EtaTf_t,Eta_t]...
          =MainBetweenhost(Tmax,Tau,dtau,Ntau,PropTreated,qs,PropTreatB,...
             T_treatmentStart,T_Recovery,muh,Lambdah,InitSR,InitInfected,N0BEtaStar,N0AEtaStar,EtaStar)
%% Fixed Parameters of the between-host model
    NbTreat=2; % Number of treatment
    QT=PropTreated.*[1-PropTreatB PropTreatB];
    qT=cell(1,NbTreat);
    qT{1}=QT(1); qT{2}=QT(2); 
   
%% Parameters coming from the within-host scale
    VectDrugGroup={'A','B'};
    betaTF=cell(1,NbTreat); betaTS=cell(1,NbTreat);
    alphaTF=cell(1,NbTreat); alphaTS=cell(1,NbTreat);
    gammaTS=cell(1,NbTreat);
    omegaTF_U=cell(1,NbTreat); omegaTS_U=cell(1,NbTreat); 
    xBarTF=cell(1,NbTreat); xBarTS=cell(1,NbTreat);
    N0TF=cell(1,NbTreat); N0TS=cell(1,NbTreat);
   
    for l=1:NbTreat
        Drug=VectDrugGroup{l};
        [x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,beta,alpha,xBar]=...
        MainWithinhost(Tau,dtau,Ntau,'U',Drug,T_treatmentStart,T_Recovery,N0BEtaStar,N0AEtaStar);
        omegaU_T=omegaUT; betaU=beta; alphaU=alpha; xBarU=xBar; N0U=N0;

        [x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,beta,alpha,xBar]=...
        MainWithinhost(Tau,dtau,Ntau,'TF',Drug,T_treatmentStart,T_Recovery,N0BEtaStar,N0AEtaStar,EtaStar);

        omegaTF_U{l}=omegaTU; betaTF{l}=beta; 
        alphaTF{l}=alpha;  xBarTF{l}=xBar; N0TF{l}=N0;

        [x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,beta,alpha,xBar]=...
        MainWithinhost(Tau,dtau,Ntau,'TS',Drug,T_treatmentStart,T_Recovery,N0BEtaStar,N0AEtaStar,EtaStar);
        omegaTS_U{l}=omegaTU; gammaTS{l}=gamma; betaTS{l}=beta; 
        alphaTS{l}=alpha;  xBarTS{l}=xBar; N0TS{l}=N0;
    end
    % Vector of N0 within-host {1}=A & {2}=B
    N0_WH=[N0TF{1} N0TS{1} N0TF{2} N0TS{2} N0U];
    
%% Basic Reproduction Number R0
    Phi=zeros(5,5,Ntau); Chi=zeros(1,Ntau); 
    QQ=[qs*qT{1};qs*qT{2};(1-qs)*qT{1};(1-qs)*qT{2};1-qT{1}-qT{2}];
    for ktau=1:Ntau
        Phi11=-((1-qT{1})*omegaTS_U{1}(ktau)+alphaTS{1}(ktau)+gammaTS{1}(ktau)+muh);
        Phi15=qs*qT{1}*omegaU_T(ktau);
        %
        Phi22=-((1-qT{2})*omegaTS_U{2}(ktau)+alphaTS{2}(ktau)+gammaTS{2}(ktau)+muh);
        Phi25=qs*qT{2}*omegaU_T(ktau);
        %
        Phi33=-((1-qT{1})*omegaTF_U{1}(ktau)+alphaTF{1}(ktau)+muh);
        Phi35=(1-qs)*qT{1}*omegaU_T(ktau);
        %
        Phi44=-((1-qT{2})*omegaTF_U{2}(ktau)+alphaTF{2}(ktau)+muh);
        Phi45=(1-qs)*qT{2}*omegaU_T(ktau);
        %
        Phi51=(1-qT{1})*omegaTS_U{1}(ktau);
        Phi52=(1-qT{2})*omegaTS_U{2}(ktau);
        Phi53=(1-qT{1})*omegaTF_U{1}(ktau);
        Phi54=(1-qT{2})*omegaTF_U{2}(ktau);
        Phi55=-((qT{1}+qT{2})*omegaU_T(ktau)+alphaU(ktau)+muh);
        %
        Phi(:,:,ktau)=[Phi11 0 0 0 Phi15;...
            0 Phi22 0 0 Phi25;...
            0 0 Phi33 0 Phi35;...
            0 0 0 Phi44 Phi45;...
            Phi51 Phi52 Phi53 Phi54 Phi55];

        Id=1:ktau; 
        Chi(ktau)=[betaTS{1}(ktau) betaTS{2}(ktau) betaTF{1}(ktau) betaTF{2}(ktau)...
            betaU(ktau)]*(expm(trapz(Tau(Id),Phi(:,:,Id),3))*QQ);
    end
    CHI=trapz(Tau,Chi); 
    R0=(Lambdah/muh)*CHI;  
%% Interval and discretization of time t
    dt=1; % discretization step
    T=0:dt:Tmax;  % Time Vector 
    Nt=length(T);  % vector length of t

%% THE BETWEEN-HOST MODEL
% Initialization of state variables of the between-host model
    S=zeros(Nt,1);  IU=zeros(Ntau,Nt); 
    ITF=cell(1,NbTreat); ITS=cell(1,NbTreat); R=cell(1,NbTreat);
     
    for l=1:NbTreat
        ITF{l}=zeros(Ntau,Nt); ITS{l}=zeros(Ntau,Nt);
        R{l}=zeros(Nt,1);
    end
        
% Initial conditions for S, RA and RB
    S(1)=InitSR(1); R{1}(1)=InitSR(2); R{2}(1)=InitSR(3);

% Initial conditions for infectious states [ITFA ITSA ITFB ITSB IU]
    ITF{1}(:,1)=InitInfected(1:Ntau); 
    ITS{1}(:,1)=InitInfected(Ntau+1:2*Ntau);
    ITF{2}(:,1)=InitInfected(2*Ntau+1:3*Ntau);
    ITS{2}(:,1)=InitInfected(3*Ntau+1:4*Ntau);
    IU(:,1)=InitInfected(4*Ntau+1:5*Ntau);
        
% Solve of the between-host model
    for t=1:Nt-1           
        lambda=betaTS{1}(:).*ITS{1}(:,t)+betaTF{1}(:).*ITF{1}(:,t)+...
             betaTS{2}(:).*ITS{2}(:,t)+betaTF{2}(:).*ITF{2}(:,t)+...
             betaU(:).*IU(:,t);

        lambda_t=trapz(Tau,lambda); 
        S(t+1)=(S(t)+dt*Lambdah)/(1+dt*muh+dt*lambda_t);
       for l=1:NbTreat
          Newinfect=S(t)*lambda_t;
          ITS{l}(1,t+1)=qs*qT{l}*Newinfect; 
          ITF{l}(1,t+1)=(1-qs)*qT{l}*Newinfect;  
          IU(1,t+1)=(1-qT{1}-qT{2})*Newinfect; 

          RR=zeros(Ntau,1);
          RR(1)=gammaTS{l}(1)*ITS{l}(1,t);

          for n=2:Ntau

            ITS{l}(n,t+1)=((ITS{l}(n,t)/dt)+(ITS{l}(n-1,t+1)/dtau)+qs*qT{l}*omegaU_T(n)*IU(n,t))/...
                ( 1/dt + 1/dtau + muh+alphaTS{l}(n)+(1-qT{l})*omegaTS_U{l}(n)+gammaTS{l}(n));

            ITF{l}(n,t+1)=( (ITF{l}(n,t)/dt) + (ITF{l}(n-1,t+1)/dtau) + (1-qs)*qT{l}*omegaU_T(n)*IU(n,t))/...
                ( 1/dt + 1/dtau + muh+alphaTF{l}(n)+(1-qT{l})*omegaTF_U{l}(n));

            IU(n,t+1)=( (IU(n,t)/dt) + (IU(n-1,t+1)/dtau) + (1-qT{1})*omegaTS_U{1}(n)*ITS{1}(n,t)...
                +(1-qT{1})*omegaTF_U{1}(n)*ITF{1}(n,t)+(1-qT{2})*omegaTS_U{2}(n)*ITS{2}(n,t)...
                +(1-qT{2})*omegaTF_U{2}(n)*ITF{2}(n,t))/( 1/dt + 1/dtau + muh+alphaU(n)+(qT{1}+qT{2})*omegaU_T(n));

            RR(n)=gammaTS{l}(n)*ITS{l}(n,t);
         end
            R{l}(t+1)=(R{l}(t)+dt*trapz(Tau,RR))/(1+dt*muh); 
       end
    end

    SR=[S R{1} R{2}];
    Infected=[ITF{1} ITS{1} ITF{2} ITS{2} IU];
%% Proportion of individuals 
    TotPopS=zeros(Nt,1); TotPopR=zeros(Nt,1); 
    TotPopITS=zeros(Nt,1); TotPopITF=zeros(Nt,1); TotPopIU=zeros(Nt,1); 
    TotPopI=zeros(Nt,1); TotPop=zeros(Nt,1); 

    TotPopITSA=zeros(Nt,1); TotPopITSB=zeros(Nt,1);
    TotPopITFA=zeros(Nt,1); TotPopITFB=zeros(Nt,1);
    TotPopRA=zeros(Nt,1); TotPopRB=zeros(Nt,1);
     %
     
   for t=1:Nt
       TotPopS(t)=TotPopS(t)+S(t);
       TotPopRA(t)=TotPopRA(t)+R{1}(t);
       TotPopRB(t)=TotPopRB(t)+R{2}(t);
       TotPopR(t)=TotPopRA(t)+TotPopRB(t);
       %
       TotPopITSA(t)=trapz(Tau,ITS{1}(:,t));
       TotPopITSB(t)=trapz(Tau,ITS{2}(:,t));
       TotPopITFA(t)=trapz(Tau,ITF{1}(:,t));
       TotPopITFB(t)=trapz(Tau,ITF{2}(:,t));
       TotPopIU(t)=trapz(Tau,IU(:,t));    
       %       
       TotPopITS(t)=TotPopITSA(t)+TotPopITSB(t);
       TotPopITF(t)=TotPopITFA(t)+TotPopITFB(t);
       %
       TotPopI(t)=TotPopITS(t)+TotPopITF(t)+TotPopIU(t);
       TotPop(t)=TotPopS(t)+TotPopR(t)+TotPopI(t);
   end
     
  % Proportion of each compartment (in the total population)
      PropS=TotPopS./TotPop;
      PropR=TotPopR./TotPop;
      
      PropTS=TotPopITS./TotPop;
      PropTF=TotPopITF./TotPop;
      PropU=TotPopIU./TotPop;
      
      PropCompartment=[PropS PropR PropTS PropTF PropU];
      
      % Proportion of treatment failure (in the total treated individuals)
      PropITFA=TotPopITFA./(TotPopITSA+TotPopITFA+TotPopRA+1);
      PropITFB=TotPopITFB./(TotPopITSB+TotPopITFB+TotPopRB+1);
      PropITF=TotPopITF./(TotPopITS+TotPopITF+TotPopR+1);
      PropTreatFail=[PropITF PropITFA PropITFB];
      PropITF_equi=PropITF(length(PropITF));
      
%%  Average level of resistance in the host population at final time T
      EtaTf_t=zeros(1,Nt); 
      for t=1:Nt
        EtaTf_t(t)=trapz(Tau,xBarTF{1}(Ntau).*ITF{1}(:,t)+xBarTF{2}(Ntau).*ITF{2}(:,t)...
        +xBarTS{1}(Ntau).*ITS{1}(:,t)+xBarTS{2}(Ntau).*ITS{2}(:,t)...
        +xBarU(Ntau).*IU(:,t))/TotPopI(t);
      end
      EtaTf=EtaTf_t(Nt);
end  


 %% Fonction within-host Model
function[x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,beta,alpha,xBar]=...
     MainWithinhost(Tau,dtau,Ntau,TreatStatus,Drug,T_treatmentStart,T_Recovery,N0BEtaStar,N0AEtaStar,EtaStar)
% Discretization
    dx=0.01; % discretization step
    xmin=-0.5; xmax=2.5; % Resistance space 
    x=xmin:dx:xmax; % discretization of resistance level x
    Nx=length(x); % vector length of x

% Fix parameters 
    pm=10; p0=0.95*pm; kappa=1; m0=0.05; sigma0=0.05; p1_0=0.5;
    p1=p0*p1_0;
    r1=9e3; beta0=1.848e-3/1.5; alpha0=7.5e-2; % For alpha and beta

% immune effect
    N00=21;
    mu=1.81*p0/N00;
  
%functional parameters p
    pNum=1+(pm/p0-1)*(p0*(pm-p1)/(p1*(pm-p0))).^x;
    pp=pm./pNum;
 
%functional parameter k
    if strcmp(TreatStatus,'U')
         eps=0.05/5;
         kk=zeros(1,Nx); 
    elseif strcmp(TreatStatus,'TS')
        eps=0.05;
        if strcmp(Drug,'A')
            k0=15; k1_0=0.3;  k1=k0*k1_0; % TS
            kk=k0*(k1/k0).^(x);
        elseif strcmp(Drug,'B')
            k0=15; k1_0=0.3;  k1=k0*k1_0; % TS
            kk=((mu*(1-N0BEtaStar/N0AEtaStar)+(k0*(k1/k0)^EtaStar))/(N0BEtaStar/N0AEtaStar))*...
                (k1/k0).^(x-EtaStar);
        end
    elseif strcmp(TreatStatus,'TF')
        eps=0.05;
        if strcmp(Drug,'A')
            k0=3; k1_0=0.01;  k1=k0*k1_0;  
            kk=k0*(k1/k0).^(x);
        elseif strcmp(Drug,'B')
            k0=3; k1_0=0.01;  k1=k0*k1_0;
            kk=((mu*(1-N0BEtaStar/N0AEtaStar)+(k0*(k1/k0)^EtaStar))/(N0BEtaStar/N0AEtaStar))*...
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

%% Solve of the within-host model
    for tau=1:Ntau-1
      B=(1+totpopb(tau))^(-kappa);
      LT=B*dx*conv(pp.*b(tau,:),J0,'same');    
      b(tau+1,:)=(b(tau,:)+dtau*LT)./(1+dtau*mu+dtau*kk); 
      bb=b(tau+1,:);
      totpopb(tau+1)=sum(bb);
    end
   
% Individual average level of resistance
    xBar=zeros(1,Ntau);
    for n=1:Ntau
    xBar(n)=sum(x.*b(n,:))/totpopb(n);
    end
  
    IndexT_Recovery=find(Tau>0 & abs(Tau-T_Recovery)<0.0001);
    IndexT_treatmentStart=find(Tau>0 & abs(Tau-T_treatmentStart)<0.0001);
    ds=totpopb(IndexT_treatmentStart)/totpopb(1);
    SeuilGuerison=totpopb(IndexT_Recovery)/totpopb(1);

% Transition rates of between-host model
    kz=100;
    omegaTU=kz*(totpopb>totpopb(1));
    omegaUT=kz*(totpopb>ds*totpopb(1));
    gamma=kz*(totpopb<SeuilGuerison*totpopb(1));
    
    % Transmission and disease-induced mortality rates
    beta=(beta0*totpopb)./(r1+totpopb);
    alpha=(alpha0*totpopb)./(r1+totpopb);
end


end

