function[R0_AB,Re,PropITF_equiAB,R0_A,PropITF_equiA,x,T_AB,...
    PropCompartment_AB,PropTreatFail_AB,N0TFA,N0TFB,...
    EtaStar,xEtaStar,SR_AB,Infected_AB,TotPopDynamics_AB,bart_year,tTF_year]=...
    MainFunctions...
    (qs,T_treatmentStart,T_Recovery,PropTreated,RelEfficacy,PropTreatB,TmaxAB,...
    DImmunity)
%% Fixed parameters
    muh = 1/(70*365);          % natural human mortality rate (/day)
    Lambdah = 30754/365;       % recruitment rate (/day)
    DFE = Lambdah/muh;

%% Dynamics with only one treatment qTA=1

    % Age vector to define initial conditions for infectious states
    dtau=1;%0.1/20; 
    Taumax=30; % Maximal time of infection age
    Tau=0:dtau:Taumax;
    Ntau=length(Tau); % vector length of infection age
    
    deta=dtau;
    Teta=0:deta:DImmunity;      % time since recovery (immune age)
    Neta=length(Teta);

%% Equilibrium with only A
    
    [Sstar,Rstar,ITSAstar,ITFAstar,IUstar,...
        Jstar,lambdaStar,R0_A,EtaStar,EtaMean,...
        PropITF_equiA,N0TFA,x] = ...
        Equilibrium_A(Tau,dtau,PropTreated,qs,...
        T_treatmentStart,T_Recovery,...
        muh,Teta,deta,Lambdah,Neta,Ntau,DFE);
    
    nx = length(x);

%% Dynamics with two treatments AB
    % Initial conditions (correspond to the equilibrium with A)
    
    InitSR_AB = ...
    [Sstar Rstar zeros(1,Neta)];

    InitInfected_AB = ...
        [ITFAstar' ITSAstar' ...
         zeros(1,Ntau) zeros(1,Ntau) ...
         IUstar'];
    
    xEtaStar = find(abs(x-round(EtaStar,2))<0.0001);
    N0AEtaStar = N0TFA(xEtaStar);
    N0BEtaStar = (1-RelEfficacy)*N0AEtaStar;

    [x,T_AB,N0_AB,R0_AB,Re,SR_AB,Infected_AB,TotPopDynamics_AB,PropCompartment_AB,...
        PropTreatFail_AB,EtaTf_AB,PropITF_equiAB,EtaTf_tAB,Eta_tAB,bart_year,tTF_year] = ...
        MainBetweenhost(TmaxAB,Tau,dtau,Ntau,PropTreated,qs,PropTreatB,...
        T_treatmentStart,T_Recovery,muh,Teta,Neta,deta,Lambdah,InitSR_AB,InitInfected_AB,...
        N0BEtaStar,N0AEtaStar,EtaStar);
   
   %N0_AB=[N0TF,A​,N0TS,A​,N0TF,B​,N0TS,B​,N0U​] == structure of N0_AB,
   N0TFB=N0_AB(2*nx+1:3*nx);
     
    %% Function Between-host model

function [x,T,N0_WH,R0,Re,SR,Infected,TotPopDynamics,...
          PropCompartment,PropTreatFail,EtaTf,PropITF_equi,...
          EtaTf_t,Eta_t,bart_year,tTF_year] = ...
    MainBetweenhost(Tmax,Tau,dtau,Ntau,PropTreated,qs,PropTreatB,...
                    T_treatmentStart,T_Recovery,muh,Teta,Neta,deta,Lambdah,...
                    InitSR,InitInfected,...
                    N0BEtaStar,N0AEtaStar,EtaStar)

VectDrugGroup = {'A','B'};
NbTreat = length(VectDrugGroup);

QT = PropTreated .* [1-PropTreatB PropTreatB];

qT = cell(1,NbTreat);
qT{1} = QT(1);
qT{2} = QT(2);

qTot = qT{1} + qT{2};
qU   = 1 - qTot;


%% ================================================================
%  PARAMETERS FROM WITHIN-HOST MODEL
% ================================================================

betaTF   = cell(1,NbTreat);
betaTS   = cell(1,NbTreat);

alphaTF  = cell(1,NbTreat);
alphaTS  = cell(1,NbTreat);

gammaTS  = cell(1,NbTreat);

omegaTF_U = cell(1,NbTreat);
omegaTS_U = cell(1,NbTreat);

xBarTF = cell(1,NbTreat);
xBarTS = cell(1,NbTreat);

N0TF = cell(1,NbTreat);
N0TS = cell(1,NbTreat);


%% ------------------------------------------------
% Untreated infection U
% IMPORTANT:
% U does not depend on Drug A/B in MainWithinhost.
%% ------------------------------------------------

[x,dx,Nx,b_u,totpopb_u,N0_u,...
 omegaTU_u,omegaUT_u,gamma_u,...
 beta_u,alpha_u,xBar_u] = ...
    MainWithinhost(Tau,dtau,Ntau,'U','NA',...
                   T_treatmentStart,T_Recovery,...
                   N0BEtaStar,N0AEtaStar,EtaStar);

omegaU_T = omegaUT_u;

betaU  = beta_u;
alphaU = alpha_u;
xBarU  = xBar_u;
N0U    = N0_u;


%% ------------------------------------------------
% Treated infections: A and B
%% ------------------------------------------------

for l = 1:NbTreat

    Drug = VectDrugGroup{l};

    % Treatment failure
    [~,~,~,~,~,N0_nc,...
     omegaTU_nc,~,~,beta_nc,alpha_nc,xBar_nc] = ...
        MainWithinhost(Tau,dtau,Ntau,'TF',Drug,...
                       T_treatmentStart,T_Recovery,...
                       N0BEtaStar,N0AEtaStar,EtaStar);

    omegaTF_U{l} = omegaTU_nc;
    betaTF{l}    = beta_nc;
    alphaTF{l}   = alpha_nc;
    xBarTF{l}    = xBar_nc;
    N0TF{l}      = N0_nc;


    % Treatment success
    [~,~,~,~,~,N0_c,...
     omegaTU_c,~,gamma_c,beta_c,alpha_c,xBar_c] = ...
        MainWithinhost(Tau,dtau,Ntau,'TS',Drug,...
                       T_treatmentStart,T_Recovery,...
                       N0BEtaStar,N0AEtaStar,EtaStar);

    omegaTS_U{l} = omegaTU_c;
    gammaTS{l}   = gamma_c;
    betaTS{l}    = beta_c;
    alphaTS{l}   = alpha_c;
    xBarTS{l}    = xBar_c;
    N0TS{l}      = N0_c;

end


% Vector of within-host N0
% {TF-A, TS-A, TF-B, TS-B, U}
N0_WH = [N0TF{1} N0TS{1} N0TF{2} N0TS{2} N0U];


%% ================================================================
% BASIC REPRODUCTION NUMBER R0
%% ================================================================

QQ = [qs*qT{1};...
      qs*qT{2};...
      (1-qs)*qT{1};...
      (1-qs)*qT{2};...
      qU];

Chi = zeros(1,Ntau);

% Integral_0^tau Phi(s) ds
IntPhi = zeros(5,5);

Phi_previous = zeros(5,5);


for ktau = 1:Ntau

    % ------------------------------------------------------------
    % Phi(tau)
    % ------------------------------------------------------------

    Phi11 = -((1-qT{1})*omegaTS_U{1}(ktau) ...
              + alphaTS{1}(ktau) ...
              + gammaTS{1}(ktau) ...
              + muh);

    Phi15 = qs*qT{1}*omegaU_T(ktau);


    Phi22 = -((1-qT{2})*omegaTS_U{2}(ktau) ...
              + alphaTS{2}(ktau) ...
              + gammaTS{2}(ktau) ...
              + muh);

    Phi25 = qs*qT{2}*omegaU_T(ktau);


    Phi33 = -((1-qT{1})*omegaTF_U{1}(ktau) ...
              + alphaTF{1}(ktau) ...
              + muh);

    Phi35 = (1-qs)*qT{1}*omegaU_T(ktau);


    Phi44 = -((1-qT{2})*omegaTF_U{2}(ktau) ...
              + alphaTF{2}(ktau) ...
              + muh);

    Phi45 = (1-qs)*qT{2}*omegaU_T(ktau);


    Phi51 = (1-qT{1})*omegaTS_U{1}(ktau);
    Phi52 = (1-qT{2})*omegaTS_U{2}(ktau);
    Phi53 = (1-qT{1})*omegaTF_U{1}(ktau);
    Phi54 = (1-qT{2})*omegaTF_U{2}(ktau);

    Phi55 = -(qTot*omegaU_T(ktau) ...
              + alphaU(ktau) ...
              + muh);


    Phi_current = ...
        [Phi11  0      0      0      Phi15;...
         0      Phi22  0      0      Phi25;...
         0      0      Phi33  0      Phi35;...
         0      0      0      Phi44  Phi45;...
         Phi51  Phi52  Phi53  Phi54  Phi55];


    % ------------------------------------------------------------
    % Cumulative trapezoidal integration
    % ------------------------------------------------------------

    if ktau > 1

        hTau = Tau(ktau) - Tau(ktau-1);

        IntPhi = IntPhi ...
               + 0.5*hTau*(Phi_previous + Phi_current);

    end


    BetaVector = ...
        [betaTS{1}(ktau),...
         betaTS{2}(ktau),...
         betaTF{1}(ktau),...
         betaTF{2}(ktau),...
         betaU(ktau)];


    Chi(ktau) = BetaVector * (expm(IntPhi)*QQ);


    Phi_previous = Phi_current;

end


CHI = trapz(Tau,Chi);

R0 = (Lambdah/muh)*CHI;


%% ================================================================
% TIME DISCRETIZATION and STATE VARIABLES
%% ================================================================
dt = dtau;
T  = 0:dt:Tmax;
Nt = length(T);

S  = zeros(Nt,1);
IU = zeros(Ntau,Nt);

ITF = cell(1,NbTreat);
ITS = cell(1,NbTreat);
R   = cell(1,NbTreat);

for l = 1:NbTreat

    ITF{l} = zeros(Ntau,Nt);
    ITS{l} = zeros(Ntau,Nt);
    R{l} = zeros(Nt,Neta);

end


%% Population totals

TotPopR   = zeros(Nt,1);
TotPopITS = zeros(Nt,1);
TotPopITF = zeros(Nt,1);
TotPopIU  = zeros(Nt,1);

TotPopI   = zeros(Nt,1);
TotPop    = zeros(Nt,1);


%% ================================================================
% INITIAL CONDITIONS
%% ================================================================

S(1) = InitSR(1);

% Initial immune-age distributions
R{1}(1,:) = InitSR(2:Neta+1);
R{2}(1,:) = InitSR(Neta+2:2*Neta+1);


% [ITFA ITSA ITFB ITSB IU]

ITF{1}(:,1) = InitInfected(1:Ntau);

ITS{1}(:,1) = ...
    InitInfected(Ntau+1:2*Ntau);

ITF{2}(:,1) = ...
    InitInfected(2*Ntau+1:3*Ntau);

ITS{2}(:,1) = ...
    InitInfected(3*Ntau+1:4*Ntau);

IU(:,1) = ...
    InitInfected(4*Ntau+1:5*Ntau);

TotPopR(1) = ...
    trapz(Teta,R{1}(1,:)) + trapz(Teta,R{2}(1,:));

TotPopITS(1) = ...
    trapz(Tau,ITS{1}(:,1)) ...
    + trapz(Tau,ITS{2}(:,1));

TotPopITF(1) = ...
    trapz(Tau,ITF{1}(:,1)) ...
    + trapz(Tau,ITF{2}(:,1));

TotPopIU(1) = ...
    trapz(Tau,IU(:,1));


%% ================================================================
% PRECOMPUTE CONSTANT DENOMINATORS
% They do not depend on time t.
%% ================================================================

DenITS = cell(1,NbTreat);
DenITF = cell(1,NbTreat);

for l = 1:NbTreat

    DenITS{l} = ...
        1/dt + 1/dtau ...
        + muh ...
        + alphaTS{l}(:) ...
        + qU*omegaTS_U{l}(:) ...
        + gammaTS{l}(:);

    DenITF{l} = ...
        1/dt + 1/dtau ...
        + muh ...
        + alphaTF{l}(:) ...
        + qU*omegaTF_U{l}(:);

end

DenIU = ...
    1/dt + 1/dtau ...
    + muh ...
    + alphaU(:) ...
    + qTot*omegaU_T(:);

DenR = 1/dt + 1/deta + muh;

%% ================================================================
% BETWEEN-HOST DYNAMICS
%% ================================================================

for t = 1:Nt-1

    %% ------------------------------------------------------------
    % Force of infection
    %% ------------------------------------------------------------

    lambda = ...
          betaTS{1}(:).*ITS{1}(:,t) ...
        + betaTF{1}(:).*ITF{1}(:,t) ...
        + betaTS{2}(:).*ITS{2}(:,t) ...
        + betaTF{2}(:).*ITF{2}(:,t) ...
        + betaU(:).*IU(:,t);


    lambda_t = trapz(Tau,lambda);


    %% ------------------------------------------------------------
    % New infections generated at time t
    %% ------------------------------------------------------------
    Newinfect = S(t)*lambda_t;

    %% ------------------------------------------------------------
    % Recovery fluxes
    %% ------------------------------------------------------------

    RecoveryFluxA = ...
        trapz(Tau,...
              gammaTS{1}(:).*ITS{1}(:,t));

    RecoveryFluxB = ...
        trapz(Tau,...
              gammaTS{2}(:).*ITS{2}(:,t));


    %% ------------------------------------------------------------
    % Infection-age dynamics
    %% ------------------------------------------------------------
    
    %% Bundary 
    n=1; 
    %Treatment A - success
    ITS{1}(n,t+1) = ...
            ( ITS{1}(n,t)/dt ...
              + qs*qT{1}*Newinfect/dtau ...
              + qs*qT{1}*omegaU_T(n)*IU(n,t) ) ...
            / DenITS{1}(n);
        
    % Treatment B - success
    ITS{2}(n,t+1) = ...
        ( ITS{2}(n,t)/dt ...
          + qs*qT{2}*Newinfect/dtau ...
          + qs*qT{2}*omegaU_T(n)*IU(n,t) ) ...
        / DenITS{2}(n);
    
    % Treatment A - failure
    ITF{1}(n,t+1) = ...
        ( ITF{1}(n,t)/dt ...
          + (1-qs)*qT{1}*Newinfect/dtau ...
          + (1-qs)*qT{1}*omegaU_T(n)*IU(n,t) ) ...
        / DenITF{1}(n);
    
    % Treatment B - failure
    ITF{2}(n,t+1) = ...
        ( ITF{2}(n,t)/dt ...
          + (1-qs)*qT{2}*Newinfect/dtau ...
          + (1-qs)*qT{2}*omegaU_T(n)*IU(n,t) ) ...
        / DenITF{2}(n);
    
    % Untreated
    IU(n,t+1) = ...
            ( IU(n,t)/dt ...
              + qU*Newinfect/dtau ...
              + qU*omegaTS_U{1}(n)*ITS{1}(n,t) ...
              + qU*omegaTF_U{1}(n)*ITF{1}(n,t) ...
              + qU*omegaTS_U{2}(n)*ITS{2}(n,t) ...
              + qU*omegaTF_U{2}(n)*ITF{2}(n,t) ) ...
            / DenIU(n);
    

    for n = 2:Ntau

        %% Treatment A - success

        ITS{1}(n,t+1) = ...
            ( ITS{1}(n,t)/dt ...
              + ITS{1}(n-1,t+1)/dtau ...
              + qs*qT{1}*omegaU_T(n)*IU(n,t) ) ...
            / DenITS{1}(n);


        %% Treatment B - success

        ITS{2}(n,t+1) = ...
            ( ITS{2}(n,t)/dt ...
              + ITS{2}(n-1,t+1)/dtau ...
              + qs*qT{2}*omegaU_T(n)*IU(n,t) ) ...
            / DenITS{2}(n);


        %% Treatment A - failure

        ITF{1}(n,t+1) = ...
            ( ITF{1}(n,t)/dt ...
              + ITF{1}(n-1,t+1)/dtau ...
              + (1-qs)*qT{1}*omegaU_T(n)*IU(n,t) ) ...
            / DenITF{1}(n);


        %% Treatment B - failure

        ITF{2}(n,t+1) = ...
            ( ITF{2}(n,t)/dt ...
              + ITF{2}(n-1,t+1)/dtau ...
              + (1-qs)*qT{2}*omegaU_T(n)*IU(n,t) ) ...
            / DenITF{2}(n);


        %% Untreated

        IU(n,t+1) = ...
            ( IU(n,t)/dt ...
              + IU(n-1,t+1)/dtau ...
              + qU*omegaTS_U{1}(n)*ITS{1}(n,t) ...
              + qU*omegaTF_U{1}(n)*ITF{1}(n,t) ...
              + qU*omegaTS_U{2}(n)*ITS{2}(n,t) ...
              + qU*omegaTF_U{2}(n)*ITF{2}(n,t) ) ...
            / DenIU(n);

    end


    %% ------------------------------------------------------------
    % Recovered classes structured by immune age
    %% ------------------------------------------------------------    
    eta=1;
    R{1}(t+1,eta)= ...
            (R{1}(t,eta)/dt+ RecoveryFluxA/deta)...
            / DenR;
        
    R{2}(t+1,eta)= ...
        (R{2}(t,eta)/dt+ RecoveryFluxB/deta)...
        / DenR; 
    
    for eta=2:Neta
        R{1}(t+1,eta)= ...
            (R{1}(t,eta)/dt+ R{1}(t+1,eta-1)/deta)...
            / DenR;
        
        R{2}(t+1,eta)= ...
            (R{2}(t,eta)/dt+ R{2}(t+1,eta-1)/deta)...
            / DenR;
    end

    %% ------------------------------------------------------------
    % Susceptibles
    %% ------------------------------------------------------------
    ImmunityLossFlux = R{1}(t+1,end) + R{2}(t+1,end);

    S(t+1) = ...
        ( S(t) ...
          + dt*Lambdah ...
          + dt*ImmunityLossFlux ) ...
        / ...
        (1 + dt*muh + dt*lambda_t);


    %% ------------------------------------------------------------
    % Totals
    %% ------------------------------------------------------------

    TotPopR(t+1) = ...
       trapz(Teta,R{1}(t+1,:)) + trapz(Teta,R{2}(t+1,:));


    TotPopITS(t+1) = ...
        trapz(Tau,ITS{1}(:,t+1)) ...
        + trapz(Tau,ITS{2}(:,t+1));


    TotPopITF(t+1) = ...
        trapz(Tau,ITF{1}(:,t+1)) ...
        + trapz(Tau,ITF{2}(:,t+1));


    TotPopIU(t+1) = ...
        trapz(Tau,IU(:,t+1));

end


%% ================================================================
% OUTPUT STATE VARIABLES
%% ================================================================

SR = [S R{1} R{2}];

Re = S*CHI;

Infected = ...
    [ITF{1} ITS{1} ITF{2} ITS{2} IU];


TotPopDynamics = cell(1,5);

TotPopDynamics{1} = S;
TotPopDynamics{2} = TotPopR;
TotPopDynamics{3} = TotPopITS;
TotPopDynamics{4} = TotPopITF;
TotPopDynamics{5} = TotPopIU;


%% ================================================================
% PROPORTIONS
%% ================================================================

TotPopITSA = zeros(Nt,1);
TotPopITSB = zeros(Nt,1);

TotPopITFA = zeros(Nt,1);
TotPopITFB = zeros(Nt,1);

TotPopRA = zeros(Nt,1);
TotPopRB = zeros(Nt,1);


for t = 1:Nt

    TotPopITSA(t) = ...
        trapz(Tau,ITS{1}(:,t));

    TotPopITSB(t) = ...
        trapz(Tau,ITS{2}(:,t));

    TotPopITFA(t) = ...
        trapz(Tau,ITF{1}(:,t));

    TotPopITFB(t) = ...
        trapz(Tau,ITF{2}(:,t));

    TotPopRA(t) = ...
        trapz(Teta,R{1}(t,:));

    TotPopRB(t) = ...
        trapz(Teta,R{2}(t,:));


    TotPopI(t) = ...
        TotPopITS(t) ...
        + TotPopITF(t) ...
        + TotPopIU(t);


    TotPop(t) = ...
        S(t) ...
        + TotPopR(t) ...
        + TotPopI(t);

end


%% Compartment proportions

PropS  = S./TotPop;

PropR  = TotPopR./TotPop;

PropTS = TotPopITS./TotPop;

PropTF = TotPopITF./TotPop;

PropU  = TotPopIU./TotPop;


PropCompartment = ...
    [PropS PropR PropTS PropTF PropU];


%% Treatment failure proportions
%% ---------------------------------------------------------------

PropITFA = ...
    TotPopITFA ./ ...
    (TotPopITSA + TotPopITFA + TotPopRA+0.05);


PropITFB = ...
    TotPopITFB ./ ...
    (TotPopITSB + TotPopITFB + TotPopRB+0.05);


PropITF = ...
    TotPopITF ./ ...
    (TotPopITS + TotPopITF + TotPopR+0.05);


PropTreatFail = ...
    [PropITF PropITFA PropITFB];


PropITF_equi = ...
    PropITF(end);


%% ================================================================
% AVERAGE RESISTANCE LEVEL
%% ================================================================

EtaTf_t = zeros(1,Nt);
Eta_t   = zeros(1,Nt);


for t = 1:Nt

    if TotPopI(t) > 0

        EtaTf_t(t) = ...
            trapz(Tau,...
                xBarTF{1}(Ntau).*ITF{1}(:,t) ...
              + xBarTF{2}(Ntau).*ITF{2}(:,t) ...
              + xBarTS{1}(Ntau).*ITS{1}(:,t) ...
              + xBarTS{2}(Ntau).*ITS{2}(:,t) ...
              + xBarU(Ntau).*IU(:,t)) ...
            / TotPopI(t);


        Eta_t(t) = ...
            trapz(Tau,...
                xBarTF{1}(:).*ITF{1}(:,t) ...
              + xBarTF{2}(:).*ITF{2}(:,t) ...
              + xBarTS{1}(:).*ITS{1}(:,t) ...
              + xBarTS{2}(:).*ITS{2}(:,t) ...
              + xBarU(:).*IU(:,t)) ...
            / TotPopI(t);

    else

        EtaTf_t(t) = NaN;
        Eta_t(t)   = NaN;

    end

end


EtaTf = EtaTf_t(end);

    %% =================================================
    %% strategy durability time OR duration of epidemiological benefit,
    %% =================================================
%     Sref = PropS(1);
%     tol = 1e-6;
%     Sdiff = PropS(:) - Sref*(1+EpiBenefThreshold); 
%     % Indices of downward crossings of Sref
%     idx_cross = find( ...
%         Sdiff(1:end-1) >= -tol & ...
%         Sdiff(2:end)   <  -tol );
% 
%     if length(idx_cross) >= 2
% 
%         idx_bart = idx_cross(2) + 1;
%         bart_year = T(idx_bart)/365;       % years
% 
%     else
%         bart_year = T(end)/365;
% 
%     end
    
    %Begin New definition
    ThresholdTF = 1e-10;
    % Last time at which PropS is still below min(0.95, PropS(1)*2)
    Sdiff = min(0.95, PropS(1)*2) -PropS(:);
    idx_last_below = find(Sdiff >= ThresholdTF,1,'last');
    
    if isempty(idx_last_below)
        bart_year = T(1)/365;
        
    elseif idx_last_below == length(PropS)
        
        % PropS does not reach and remain above the threshold
        % during the simulation
        bart_year = T(end)/365;
        
    else
        % First time after the last threshold exceedance
        idx_tTF = idx_last_below + 1;
        bart_year = T(idx_tTF)/365;        % years
        
    end
    %End New definition 

    
    %% =================================================
    %% time to quasi-elimination of treatment-failure infections
    %% =================================================
    ThresholdTF = 1e-10;
    % Last time at which PropTF is still above the threshold
    idx_last_above = find(PropTF(:) >= ThresholdTF,1,'last');
    
    if isempty(idx_last_above)
        tTF_year = T(1)/365;
        
    elseif idx_last_above == length(PropTF)
        
        % PropTF has not permanently fallen below the threshold
        % during the simulation
        tTF_year = T(end)/365;
        
    else
        % First time after the last threshold exceedance
        idx_tTF = idx_last_above + 1;
        tTF_year = T(idx_tTF)/365;        % years
        
    end
    

end


 %% Fonction within-host Model
function[x,dx,Nx,b,totpopb,N0,omegaTU,omegaUT,gamma,beta,alpha,xBar]=...
     MainWithinhost(Tau,dtau,Ntau,TreatStatus,Drug,T_treatmentStart,...
     T_Recovery,N0BEtaStar,N0AEtaStar,EtaStar)
% Discretization
    dx=0.01; % discretization step
    xmin=-0.5; xmax=2.5; % Resistance space 
    x=xmin:dx:xmax; % discretization of resistance level x
    Nx=length(x); % vector length of x

% Fix parameters 
    pm=10; p0=0.95*pm; kappa=1; m0=0.05; sigma0=0.05; p1_0=0.5;
    p1=p0*p1_0;
    r1=9e3; beta0=2*((1.848e-3/1.5)/2)/3.5; alpha0=7.5e-2; % For alpha and beta

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
  
    % Robust mapping of continuous event times to the numerical Tau grid
    if T_Recovery < Tau(1) || T_Recovery > Tau(end) || ...
       T_treatmentStart < Tau(1) || T_treatmentStart > Tau(end)
        error('T_Recovery and T_treatmentStart must lie inside the Tau interval.');
    end

    [~,IndexT_Recovery] = min(abs(Tau-T_Recovery));
    [~,IndexT_treatmentStart] = min(abs(Tau-T_treatmentStart));

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


    %% the equilibrium function
    function [Sstar,Rstar,ITSAstar,ITFAstar,IUstar,...
              Jstar,lambdaStar,R0eq,EtaStar,EtaMean,...
              PropITFstar,N0TFA,x] = ...
        Equilibrium_A(Tau,dtau,PropTreated,qs,...
                      T_treatmentStart,T_Recovery,...
                      muh,Teta,deta,Lambdah,Neta,Ntau,DFE)
        % EQUILIBRIUM_A
        % Direct computation of the positive endemic equilibrium for treatment A only.
        %
        
        % Outputs
        % -------
        % Sstar       : susceptible population at endemic equilibrium
        % Rstar       : immune-age density R_A^*(eta), row vector
        % ITSAstar    : I_TS,A^*(tau), column vector
        % ITFAstar    : I_TF,A^*(tau), column vector
        % IUstar      : I_U^*(tau), column vector
        % Jstar       : equilibrium incidence S^* lambda^*
        % lambdaStar  : equilibrium force of infection
        % R0eq        : reproduction number associated with the stationary discrete scheme
        % EtaStar     : same "terminal within-host resistance" quantity used by EtaTf_t
        % EtaMean     : resistance mean using xBar(tau), analogous to Eta_t
        % PropITFstar : equilibrium treatment-failure proportion used in current code
        % N0TFA       : within-host N0(x) for treatment failure under A
        % x           : resistance grid
        %
        % A positive equilibrium exists only if R0eq > 1.

        %% ---------------------------------------------------------------
        % Treatment probabilities: A only
        %% ---------------------------------------------------------------
        qA = PropTreated;
        qU = 1 - qA;

        %% ---------------------------------------------------------------
        % Within-host quantities
        %% ---------------------------------------------------------------

        % Untreated
        [x,~,~,~,~,~,~,omegaUT_U,~,betaU,alphaU,xBarU] = ...
            MainWithinhost(Tau,dtau,Ntau,'U','NA',...
                           T_treatmentStart,T_Recovery,1,1,0);

        % Treatment failure under A
        [~,~,~,~,~,N0TFA,omegaTF_U,~,~,betaTF,alphaTF,xBarTF] = ...
            MainWithinhost(Tau,dtau,Ntau,'TF','A',...
                           T_treatmentStart,T_Recovery,1,1,0);

        % Treatment success under A
        [~,~,~,~,~,~,omegaTS_U,~,gammaTS,betaTS,alphaTS,xBarTS] = ...
            MainWithinhost(Tau,dtau,Ntau,'TS','A',...
                           T_treatmentStart,T_Recovery,1,1,0);

        % In MainBetweenhost, omegaU_T is omegaUT from the untreated within-host run
        omegaU_T = omegaUT_U;

        % Column vectors
        omegaU_T = omegaU_T(:);
        omegaTS_U = omegaTS_U(:);
        omegaTF_U = omegaTF_U(:);
        gammaTS = gammaTS(:);

        betaTS = betaTS(:);
        betaTF = betaTF(:);
        betaU  = betaU(:);

        alphaTS = alphaTS(:);
        alphaTF = alphaTF(:);
        alphaU  = alphaU(:);

        xBarTS = xBarTS(:);
        xBarTF = xBarTF(:);
        xBarU  = xBarU(:);

        %% ---------------------------------------------------------------
        % zhat(tau) = [ITS_A, ITF_A, IU]'
        %
        % At each infection age, the three states are coupled, so solve a
        % 3 x 3 linear system.
        %% ---------------------------------------------------------------

        ITSAhat = zeros(Ntau,1);
        ITFAhat = zeros(Ntau,1);
        IUhat   = zeros(Ntau,1);

        for n = 1:Ntau

            aTS = qU*omegaTS_U(n) + alphaTS(n) + gammaTS(n) + muh;
            aTF = qU*omegaTF_U(n) + alphaTF(n) + muh;
            aU  = qA*omegaU_T(n)  + alphaU(n)  + muh;

            M = [1/dtau + aTS,           0, -qs*qA*omegaU_T(n); ...
                              0, 1/dtau + aTF, -(1-qs)*qA*omegaU_T(n); ...
                -qU*omegaTS_U(n), -qU*omegaTF_U(n), 1/dtau + aU];

            if n == 1
                % Boundary inflow corresponding to unit incidence J = 1
                rhs = [qs*qA/dtau; ...
                       (1-qs)*qA/dtau; ...
                       qU/dtau];
            else
                rhs = [ITSAhat(n-1)/dtau; ...
                       ITFAhat(n-1)/dtau; ...
                       IUhat(n-1)/dtau];
            end

            z = M\rhs;

            ITSAhat(n) = z(1);
            ITFAhat(n) = z(2);
            IUhat(n)   = z(3);
        end

        %% ---------------------------------------------------------------
        % Transmission generated by one unit of incidence
        %% ---------------------------------------------------------------

        K = trapz(Tau, ...
                  betaTS.*ITSAhat + ...
                  betaTF.*ITFAhat + ...
                  betaU .*IUhat);

%         DFE = Lambdah/muh;
        R0eq = DFE*K;

%         if ~(isfinite(R0eq) && R0eq > 1)
%             error('No positive endemic equilibrium: R0eq = %.8g <= 1.',R0eq);
%         end

        % Positive-equilibrium susceptible level
        Sstar = 1/K;

        %% ---------------------------------------------------------------
        % Recovery
        %% ---------------------------------------------------------------

        Ghat = trapz(Tau,gammaTS.*ITSAhat);

        % Stationary version of the discrete immune-age transport scheme:
        % R(eta_1) = Ghat/(1 + muh*deta)
        % R(eta_j) = R(eta_{j-1})/(1 + muh*deta)
        Rhat = zeros(1,Neta);
        survivalStep = 1/(1 + muh*deta);

        Rhat(1) = Ghat*survivalStep;

        for eta = 2:Neta
            Rhat(eta) = Rhat(eta-1)*survivalStep;
        end

        % Return-to-susceptible flux per unit incidence
        LossPerIncidence = Rhat(end);

%         if LossPerIncidence >= 1
%             error(['Invalid endemic balance: immunity-loss flux per unit incidence ' ...
%                    'is >= 1 (value %.8g).'],LossPerIncidence);
%         end

        %% ---------------------------------------------------------------
        % Equilibrium incidence from susceptible balance
        %
        % 0 = Lambda_h - muh*S* - J* + LossFlux*
        % LossFlux* = J* LossPerIncidence
        %% ---------------------------------------------------------------

        Jstar = (Lambdah - muh*Sstar)/(1 - LossPerIncidence);

%         if ~(isfinite(Jstar) && Jstar > 0)
%             error('Computed endemic incidence is non-positive: Jstar = %.8g.',Jstar);
%         end

        lambdaStar = Jstar/Sstar;

        %% ---------------------------------------------------------------
        % Scale unit-incidence profiles
        %% ---------------------------------------------------------------

        ITSAstar = Jstar*ITSAhat;
        ITFAstar = Jstar*ITFAhat;
        IUstar   = Jstar*IUhat;
        Rstar    = Jstar*Rhat;

        %% ---------------------------------------------------------------
        % Equilibrium totals and treatment-failure
        %% ---------------------------------------------------------------

        TotITS = trapz(Tau,ITSAstar);
        TotITF = trapz(Tau,ITFAstar);
        TotIU  = trapz(Tau,IUstar);
        TotI   = TotITS + TotITF + TotIU;
        TotR   = trapz(Teta,Rstar);

        PropITFstar = TotITF/(TotITS + TotITF + TotR +0.05);

        %% ---------------------------------------------------------------
        % Resistance quantities
        %% ---------------------------------------------------------------

        EtaStar = ...
            trapz(Tau, ...
                  xBarTF(end).*ITFAstar + ...
                  xBarTS(end).*ITSAstar + ...
                  xBarU(end) .*IUstar) / TotI;

        % Same structure as Eta_t:
        EtaMean = ...
            trapz(Tau, ...
                  xBarTF.*ITFAstar + ...
                  xBarTS.*ITSAstar + ...
                  xBarU .*IUstar) / TotI;

    end



end

