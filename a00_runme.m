
%% Replication files for figures in
% The Impact of Pessimistic Expectations on the Effects of Covid-19-
% induced Uncertainty in the Euro Area
% by Giovanni Pellegrino, Federico Ravenna & Gabriel Züllig
% Corresponding author for codes: Gabriel Züllig, gabrielzuellig@snb.ch.
% Oxford Bulletin of Economics and Statistics, 83(4), 2021


clear all; clc;

addpath('./model_functions/') 
addpath('./support_functions/')  % heavily relying on utils by A. Cesa Bianchi


%% Data preparation
uctyvars = {'VSTOXX','VSTOXXres','VIX','JLNFin','ConsDisp','IndDisp'}; % 
statevars = {'ConsConfEC','IndConfEC','EconSent','ConsConf','ConsConfECres'};
% These vectors describe all 'shocked' and all 'state' variables that are
% interacted with one another.
c01_prepare_data;  % load data, conducts necessary transformations


%% Generall settings across specifications
% VAR specification
pick = 1;         % variable to be shocked (Cholesky ordering)
c_case = 1;       % 1 for constant, 2 for constant and trend in our VAR 
nlags = 4;        % number of lags in VAR (max. 12)
ordend = pick;    % position of the shocked variable that is part of the interaction term
ordend2 = 3;      % position of the 'state' variable that enters the interaction term second
ordend2_lin = 0;  % for the linear model, there is no interaction term
defstate.growth = 0; % 1 if state should enter in growth rates, e.g. change of confidence
defstate.cumul = 0;
proxyvar = 0;     % can identify shock via exogenous instruments, see example in robustness
extexoshock = 0;  % index of variables that should be considered exogenous, 
                  % i.e. not load on anything else but its own lags
% computational settings 
% settings to 0 will save considerable time, but might result in
% failure
fast = 1;         % good approximation of GIRFs
comp = 1;
dolinear = 1;     % run linear model
dobootstrap = 1;  % run bootstrap of GIRFs, should be set to 1
dofevd = 0;       % compute factor error variance decomposition
% IRFs
shocksize = 0;    % 0 = 1 standard deviation shock, 1 = unit shock, else = size in standard deviation terms
nsteps = 36;      % number of months to generate GIRFs
mode = 0;         % 0 = draw from the empirical distribution of residuals for future shocks, 1 = from Gaussian distribution respecting VCV matrix
method = 0;       % 0 = draws from the empirical distribution, 1 = Monte Carlo, 2 = fixed design wild bootstrap based on a Rademacher distribution
option = 0;       % 0, 1 or 2. Use 0 for fully empirical papers     
reverseHL = 0;    % reverses high and low state (e.g. if compare states of outputgap and unemployment)
horplus = 24;     % need IRFs not to explode for nsteps+horplus periods
% bootstrap
ndraws = 500;     % number of draws
confLev = 68; 
confLev2 = 90;
irfmedian = 0;    % 0 = mean response, 1 = median 
irfmeanfromdraws = 0;



%% Baseline model
% fast = 0;       % The figure in the paper is produced using ndraws=10k and fast = 0.
draws = 10000;
perc = 20;
nlags = 4;  
pick = 1;
ordend = 1;
ordend2 = 4;
spec = strcat('m01_baseline');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','ConsConfEC','EONIA','VSTOXXxConsConfEC'};

r01_subset_data;  % Generates Figure 1

r02_ivar; 


%% Figures 2 + 3
load data/data_m_full.mat
load m01_baseline/out.mat

labels_print={'Uncertainty','Industrial production','Inflation','Consumer confidence','Policy rate'}
leg=1;  % turns off mean response and makes legend more legible
plPE=0;

figure()
VARirplot_2New_CCP(IRFs,OIRF2avg,OIRFavg,0,pick,labels_print,INFavgH2,SUPavgH2,INFavgL2,SUPavgL2)
tightfig;
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2','-loose','paper_results/fig2.eps');

figure()
VARtestplot_2New_CCP(dOIRFavg,dOIRFavg,pick,labels_print,dIRFinf,dIRFsup,dIRFinf2,dIRFsup2,leg,plPE)
tightfig;
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2','-loose','paper_results/fig3.eps');


%% Run all robustness checks
dolinear = 0;
draws = 500;
c02_robustness;   % Generates Figures 4-6




