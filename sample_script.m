%%%%%%%%%%%%%%%%%%%%%%% WPHFLM Sample Script %%%%%%%%%%%%%%%%%%%%%%%
% Simulated Peak Effect , burnin 1000, sample 1000                 %
% x(v) ~ GP(0, S) S~AR(1) estimated covariance from data           %
% N = 20, T = 2^5                                                  %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NB: This script is meant as a guide. We use a simulated dataset, but  %%
%%      simX and W can be replaced with any matrix of the desired        %%
%%      functional X and scalar effects W. Y is the matrix of functional %%
%%      outcomes.                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add paths %%
% addpath('./Code');
addpath('/Users/mjm556/Dropbox/Research/Drafts/Historical/Code')

%% set number of grid points %%
T       = 2^5;
N       = 20;

%% grid for simulation surfaces %%
sDens   = 1; % 1 (64), 0.5 (128), 0.25 (256), 0.125 (512), 0.0625 (1024)
[v, t]  = meshgrid(0:sDens:(T-sDens));

%% scale to control STNR %%
stnrs   = 4006/100;

%% bivariate normal params %%
mv      = 19;
sv      = 2.5;
mt      = 7.5;
st      = 2.5;
rho     = 0;

%% redefine T %%
T           = size(t,1);

%% Simulate Coefficients %%
b1      = zeros(T,T);

for i = 1:T
    ti  = t(i,1);
    for j = 1:T
        vj          = v(1,j);
        b1(i,j)     = stnrs*(1/(2*pi*sv*st*sqrt(1-rho^2)))*exp((-1/(2*(1-rho^2)))*(((vj-mv)^2)/(sv^2) + ((ti-mt)^2)/(st^2) - (2*rho*(vj-mv)*(ti-mt))/(sv*st)));
    end
end

%% constrain true surface %%
bh  = zeros(size(b1));
for i = 1:size(b1,1)
    for j = i:size(b1,2)
        bh(i,j) = b1(i,j);
    end
end

%% generate ar(1) covariance pattern %%
ar1Corr     = eye(T);
sigma       = 3.5; % get from Journeyman data
rho         = 0.75; % get from Journeyman data
for i = 1:T
    for j = (i+1):T
        ar1Corr(i,j) = rho^(j-i);
        ar1Corr(j,i) = rho^(j-i);
    end
end

ar1Cov      = sigma*ar1Corr;

%% generate x(v) based on AR(1) estimated cov %%
simX        = NaN(N,T);
muX         = zeros(T,1);
for i = 1:N
    simX(i,:)   = mvnrnd(muX,ar1Cov);
end

%% Generate Covariance structure for E ~ GP(0, Sigma_E)
CovStr      = eye(T);
sigmae      = 0.1;
rho         = 0.5;
for i = 1:T
    for j = (i+1):T
        CovStr(i,j) = rho^(j-i);
        CovStr(j,i) = rho^(j-i);
    end
end

Sigma_E     = sigmae*CovStr;

%% specs for packets %%
wpspecs.wavelet     = 'db3';
wpspecs.wtmode      = 'zpd';
wpspecs.nlevels     = 3;

%% decompose X(v) %%
[ Dx, wpspecsx ]    = dwpt_rows(simX,wpspecs); %% replace simX
wpspecs.wpspecsx    = wpspecsx;

%% set lag based on number of levels %%
perlagback          = 1.1;    % 1.1 for full surface
wpspecs.lag         = round(perlagback*size(Dx,2)/(2^(wpspecs.nlevels)));
wpspecs.perlagback  = perlagback;

%% update Dx based on threshold %%
model.Tx            = size(Dx,2);
model.thresh        = model.Tx*(6/8); % model.Tx*(6/8);
model.keep          = 1:(model.Tx-model.thresh);
Dx                  = Dx(:,model.keep);

%% MCMCspecs %%
MCMCspecs.B                 = 1000;
MCMCspecs.burnin            = 1000;
MCMCspecs.thin              = 1;
MCMCspecs.propsdTheta       = 1;
MCMCspecs.nj_nosmooth       = 1;        % can be 0, if is 0, don't do any constraint for pi_{ij}, bigTau_{ij}.
MCMCspecs.minp              = 1e-14;
MCMCspecs.maxO              = 1e20;
MCMCspecs.minVC             = 1e-20;
MCMCspecs.VC0_thresh        = 1e-4;
MCMCspecs.delta_theta       = 1e-4;
MCMCspecs.thetaMLE_maxiter  = 10^6;
MCMCspecs.EmpBayes_maxiter  = 10^6;
MCMCspecs.time_update       = 100;
MCMCspecs.tau_prior_var     = 1e3;      % the variance of tau_{ijk} when finding prior parameters for tau_{ijk}.
MCMCspecs.tau_prior_idx     = 1;        % 1 indicate that a_tau and b_tau depend on ij, 0 indicate that they depend on jk. 
MCMCspecs.PI_prior_var      = 0.06;     % this range should be in [0.02 0.09].

%% Use Sigma_E to generate matrix of model errors %%
E           = zeros(N,T);
muE         = zeros(T,1);
for i = 1:N
    E(i,:)  = mvnrnd(muE,Sigma_E);
end

%% construct and decompose Y(t) %%
Y                   = simX*bh + E;
[ Dy, wpspecsy ]    = dwpt_rows(Y,wpspecs); %% replace Y
wpspecs.wpspecsy    = wpspecsy;
    
%% set up function specs %%
a   = ones(N,1);
W   = []; % add scalar covariates here
Z   = []; % use fixed effects only
model.X     = [ Dx a W ];
model.Z{1}  = Z;
model.W     = W;
model.Dx    = Dx;
model.C     = ones(size(Y,1),1);
model.H     = 1;
model.Hstar = 0;

%% run hwfmm %%
tic;
res             = wphflm(Dy,model,wpspecs,MCMCspecs);
res.MCMCrun     = toc;

%% extract results for post-processing %%
MCMC_beta           = res.MCMC_beta;
MCMC_zeta           = res.MCMC_zeta;
MCMC_alpha          = res.MCMC_alpha;
MCMC_flag_theta     = res.MCMC_flag_theta;
MCMC_tau            = res.MCMC_tau;
MCMC_pi             = res.MCMC_pi;
MCMC_theta          = res.MCMC_theta;
theta               = res.theta;
model               = res.model;
wpspecs             = res.wpspecs;
model.delt          = 0.05;
model.alf           = 0.05;

%% run post-processor %%
postout             = PostProcess(MCMC_beta,MCMC_zeta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wpspecs);
postout.runtime     = toc;
postout.res         = res;
    
%% view posterior results %%

%% set-up for graphing %%
T       = size(postout.bhat,1);
bmath   = nan(T,T);
for i = 1:T
    for j = i:T
        if j < i+round(T*1.1)
            bmath(i,j) = postout.bhat(i,j);
        end
    end
end

%% define bfFlag %%
%  indexes hist coefs  %
bf      = zeros(size(postout.bhat));
for i = 1:size(bf,1)
    for j = i:size(bf,2)
        bf(i,j) = 1;
    end
end
T       = size(bf,1);
bfFlag	= reshape(bf,1,T*T);

%% generate joint intervals and SimBa Scores %%
%  NB: uses alpha = 0.05, change if desired   %
[SBSb, upper_CIb, lower_CIb]	= jointband_sbs(postout.bINF,0.05);
SBS                             = zeros(size(bfFlag));
SBS(bfFlag == 1)                = SBSb;
sbs1                            = reshape(SBS,T,T);

USBS                            = zeros(size(bfFlag));
USBS(bfFlag == 1)               = upper_CIb;
usbs                            = reshape(USBS,T,T);

LSBS                            = zeros(size(bfFlag));
LSBS(bfFlag == 1)               = lower_CIb;
lsbs                            = reshape(LSBS,T,T);

%% set-up for graphing %%
sbs     = nan(T,T);
sbsf    = nan(T,T);
sbsu    = nan(T,T);
sbsl    = nan(T,T);
for i = 1:T
    for j = i:T
        if j < i+round(T*1.1)
            sbs(i,j)    = sbs1(i,j);
            sbsu(i,j)   = usbs(i,j);
            sbsl(i,j)   = lsbs(i,j);
            sbsf(i,j)   = 1*(sbs(i,j) < 0.05);
        end
    end
end
    
%% estimated surface
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(bmath','FaceColor','interp','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14);
ylabel('t','FontSize',14);
title({'\beta(v,t)'});

%% simba scores
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(sbs','FaceColor','interp','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
colormap(copper)
caxis([0 0.5])
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14)
ylabel('t','FontSize',14);
title({'SimBa Scores'});

%% significant coef
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(sbsf','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
colormap(copper)
caxis([0 1])
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14);
ylabel('t','FontSize',14);
title({'Significant Coefficients'});

%% joint lower
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(sbsl','FaceColor','interp','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14);
ylabel('t','FontSize',14);
title({'Joint Interval, Lower Bound'});

%% joint upper
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(sbsu','FaceColor','interp','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14);
ylabel('t','FontSize',14);
title({'Joint Interval, Upper Bound'});


