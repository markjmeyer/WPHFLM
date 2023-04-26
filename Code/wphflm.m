function postout = wphflm(Y, X, model, wpspecs, MCMCspecs)

%% decompose X(v) %%
[ Dx, wpspecsx ]    = dwpt_rows(X,wpspecs); %% replace simX
wpspecs.wpspecsx    = wpspecsx;

%% set lag based on number of levels %%
perlagback          = 1.1;    % 1.1 for full surface
wpspecs.lag         = round(perlagback*size(Dx,2)/(2^(wpspecs.nlevels)));
wpspecs.perlagback  = perlagback;

%% update Dx based on threshold %%
model.Tx            = size(Dx,2);
model.thresh        = model.Tx*(1-model.wpKeep); % model.Tx*(6/8);
model.keep          = 1:(model.Tx-model.thresh);
Dx                  = Dx(:,model.keep);

%% MCMCspecs %%
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

%% construct and decompose Y(t) %%
[ Dy, wpspecsy ]    = dwpt_rows(Y,wpspecs); %% replace Y
wpspecs.wpspecsy    = wpspecsy;
    
%% set up function specs %%
N   = size(Y,1);
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
res             = wphflm_fit(Dy,model,wpspecs,MCMCspecs);
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
postout.bmath       = bmath;

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
postout.bfFlag       = bfFlag;

%% generate joint intervals and SimBa Scores %%
%  NB: uses alpha = 0.05, change if desired   %
[SBSb, upper_CIb, lower_CIb]	= jointband_sbs(postout.bINF,model.alf);
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
postout.sbs       = sbs;
postout.sbsu      = sbsu;
postout.sbsl      = sbsl;
postout.sbsf      = sbsf;



