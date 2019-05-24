%%%%%%%%%%%%%%%%%%%%%%%%%%% Odysssey SIM %%%%%%%%%%%%%%%%%%%%%%%%%%%
% HWFM Simulation based on Generated data for Odyssey              %
% Depressed Effect , burnin 1000, sample 1000                      %
% x(v) ~ GP(0, S) S~AR(1) estimated covariance from data           %
% sample tau and pi                                                %
% N = 50, T = 2^6                                                  %
%                                                                  %
% Created:      03/20/2014                                         %
% Modified:     04/04/2019                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add paths %%
% addpath('./CODE');
addpath('/Users/mjm556/Dropbox/Research/Drafts/Historical/WPHFLM Code')
% seed    = 1;

%% get seed %%
% seeds   = getenv('SLURM_ARRAY_TASK_ID');
% seed    = str2double(seeds);

%% set seed %%
% rng(seed);

%% set number of grid points %%
T       = 2^6;
N       = 50;

%% grid for simulation surfaces %%
sDens   = 0.5; % 1 (64), 0.5 (128), 0.25 (256), 0.125 (512), 0.0625 (1024)
[v, t]  = meshgrid(0:sDens:(T-sDens));

%% scale to control STNR %%
stnrs   = 1257/4; %5000/557;

%% bivariate normal params %%
mv      = 45;
sv      = 10;
mt      = 10;
st      = 5;
rho     = 0;

%% redefine T %%
T           = size(t,1);

%% Simulate Coefficients %%
b1      = zeros(T,T);

for i = 1:T
    ti  = t(i,1);
    for j = 1:T
        vj          = v(1,j);
        b1(i,j)     = -stnrs*(1/(2*pi*sv*st*sqrt(1-rho^2)))*exp((-1/(2*(1-rho^2)))*(((vj-mv)^2)/(sv^2) + ((ti-mt)^2)/(st^2) - (2*rho*(vj-mv)*(ti-mt))/(sv*st)));
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
sigma       = 3.5; %1.5252;   % get from Journeyman data
rho         = 0.75; %0.7697;   % get from Journeyman data
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

%% shell specs for hacketts
wpspecs.wavelet     = 'db3';
wpspecs.wtmode      = 'zpd';
wpspecs.nlevels     = 3;

%% decompose X(v) %%
[ Dx, wpspecsx ]    = dwpt_rows(simX,wpspecs);
wpspecs.wpspecsx    = wpspecsx;

%% set lag based on number of levels %%
perlagback          = 1.1;    % 1.1 for full surface
wpspecs.lag         = round(perlagback*size(Dx,2)/(2^(wpspecs.nlevels)));
wpspecs.perlagback  = perlagback;

%% perform hard thresholding on Dx %%
cDx2    = cumsum(Dx.^2,2);
pcDx2   = zeros(size(cDx2));
for i = 1:size(pcDx2,1)
    pcDx2(i,:) = cDx2(i,:)./sum(Dx(i,:).^2);
end
% vDx                 = var(Dx);
% nu                  = size(vDx,2);
% unth                = sqrt(2*log(nu)); % universal threshold
% thresh              = find(vDx > unth, 1, 'last' );

%% update Dx based on threshold %%
% model.Tx            = size(Dx,2);
% model.thresh        = model.Tx-thresh;
% model.keep          = 1:(model.Tx-model.thresh);
% Dx                  = Dx(:,model.keep);
model.Tx            = size(Dx,2);
model.thresh        = model.Tx*(4/8); % model.Tx*(6/8);
model.keep          = 1:(model.Tx-model.thresh);
Dx                  = Dx(:,model.keep);

%% MCMCspecs
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
% if isempty(Z)
%     MCMCspecs.sampleU           = 0;
% else
%     MCMCspecs.sampleU           = 1;
% end

% sampleU, get_sigma
% sampleU     = MCMCspecs.sampleU;
% get_sigma   = sampleU;

%% set seed %%
rng(2017) % 2020
count   = 1; %0

%% simulate 200 datasets %%
while count < 200
    %% set seed %%
%     rng(seed);
    
    %% Use Sigma_E to generate matrix of model errors
    E           = zeros(N,T);
    muE         = zeros(T,1);
    for i = 1:N
        E(i,:)  = mvnrnd(muE,Sigma_E);
    end

    %% construct and decompose Y(t) %%
    Y                   = simX*bh + E;
    [ Dy, wpspecsy ]    = dwpt_rows(Y,wpspecs);
    wpspecs.wpspecsy    = wpspecsy;
    
    %% set up function specs %%
    a   = ones(N,1);
    W   = []; % scalar covariates
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
    try
        res             = wphflm(Dy,model,wpspecs,MCMCspecs);
    catch ME
        continue
    end
    res.MCMCrun     = toc;
    count           = count + 1;
    
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

    %% save output %%
    fname               = sprintf('/Users/mjm556/Documents/Research/Historical/Simulation/Depressed/Sparse/Half/N50/n50t64hd%d.mat',count);
    save(fname,'postout');

    %% clear large output %%
    clear res postout
end
    
%% end




