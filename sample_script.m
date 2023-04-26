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
% addpath('~/Code');
addpath('~/Documents/Code/WPHFLM/Code')

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

%% Use Sigma_E to generate matrix of model errors %%
E           = zeros(N,T);
muE         = zeros(T,1);
for i = 1:N
    E(i,:)  = mvnrnd(muE,Sigma_E);
end

%% construct and decompose Y(t) %%
Y                   = simX*bh + E;

%% model specs %%
model.alf           = 0.05; % global alpha
model.delt          = 0.05; % threshold for Bayesian FDR
model.wpKeep        = 0.25; % keep wpKeep percent of x-space wavelet packets

%% specs for packets %%
wpspecs.wavelet     = 'db3';
wpspecs.wtmode      = 'zpd';
wpspecs.nlevels     = 3;

%% MCMCspecs %%
MCMCspecs.B                 = 1000;
MCMCspecs.burnin            = 1000;
MCMCspecs.thin              = 1;

%% run model %%
postout = wphflm(Y, simX, model, wpspecs, MCMCspecs);

%% estimated surface
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(postout.bmath','FaceColor','interp','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90);
colorbar(axes1,'Ticks',[0 0.5 1 1.5 2],'FontSize',14);
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14);
ylabel('t','FontSize',14);
title({'Meyer, Malloy, & Coull (2021)'});

%% simba scores
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(postout.sbs','FaceColor','interp','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
colormap(copper)
caxis([0 0.5])
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14)
ylabel('t','FontSize',14);
title({'SimBa Scores'});

%% significant coef
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(postout.sbsf','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
colormap(copper)
caxis([0 1])
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14);
ylabel('t','FontSize',14);
title({'Significant Coefficients'});

%% joint lower
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(postout.sbsl','FaceColor','interp','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14);
ylabel('t','FontSize',14);
title({'Joint Interval, Lower Bound'});

%% joint upper
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');
surf(postout.sbsu','FaceColor','interp','EdgeColor','none'); xlim([1 T]); ylim([1 T]); view(0, 90); colorbar
set(axes1,'FontSize',14,'XGrid','on','YGrid','on');
xlabel('v','FontSize',14);
ylabel('t','FontSize',14);
title({'Joint Interval, Upper Bound'});


