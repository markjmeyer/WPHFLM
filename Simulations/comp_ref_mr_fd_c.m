%%%%%%%%%%%%%%%%%%%%%%%% Cumulative Comparison Sim %%%%%%%%%%%%%%%%%%%%%%%%
% x(v) ~ GP(0, S) S~AR(1) estimated covariance from data                  %
% N = 50/200, T = 2^6/2^7                                                 %
%                                                                         %
% Created:      02/26/2020                                                %
% Modified:     07/08/2020                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% add paths %%
addpath('/Users/mjm556/Dropbox/Research/Matlab/fdaM')
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1)

%% set number of grid points %%
T       = 2^7; % 2^6
N       = 200; % 50
ndata   = 200;

%% grid for simulation surfaces %%
sDens   = 1; % 1 (64), 0.5 (128), 0.25 (256), 0.125 (512), 0.0625 (1024)
[v, t]  = meshgrid(0:sDens:(T-sDens));

%% scale to control STNR %%
% lagged %
b1var   = 0.01;
b1  = (1/(60*sqrt(2*pi*b1var)))*exp(-1/(2*b1var)*(t./T-v./T-0.5).^2);
b1  = b1';
% cumulative %
b2  = b1;
for i = 1:size(b1,2)
%     b2(i,:) = 2*(1-(size(b1,2)-i)/T)*b1(i,:);
    b2(i,:) = 2*(1-(T-i)/T)*b1(i,:);
end

%%
b3      = (1/max(max(b2)))*b2;

%% coefficients for evaluation %%
bh      = zeros(size(b3));
bf      = zeros(size(b3));
for i = 1:size(b3,1)
    for j = i:size(b3,2)
        bh(i,j) = b3(i,j);
        bf(i,j) = 1;
    end
end
T       = size(b3,1);
bhr     = reshape(bh,1,T*T);
bFlag     = reshape(bf,1,T*T);

%% subset historical coefs %%
bht     = bhr(bFlag == 1);

%% set number of grid points for coarse comparison %%
Tc        = 14;

%% grid for simulation surfaces %%
sDensc   = 9; % 4.5
[vc, tc]  = meshgrid(0:sDensc:(T-sDensc));

%% scale to control STNR %%
b1c  = (1/(60*sqrt(2*pi*b1var)))*exp(-1/(2*b1var)*(tc./Tc-vc./Tc-0.5).^2);
b1c  = b1c';
% cumulative %
b2c  = b1c;
for i = 1:size(b1c,2)
%     b2(i,:) = 2*(1-(size(b1,2)-i)/T)*b1(i,:);
    b2c(i,:) = 2*(1-(Tc-i)/Tc)*b1c(i,:);
end

%% coefficients for evaluation %%
bhc      = zeros(size(b2c));
bfc      = zeros(size(b2c));
for i = 1:size(b2c,1)
    for j = i:size(b2c,2)
        bhc(i,j) = b2c(i,j);
        bfc(i,j) = 1;
    end
end
bhrc     = reshape(bhc,1,Tc*Tc);
bFlagc     = reshape(bfc,1,Tc*Tc);

%% subset historical coefs %%
bhtc     = (1/max(bhrc))*bhrc(bFlagc == 1);

%% set number of grid points for Refund comparison %%
Tr        = 40;

%% grid for simulation surfaces %%
sDensr   = 3.2; % 1.6
[vr, tr]  = meshgrid(0:sDensr:(T-sDensr));

%% scale to control STNR %%
b1r  = (1/(60*sqrt(2*pi*b1var)))*exp(-1/(2*b1var)*(tr./Tr-vr./Tr-0.5).^2);
b1r  = b1r';
% cumulative %
b2r  = b1r;
for i = 1:size(b1r,2)
%     b2(i,:) = 2*(1-(size(b1,2)-i)/T)*b1(i,:);
    b2r(i,:) = 2*(1-(Tr-i)/Tr)*b1r(i,:);
end

%% coefficients for evaluation %%
bhr      = zeros(size(b2r));
bfr      = zeros(size(b2r));
for i = 1:size(b2r,1)
    for j = i:size(b2r,2)
        bhr(i,j) = b2r(i,j);
        bfr(i,j) = 1;
    end
end
bhrr     = reshape(bhr,1,Tr*Tr);
bFlagr     = reshape(bfr,1,Tr*Tr);

%% subset historical coefs %%
bhtr     = (1/max(bhrr))*bhrr(bFlagr == 1);

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

%% Generate Covariance structure for E ~ GP(0, Sigma_E) %%
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

%% empty matrices %%
miseS50     = NaN(ndata,3);
pwcoS50     = NaN(ndata,1);

%% MR model set up %%
ts          = (T-1)/100;
timevec     = linspace(0,ts,T);
tfine       = (0:0.0025:ts)';
nfine       = length(tfine);

%% simulate 200 datasets %%
tic
for seed = 1:ndata
    %% set seed %%
    rng(seed)

    %% Use Sigma_E to generate matrix of model errors %%
    E           = zeros(N,T);
    muE         = zeros(T,1);
    for i = 1:N
        E(i,:)  = mvnrnd(muE,Sigma_E);
    end

    %% construct and decompose Y(t) %%
    Y           = simX*bh + E;
        
    %% convert to func data objs %%
    ymatlag   = zeros(nfine, N);
    xmatlag   = zeros(nfine, N);

    for i=1:N
        ytemp = interp1(timevec,Y(i,:),tfine);
        xtemp = interp1(timevec,simX(i,:),tfine);
        ymatlag(:,i) = ytemp;
        xmatlag(:,i) = xtemp;
    end

    tfine = tfine(1:nfine);

    nbasis = 93;
    norder =  6;
    basis  = create_bspline_basis([0,ts], nbasis, norder);

    yfd = data2fd(ymatlag, tfine, basis);
    xfd = data2fd(xmatlag, tfine, basis);

    yfd0 = center(yfd);
    xfd0 = center(xfd);

    %% define finite element basis %%
%     M = 63;
    M = 13; %39
    lambda = ts/M;    
    B = M;
    eleNodes = NodeIndexation(M, B);
    [Si, Ti] = ParalleloGrid(M, ts, B);

    %% estimating the regression function %%
    npts  = 4; %4
    ntpts = M*npts;
    delta = lambda/(2*npts);
    tpts  = linspace(delta, ts-delta, M*npts)';

    psiMat  = DesignMatrixFD(xfd0, npts, M, eleNodes, Si, Ti, B);

    singvals  = svd(full(psiMat));
    condition = max(singvals)/min(singvals);
%     disp(['Condition number = ',num2str(condition)])

    yMat  = eval_fd(yfd0,tpts)';
    yVect = reshape(yMat, N*M*npts, 1);

    bHatc = psiMat\yVect;
    
    %% mise for M&R %%
    miseS50(seed,1)	= mean((bHatc' - bhtc).^2);
    
    %% call R, run Refund %%
    cd('/Users/mjm556/Documents/MATLAB')
    save('Y.mat', 'Y');
    save('simX.mat', 'simX');
    
    %% call R from matlab %%
    ! R CMD BATCH refFDB.R outfile.txt
    
    %% load Refund and FDBoost output %%
    bHref   = load('bhmatR.mat');
    bHRef   = bHref.bhmat;
    bHfdb   = load('bhmatF.mat');
    bHFDB   = bHfdb.bhmatF;
    lowr    = load('lowerR.mat');
    lowmat  = lowr.lower;
    uppr    = load('upperR.mat');
    uppmat  = uppr.upper;
    
    %% process Refund output %%
    bhre      = zeros(size(bHRef));
    bhfd      = zeros(size(bHFDB));
    lowh      = zeros(size(bHRef));
    upph      = zeros(size(bHRef));
    for i = 1:size(bHRef,1)
        for j = i:size(bHRef,2)
            bhre(i,j)   = bHRef(i,j);
            bhfd(i,j)   = bHFDB(i,j);
            lowh(i,j)   = lowmat(i,j);
            upph(i,j)   = uppmat(i,j);
        end
    end
    Trf      = size(bHRef,1);
    bhrf     = reshape(bhre,1,Trf*Trf);
    Tfd      = size(bHFDB,1);
    bhfb     = reshape(bhfd,1,Tfd*Tfd);
    lorf     = reshape(lowh,1,Trf*Trf);
    uprf     = reshape(upph,1,Trf*Trf);
    bHatr    = bhrf(bFlagr == 1);
    bHatf    = bhfb(bFlagr == 1);
    lower    = lorf(bFlagr == 1);
    upper    = uprf(bFlagr == 1);
    
    %% mise for Refund %%
    miseS50(seed,2)	= mean((bHatr - bhtr).^2);
    miseS50(seed,3)	= mean((bHatf - bhtr).^2);

    %% coverage for refund %%
    pwcoS50(seed)   = mean(upper > bhtr & lower < bhtr);    
    
    %%
    if mod(seed, 5) == 0
        fprintf('.')
    end
    if mod(seed, 50) == 0
       fprintf('\n %d \n',seed),toc;
    end    
    
end
% toc

%% MISE evaluation %%
mean(sqrt(miseS50))
mean(pwcoS50)

%% save output %%
% cd('/Users/mjm556/Dropbox/Research/Drafts/Historical/Comp Sim')
% save('mise50cT64.mat', 'miseS50')
% save('pwco50cT64.mat', 'pwcoS50')

%%
% cd('/Users/mjm556/Dropbox/Research/Drafts/Historical/Comp Sim')
% save('mise50cT128.mat', 'miseS50')
% save('pwco50cT128.mat', 'pwcoS50')

%%
% cd('/Users/mjm556/Dropbox/Research/Drafts/Historical/Comp Sim')
% save('mise200cT64.mat', 'miseS50')
% save('pwco200cT64.mat', 'pwcoS50')

%%
% cd('/Users/mjm556/Dropbox/Research/Drafts/Historical/Comp Sim')
% save('mise200cT128.mat', 'miseS50')
% save('pwco200cT128.mat', 'pwcoS50')

%%
