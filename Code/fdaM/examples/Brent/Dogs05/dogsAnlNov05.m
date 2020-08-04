%% Created: 25Nov2005
%% Modified: 25Nov2005
%%   Based on the analysis performed for the Shannon's HRV data 
%%      File - "HeartAnlRev.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do analysis on the original data using both "anlFM" and "anlTwoFM", 
%%  Provide summaries using "SummaFM" and "compTwoFM" respectively
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and possibly later on 15-min averages
%%

addpath ('c:\matlab6p5\fdaM')
addpath ('c:\matlab6p5\fdaM\examples\Brent\Dogs05')
addpath ('c:\Jarek\research\Thesis\Matlab\PHFLM')
addpath ('c:\Jarek\research\Thesis\Matlab\Mixed');
addpath ('c:\Jarek\research\Thesis\Matlab\Grandvalet');
addpath ('c:\Jarek\research\Thesis\Matlab\SimPaper');


%  -----------------------------------------------------------------------
%       Journeyman air pollution data
%  -----------------------------------------------------------------------

%  ----------------  input the data  ------------------------

% DATA (to be worked on as of 25Nov2005)
 mfmat = load ('mf.dat'); 
 mpefmat = load ('mpef.dat'); 
 pm25mat = load ('pm25.dat'); 


%% Redefining mfmat and pm25mat to be on the log scale
xMat = log(pm25mat);
% yMat = log(mfmat); %% R^2 below 20% with the best method
yMat = log(mpefmat); %% Slightly better R^2 around 30% with the best method


[K N] = size(xMat); 

T = 335;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% figure 
% plot(time,xMat)
% 
% figure 
% plot(time,yMat)

% figure 
% plot(time,y2Mat)

%  --------------- set up data
timevec = linspace(0,T,K)';
tfine = (0:5:T)'; 
nfine = length(tfine); 

 yMatlag = zeros(nfine,N); 
 xMatlag = zeros(nfine,N); 

for i=1:N 
    mftemp = interp1(timevec,yMat(:,i),tfine); 
    pm25temp = interp1(timevec,xMat(:,i),tfine); 
    yMatlag(:,i) = mftemp(1:nfine);
    xMatlag(:,i) = pm25temp(1:nfine); 
end 



%  ----------------------- Now we convert these discrete data 
%  ----------------------- to functional data objects using a B­spline basis. 

% Trying different basis functions: in the original analysis
% B-spline basis were used
nbasis = 68; 
norder = 2; 
basis = create_bspline_basis([0,T], nbasis, norder); 
mffd = data2fd(yMatlag, tfine, basis); 
pm25fd = data2fd(xMatlag, tfine, basis); 


% ----------    Does NOT work at the moment: 26 June 2003
% Trying polygonal basis (preserving raw data)
% basis = create_polygon_basis(tfine);
% mffd = data2fd(yMatlag, tfine, basis); 
% pm25fd = data2fd(xMatlag, tfine, basis); 

%  ----------------  We'll also need the mean function for each variable. 
pm25meanfd = mean(pm25fd); 
mfmeanfd = mean(mffd); 

%  ---------------  centered functions
mffd0 = center(mffd); 
pm25fd0 = center(pm25fd); 


%   ---------------------------------------------------------------
%   ------------------  Plotting the Bivariate Correlation Function 

%  --------------  First we define a fine mesh of time values. 
nfiney = 68;


nfinSt = nfiney + 1;
nfinEnd = 2*nfiney;

tmesh = linspace(0,T,nfiney)'; 

%   --------    get the discrete data from the functions. 
mffdDisc = eval_fd(mffd, tmesh); 
pm25fdDisc = eval_fd(pm25fd, tmesh); 

%   -------- plot the discretized version of the functions
subplot(2,1,1) 
plot(tmesh, pm25fdDisc) 
axis([0,T,3,8]) 
subplot(2,1,2) 
plot(tmesh, mffdDisc) 
axis([0,T,1,6])


%   --------- compute the correlations between the measures across curves. 
mfpm25corr = corrcoef([mffdDisc',pm25fdDisc']); 
mfpm25corr = mfpm25corr(nfinSt:nfinEnd,1:nfiney); 
for i=2:nfiney, mfpm25corr(i,1:i-1) = 0; end 


%  ------------- display the correlation surface
%  ----------------Use the rotation and zoom features of Matlab's 
%  -----surface display function to examine the surface from various angles. 
%  -- COLOR plotting: a) save as encapsulated color postscript file
%                     b) in GhostView use "PS printing" option
subplot(1,1,1) 
colormap(hot) 
surf(tmesh, tmesh, mfpm25corr') 
xlabel('\fontsize{16} s') 
ylabel('\fontsize{16} t') 
axis([0,T,0,T,-1,1]) 
axis('square') 
colorbar

%  ------------------ Defining the Finite Element Basis 
%  M - number of intervals to split the time range into
%  lambda - width of the time interval
M = 17; 
lambda = T/M; 

% --- Number of time blocks going back
% B = 6; 

% For non-smokers we have only 5 subjects
B = 9;

% -- got the functions "NodeIndexation" and "ParalleloGrid" from J. Ramsay
eleNodes = NodeIndexation(M, B); 

[Si, Ti] = ParalleloGrid(M, T, B); 

%   -- Estimating the Regression Function 
%  - The actual computation requires a discretization of the continuous 
%  - variable t. We define the spacing between these discrete values by 
%  - specifying the number of discrete values within each of the M intervals. 
%  - A value of two or four is usually sufficient to ensure a reasonably 
%  - accurate approximation. 

npts = 2; 
ntpts = M*npts; 
delta = lambda/(2*npts); 
tpts = linspace(delta, T - delta, M*npts)';

%  - Set up the design matrix that will be used in the discrete version 
%  - of the regression analysis. 
%%%%%   ------  This is a fairly length calculation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psiMat = DesignMatrixFD(pm25fd0, npts, M, eleNodes, Si, Ti, B); 

%  - check design matrix is NOT singular

singvals = svd(full(psiMat)); 
condition = max(singvals)/min(singvals); 
disp(['Condition number = ',num2str(condition)]) 

%  - vector of dependent variable values.  
yMat = eval_fd(mffd0,tpts)'; 
yVect = reshape(yMat, N*M*npts, 1);

%  - Least squares approximation that gives us our vector of regression 
%  - coefficients multiplying our basis functions. 

bHat = psiMat\yVect; 


%%%       ----------    Plotting the Regression Function 
subplot(1,1,1) 
colormap(hot) 
trisurf(eleNodes, Si, Ti, bHat) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t'); 
title(['mf/pm25 Data ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar

%% --->     GO TO Lasso file for comparisons: "JourneymanLasso.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% - More refined plot (color): use special function BetaEvalFD. 

svec = linspace(0,T,50*4+1)'; 
tvec = linspace(0,T,50*4+1)'; 
betaMat = BetaEvalFD(svec, tvec, bHat, M, T, lambda, ... 
eleNodes, Si, Ti, B); 

subplot(1,1,1) 
colormap(hot) 
H = imagesc(tvec, tvec, betaMat); 
xlabel('\fontsize{16} s (min)') 
ylabel('\fontsize{16} t (min)') 
axis([0,T,0,T]) 
axis('square') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 
title(['Coefficients b(s,t) ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar

%%%%    -------  Computing the Fit to the Data 
%  - Set up a large super­matrix containing the approximated mf 
%  - using the special function XMatrix. 
%%%         --- This is also a lengthy calculation. 

psiArray = XMatrix(pm25fd, tmesh, M, eleNodes, Si, Ti, B); 

%  - Matrix of approximation values for the mf acceleration curves 

yHat1 = zeros(nfiney,N); 

for i=1:N 
Xmati = squeeze(psiArray(i,:,:))'; 
yHat1(:,i) = Xmati*bHat; 
end 

%  - Approximation is based only on the estimated regression function b(s,t). 
%  - To complete the approximation, we must get the intercept function a(t). 
%  - This requires using the mean pm25 curve as a model, and subtracting 
%  - the fit that this gives from the mean mf acceleration. 

psimeanArray = XMatrix(pm25meanfd, tmesh, M, eleNodes, Si, Ti, B); 
yHatMean = squeeze(psimeanArray)'*bHat; 
mfmeanvec = eval_fd(mfmeanfd, tmesh); 
alphaHat = mfmeanvec - yHatMean; 

%  - Plot the intercept function. 
subplot(2,1,1)
plot(tmesh,alphaHat) 
axis([0,T,0,8]) 
title('Estimated intercept function'); 

%  - Plot the fit at the mean pm25 concentrations
subplot(2,1,2)
plot(tmesh,yHatMean) 
axis([0,T,-4,4]) 
title('Estimated predictions at the mean PM2.5 concentrations'); 

%  - Final fit to the data. 
yHat = alphaHat*ones(1,N) + yHat1; 

%  - Plot the data and the fit to the data. 
subplot(2,1,1) 
plot(tmesh, yHat) 
axis([0,T,2,8]) 
title('Fit to the data')
subplot(2,1,2) 
plot(tmesh, mffdDisc) 
axis([0,T,2,8])
title('Original data')

%  - Plot the residuals. 
subplot(1,1,1) 
resmat = mffdDisc - yHat; 
plot(tmesh, resmat) 
axis([0,T,-2,2]) 


%%%%%       ----- Assessing the Fit 
%  - Error sum of squares function. 
SSE = sum(resmat.^2,2); 

%  - Benchmark against which we can assess this fit, we need to get the 
%  - corresponding error sum of squares function when the model is simply 
%  - the mean mf acceleration curve. 

yMat0 = eval_fd(mffd0, tmesh); 
SSY = sum(yMat0.^2,2); 

%  - compute a squared multiple correlation function and plot it. 
%  - Don't be suprised, though, to see it go below zero; 
%  - the fit from the mean is not embedded within the fit by the model. 
RSQ = (SSY - SSE)./SSY; 
subplot(1,1,1) 
plot(tmesh,RSQ); 
axis([0,T,-0.2,1]) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%           1. Analysis with M fixed at 17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mvec = [17 17 17 17 17];
Bvec = [3 5 7 9 11];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part
[resAll, resAdd] = ...
    anlFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1);

%% elapsed_time =
%%  37.25
  
[yHatArray, SSArray] = compFM(xMat, yMat, resAll, ...
                                       time, timeLen, Mvec, Bvec);
% elapsed_time =
%   196.5310

[summa_all] = SummaFM(resAll, yHatArray, xMat, SSArray, ...
        N, Mvec, Bvec, time, timeLim, timeLen, 1, 1);

    
save('c:\matlab6p5\fdaM\examples\Brent\Dogs05\mfAnlM17.mat', ... 
    '*Array', 'res*', 'summa_all');

clear


%%%%%%%%%       Running analysis using anlTwoFM         %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mfmat = load ('mf.dat'); 
pm25mat = load ('pm25.dat'); 

%% Redefining mfmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(mfmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

Mvec = [17 17 17 17 17];
Bvec = [3 5 7 9 11];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;


[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);

% elapsed_time =
% 
%   634.9230

[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);
% elapsed_time =
% 
%   1.1937e+003

    
save('c:\JarekSim\Nan\HRVres\hrvAnlM17Two.mat','*Array', 'res*', 'summa_all');

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%     Summaries of fits %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% anlFM
load('c:\JarekSim\Nan\HRVres\hrvAnlM17.mat');
for (i=1:11),
    summa_all{i}
end;

%% summa_all{7}
% ans =
% 
%   1.0e+003 *
% 
%    -2.1845   -1.6238   -1.6238   -2.0408   -1.9653
%    -2.4060   -1.6238   -1.6238   -2.3239   -2.1461
%    -2.6114   -1.6238   -1.9620   -2.5068   -2.3468
%    -2.7729   -1.6238   -2.2691   -2.6152   -2.3908
%    -2.9002   -1.6238   -2.5121   -2.7322   -2.7465

clear

%% anlTwoFM
load('c:\JarekSim\Nan\HRVres\hrvAnlM17Two.mat');
for (i=1:11),
    summa_all{i}
end;

%% summa_all{7}
% ans =
% 
%   1.0e+004 *
% 
%       NaN +    NaNi  -0.1624            -0.1832            -0.2097             1.6043          
%       NaN +    NaNi  -0.1624            -0.1408            -0.2309             1.7195          
%       NaN +    NaNi  -0.1624            -0.1307            -0.2477             1.7930          
%       NaN +    NaNi  -0.1624            -0.0779            -0.2618             1.8401          
%       NaN +    NaNi  -0.1624            -0.1558            -0.2725             1.8742          
% 


M=17;
B=7;
T=355;
eleNodes = NodeIndexation(M, B); 
[Si, Ti] = ParalleloGrid(M, T, B);                                    
for i=1:5
figure
subplot(1,1,1) 
colormap(hot) 
trisurf(eleNodes, Si, Ti, resAll{3,1}(:,i)) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t');  
title(['HRV data ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%            2. Analysis with M fixed at 34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [34 34 34 34 34];
Bvec = [5 8 11 14 17];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1);


[yHatArray, SSArray] = compFM(xMat, yMat, resAll, ...
                                       time, timeLen, Mvec, Bvec);
% elapsed_time =
% 
%   1.2803e+003
save('c:\JarekSim\Nan\HRVres\hrvAnlM34.mat','*Array', 'res*');
clear


load('c:\JarekSim\Nan\HRVres\hrvAnlM34.mat');

[summa_all] = SummaFM(resAll, yHatArray, xMat, SSArray, ...
        N, Mvec, Bvec, time, timeLim, timeLen, 1, 1);
   
save('c:\JarekSim\Nan\HRVres\hrvAnlM34.mat','*Array', 'res*', 'summa_all');

% summa_all{7}
% 
% ans =
% 
%   1.0e+003 *
% 
%    -2.8883   -1.8216   -1.9951   -2.3392   -2.1736
%    -3.4498   -1.6238   -2.5797   -3.2028   -3.0035
%    -3.8951   -1.6238   -3.1328   -3.7178   -3.4304
%    -3.8687   -1.6238   -3.6028   -3.4497   -4.6401
%    -4.3464   -1.6238   -3.9013   -4.4409   -4.8233

clear



%%%%%%%%%%%%%%%%%%      anlTwoFM            %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

Mvec = [34 34 34 34 34];
Bvec = [5 8 11 14 17];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);
% elapsed_time =
% 
%   1.2569e+004

[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);
% elapsed_time =
% 
%   1.3716e+003

save('c:\JarekSim\Nan\HRVres\hrvAnlM34Two.mat','*Array', 'res*', 'summa_all');

clear

summa_all{7}

% ans =
% 
%   1.0e+004 *
% 
%       NaN +    NaNi  -0.1624            -0.1962            -0.2225             1.5649          
%   -0.2254            -0.1624            -0.0320            -0.3223             1.6708          
%   -0.2609            -0.1624             0.1391            -0.3795             1.7416          
%   -0.2911            -0.1624             0.1677            -0.3542             1.7939          
%   -0.3119            -0.1624             0.2013            -0.4444             1.8310          
% 


%%%%%%%%%%      Running summaries on the model with M=34        %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('c:\JarekSim\Nan\HRVres\hrvAnlM34Two.mat');
M=34;
% Bvec = [5 8 11 14 17];
B=14;
T=510;
eleNodes = NodeIndexation(M, B); 
[Si, Ti] = ParalleloGrid(M, T, B);                                    
for i=1:5
figure
subplot(1,1,1) 
colormap(hot) 
trisurf(eleNodes, Si, Ti, resAll{4,1}(:,i)) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t');  
title(['HRV data ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar
end;


mean(SSArray{3},3)
% 
% ans =
% 
%     0.1685    0.0087    0.0571    0.0622    0.0595
%     0.3159    0.0000    0.1933    0.2978    0.2485
%     0.4329    0.0000    0.3362    0.4413    0.3745
%     0.5369    0.0000    0.4505    0.5451    0.5738
%     0.5844    0.0000    0.5098    0.6299    0.6013




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  2b. Repeated analysis with M fixed at 34 and B's
%%%%%%%%%%%%  corresponding to twice the B's with M fixed at 17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [34 34 34 34 34];
Bvec = [6 10 14 18 22];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1);
% elapsed_time =
% 
%   2.1799e+003
% 

[yHatArray, SSArray] = compFM(xMat, yMat, resAll, ...
                                       time, timeLen, Mvec, Bvec);
% elapsed_time =
%   1.4123e+003

save('c:\JarekSim\Nan\HRVres\hrvAnlM34rep.mat','*Array', 'res*');

clear


load('c:\JarekSim\Nan\HRVres\hrvAnlM34rep.mat');

[summa_all] = SummaFM(resAll, yHatArray, xMat, SSArray, ...
        N, Mvec, Bvec, time, timeLim, timeLen, 1, 1);
% elapsed_time =
%   263.5490
   
save('c:\JarekSim\Nan\HRVres\hrvAnlM34rep.mat','*Array', 'res*', 'summa_all');

% round(summa_all{7})
% ans =
% 
%        -3090       -1624       -2140       -2465       -2787
%        -3776       -1624       -2908       -3565       -3544
%        -3869       -1624       -3603       -3450       -4640
%        -4668       -1624       -4020       -4546       -4693
%        -2145       -1624       -4160       -4646       -4802

clear










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%            3. Analysis with M fixed at 51
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);


Mvec = [51 51 51 51];
Bvec = [10 15 20 25];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part
[resAll, resAdd] = ...
    anlFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1);
% elapsed_time =
%   1.0551e+004

save('c:\JarekSim\Nan\HRVres\hrvAnlM51.mat', 'res*');

load('c:\JarekSim\Nan\HRVres\hrvAnlM51.mat');

[yHatArray, SSArray] = compFM(xMat, yMat, resAll, ...
                                       time, timeLen, Mvec, Bvec);
% elapsed_time =
%   702.2180

save('c:\JarekSim\Nan\HRVres\hrvAnlM51.mat','*Array', 'res*');


load('c:\JarekSim\Nan\HRVres\hrvAnlM51.mat');
[summa_all] = SummaFM(resAll, yHatArray, xMat, SSArray, ...
        N, Mvec, Bvec, time, timeLim, timeLen, 1, 1);
% elapsed_time =
% 
%   408.3770
   
save('c:\JarekSim\Nan\HRVres\hrvAnlM51.mat','*Array', 'res*', 'summa_all');

summa_all{7}
% ans =
%   1.0e+005 *
% 
%    -0.0413   -0.0189   -0.0303   -0.0384   -0.0391
%    -0.0458   -0.0162   -0.0441   -0.0401   -0.0579
%    -0.0572   -0.0162   -0.0523   -0.0550   -0.0669
%     1.2010   -0.0162   -0.0524   -0.0611   -0.0693

load('c:\JarekSim\Nan\HRVres\hrvAnlM51.mat');
M=51;
B=20;
T=510;
eleNodes = NodeIndexation(M, B); 
[Si, Ti] = ParalleloGrid(M, T, B); 
for i=1:5
figure
subplot(1,1,1) 
colormap(hot) 
trisurf(eleNodes, Si, Ti, resAll{3,1}(:,i,1)) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t');  
title(['HRV data ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar
end;

clear





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%              END of the analysis as of 11Mar2005 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   3c. Rerun of the Analysis with M fixed at 51 and
%%%%%%%%%%%%   corresponding B's to the run with M=17 and half-steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [51 51 51];
Bvec = [12 18 24];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part
[resAll, resAdd] = ...
    anlFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1);
% elapsed_time =
%   8.4140e+003

  [yHatArray, SSArray] = compFM(xMat, yMat, resAll, ...
                                       time, timeLen, Mvec, Bvec);
%  elapsed_time =
%   537.3600

[summa_all] = SummaFM(resAll, yHatArray, xMat, SSArray, ...
        N, Mvec, Bvec, time, timeLim, timeLen, 1, 1);
% elapsed_time =
%   127.3440

   
save('c:\JarekSim\Nan\HRVres\hrvAnlM51repHalf.mat','*Array', 'res*', 'summa_all');

summa_all{7}
ans =
  1.0e+003 *
   -4.2815   -1.6238   -3.5305   -3.7998   -5.0458
   -5.7148   -1.6238   -5.2131   -5.2688   -6.4577
    5.3793   -1.6238   -5.2250   -6.0063   -6.8732



load('c:\JarekSim\Nan\HRVres\hrvAnlM51.mat');
M=51;
B=20;
T=510;
eleNodes = NodeIndexation(M, B); 
[Si, Ti] = ParalleloGrid(M, T, B); 
for i=1:5
figure
subplot(1,1,1) 
colormap(hot) 
trisurf(eleNodes, Si, Ti, resAll{3,1}(:,i,1)) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t');  
title(['HRV data ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar
end;

clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [51 51 51 51 51];
Bvec = [9 15 21 27 33];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part
[resAll, resAdd] = ...
    anlFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1);
% elapsed_time =
%   1.0551e+004

save('c:\JarekSim\Nan\HRVres\hrvAnlM51rep.mat', 'res*');

load('c:\JarekSim\Nan\HRVres\hrvAnlM51rep.mat');

[yHatArray, SSArray] = compFM(xMat, yMat, resAll, ...
                                       time, timeLen, Mvec, Bvec);
% elapsed_time =
%  

save('c:\JarekSim\Nan\HRVres\hrvAnlM51rep.mat','*Array', 'res*');


load('c:\JarekSim\Nan\HRVres\hrvAnlM51rep.mat');
[summa_all] = SummaFM(resAll, yHatArray, xMat, SSArray, ...
        N, Mvec, Bvec, time, timeLim, timeLen, 1, 1);
% elapsed_time =
% 
   
save('c:\JarekSim\Nan\HRVres\hrvAnlM51rep.mat','*Array', 'res*', 'summa_all');

summa_all{7}




%% M = 51
%%%%%%%%%%%%%%%%%%      anlTwoFM            %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [51 51 51 51];
Bvec = [9 15 21 27];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);
%% RUNNING for 16 hours and still at B=21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);

save('c:\JarekSim\Nan\HRVres\hrvAnlM51Tworep.mat','*Array', 'res*', 'summa_all');

clear







% M=51;
% B=20;
% T=510;
% eleNodes = NodeIndexation(M, B); 
% [Si, Ti] = ParalleloGrid(M, T, B); 
% for i=1:5
% figure
% subplot(1,1,1) 
% colormap(hot) 
% trisurf(eleNodes, Si, Ti, resAll{3,1}(:,i,1)) 
% xlabel('\fontsize{12} s'); 
% ylabel('\fontsize{12} t');  
% title(['HRV data ','M = ',num2str(M),', B = ',num2str(B)]) 
% colorbar
% end;
% 
% plot(squeeze(SSArrayHrlMs{3,1}(4,3,:)))
% 
% 
% clear
                                   



%%%%%%%%%%%%%%%%%%      anlTwoFM            %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [34 34 34 34];
Bvec = [6 10 14 18];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);
% elapsed_time =
%   7.9244e+003

[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);
% elapsed_time =
%   686.7810


save('c:\JarekSim\Nan\HRVres\hrvAnlM34Tworep.mat','*Array', 'res*', 'summa_all');

clear

summa_all{7}
% round(summa_all{7}(2:4,:))
% ans =
%        -2470       -1624         907       -3593       17211
%        -2911       -1624        1677       -3542       17939
%        -3206       -1624        2119       -4558       18411
% round(summa_all{7}(:,2:4))
% ans =
%        -1624       -2062       -2264
%        -1624         907       -3593
%        -1624        1677       -3542
%        -1624        2119       -4558

%% ROUNDED 
%NaN +    NaNi       -1624       -2062       -2264
%        -2470       -1624         907       -3593       17211
%        -2911       -1624        1677       -3542       17939
%        -3206       -1624        2119       -4558       18411

%%%%%%%%%%      Running summaries on the model with M=34        %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('c:\JarekSim\Nan\HRVres\hrvAnlM34Tworep.mat');
% M=34;
% % Bvec = [6 10 14 18 22];
% B=14;
% T=510;
% eleNodes = NodeIndexation(M, B); 
% [Si, Ti] = ParalleloGrid(M, T, B);                                    
% for i=1:5
% figure
% subplot(1,1,1) 
% colormap(hot) 
% trisurf(eleNodes, Si, Ti, resAll{4,1}(:,i)) 
% xlabel('\fontsize{12} s'); 
% ylabel('\fontsize{12} t');  
% title(['HRV data ','M = ',num2str(M),', B = ',num2str(B)]) 
% colorbar
% end;
% 
% 
% mean(SSArray{3},3)
% 


%%%%%%%%%       Rerun on 22Mar2005          %%%%%%%%%%%%%%%%%%%%%%%%%%
%%      run in smaller chunks (a. Bvec=[9 12 15])
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% M = 51
%%%%%%%%%%%%%%%%%%      anlTwoFM            %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [51 51 51];
Bvec = [9 12 15];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);
% elapsed_time =
%   3.8005e+004

[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);
% elapsed_time =
%   426.0470
save('c:\JarekSim\Nan\HRVres\hrvAnlM51Tworep15.mat','*Array', 'res*', 'summa_all');

% round(summa_all{7})
% ans =
%        -2338       -1624        -734       -3756       16056
%        -2707       -1624        2355       -3716       16715
%        -3156       -1624        -233       -4094       17218
       
       
clear


%%%%%%%%%       Bvec = [18 21]              %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [51 51];
Bvec = [18 21];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);
% elapsed_time =
%   5.8720e+004

[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);
% elapsed_time =
%   450.0470

save('c:\JarekSim\Nan\HRVres\hrvAnlM51Tworep21.mat','*Array', 'res*', 'summa_all');

clear

% round(summa_all{7})
% ans =
%        -3622       -1624        3010       -5270       17613
%        -3654       -1624        4189       -5760       17946
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      RERUN in the evening of 28Mar2005           %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% M = 51
%%%%%%%%%%%%%%%%%%      anlTwoFM            %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

Mvec = [51 51];
Bvec = [23 24];
% Bvec = [24 27];
%% Did NOT converge at Bvec=27, run analysis instead at Bvec = [23 24]

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);
% elapsed_time =
%   6.6733e+004

[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);
% elapsed_time =
%   545.2820
  
save('c:\JarekSim\Nan\HRVres\hrvAnlM51Two24.mat','*Array', 'res*', 'summa_all');

% round(summa_all{7})
% ans =
%        -3647       -1624        3532       -5939       18135
%        -3628       -1624        3397       -6119       18212

clear


%% TO BE RUN at night on 29Mar2005    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% M = 51
%%%%%%%%%%%%%%%%%%      anlTwoFM            %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

Mvec = [51];
Bvec = [25];
% Bvec = [24 27];
%% Did NOT converge at Bvec=27, run analysis instead at Bvec = [23 24]

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);
% elapsed_time =
%   3.7384e+004
[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);
% elapsed_time =
%   286.5780

round(summa_all{7})
% ans =
%        -3638       -1624        3503       -6156       18283

save('c:\JarekSim\Nan\HRVres\hrvAnlM51Two26.mat','*Array', 'res*', 'summa_all');

clear


%% TO BE RUN at night on 30Mar2005    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% M = 51
%%%%%%%%%%%%%%%%%%      anlTwoFM            %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

Mvec = [51];
Bvec = [26];
% Bvec = [24 27];
%% Did NOT converge at Bvec=26 & 27

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part: anlFM
[resAll, resAdd] = ...
    anlTwoFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1, [1 30 20 5]);

[yHatArray, xHatArray, SSArray, summa_all] = ...
    compTwoFM(xMat, yMat, resAll, time, timeLen, Mvec, Bvec, 2);

round(summa_all{7})
% ans =
%        -3638       -1624        3503       -6156       18283

save('c:\JarekSim\Nan\HRVres\hrvAnlM51Two26Plus.mat','*Array', 'res*', 'summa_all');

clear

%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%       %%%%%%%%%%

hrvmat = load ('hrvJour.dat'); 
pm25mat = load ('pm25Jour.dat'); 

%% Redefining hrvmat and pm25mat to be on the log scale
xMat = log(pm25mat);
yMat = log(hrvmat);

[K N] = size(xMat); 

T = 510;
timeLim = [0,T];
time = linspace(0,T,K)';
timeLen = size(time,1);

% Mvec = [17 17 17 17 17];
% Bvec = [3 5 7 9 11];

Mvec = [51];
Bvec = [26];

finElem = size(Mvec,2);

%% Setting up the output data sets
psiCell = cell(finElem,1);
yCell = cell(finElem,1);

for g=1:finElem;
    M=Mvec(g);
    B=Bvec(g);
    totalPts = N*M*2;
    
    [psiMat,yVecfd] = ...
            DatPrepFM(xMat, yMat, M, B, time, timeLen);

    psiCell{g,1} = psiMat;
    yCell{g,1} = yVecfd;
end;

%% Analysis part
[resAll, resAdd] = ...
    anlFM(psiCell, yCell, N, Mvec, Bvec, timeLim, timeLen, 1, 1);


[yHatArray, SSArray] = compFM(xMat, yMat, resAll, ...
                                       time, timeLen, Mvec, Bvec);
[summa_all] = SummaFM(resAll, yHatArray, xMat, SSArray, ...
        N, Mvec, Bvec, time, timeLim, timeLen, 1, 1);
   
save('c:\JarekSim\Nan\HRVres\hrvAnlM51rep26plus.mat','*Array', 'res*', 'summa_all');
indAnl =
     1
elapsed_time =
  4.7584e+003
j =
     1
elapsed_time =
  253.9690
j =
     1
elapsed_time =
   67.5790


round(summa_all{7})
% ans =
%        69371       -1624       -5271       -6151       -7013



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Plot of the fit with B=26       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%        lowest hAIC = -7013         %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('c:\JarekSim\Nan\HRVres\hrvAnlM51rep26plus.mat');

M=51;
B=26;
T=510;
eleNodes = NodeIndexation(M, B); 
[Si, Ti] = ParalleloGrid(M, T, B); 
for i=1:5
figure
subplot(1,1,1) 
colormap(hot) 
trisurf(eleNodes, Si, Ti, resAll{1,1}(:,i,1)) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t');  
title(['HRV data ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar
end;
saveas(H, 'C:\Jarek\research\Thesis\Matlab\SimPaper\Plots\ManL1M25B5.eps', 'psc2');
saveas(H, 'C:\Jarek\research\Thesis\Matlab\SimPaper\Plots\ManL1M25B5.jpg');


%%%%%%%%        PLOT for the manuscript     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

svec = linspace(0,T,205)'; 
tvec = linspace(0,T,205)'; 
lambda = T/M; 
betaMat = BetaEvalFD(svec, tvec, resAll{1,1}(:,5,1), M, T, lambda, ... 
eleNodes, Si, Ti, B); 

figure
subplot(1,1,1) 
colormap(hot) 
H = imagesc(tvec, tvec, betaMat); 
xlabel('\fontsize{16} s (min)') 
ylabel('\fontsize{16} t (min)') 
axis([0,T,0,T]) 
axis('square') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 
hold on
plot([45,45], [0,510], ...
     [70,70], [0,510], ...
     [215, 215], [0,510], ...
     [250, 250], [0,510], ...
     [405, 405], [0,510], ...
     [430, 430], [0,510])
hold off
title(['Estimated regression surface \beta(s,t) ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar





















% hrvmat    - response matrix (pefmat)
% pm25mat   - predictor matrix (chamat)

%  ---------------  define number of records and time values in minutes
N = size(pm25mat,2); 
K = size(pm25mat,1);


% OLD data: Starting time is 8:20 and end time is 16:50
T = 510;

timevec = linspace(0,T,K)'; 

%%%%%%%%%%%     IMPORTANT       %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrvmat = hrvmatLog;
pm25mat = pm25matLog;


%  --------------- set up data
tfine = (0:1:T)'; 
nfine = length(tfine); 

 hrvmatlag = zeros(nfine,N); 
 pm25matlag = zeros(nfine,N); 

for i=1:N 
    hrvtemp = interp1(timevec,hrvmat(:,i),tfine); 
    pm25temp = interp1(timevec,pm25mat(:,i),tfine); 
    hrvmatlag(:,i) = hrvtemp(1:nfine);
    pm25matlag(:,i) = pm25temp(1:nfine); 
end 



%  ----------------------- Now we convert these discrete data 
%  ----------------------- to functional data objects using a B­spline basis. 

% Trying different basis functions: in the original analysis
% B-spline basis were used
%nbasis = 103; 
%norder = 2; 
nbasis = 34;
norder = 4;
basis = create_bspline_basis([0,T], nbasis, norder); 
hrvfd = data2fd(hrvmatlag, tfine, basis); 
pm25fd = data2fd(pm25matlag, tfine, basis); 


% ----------    Does NOT work at the moment: 26 June 2003
% Trying polygonal basis (preserving raw data)
% basis = create_polygon_basis(tfine);
% hrvfd = data2fd(hrvmatlag, tfine, basis); 
% pm25fd = data2fd(pm25matlag, tfine, basis); 

%  ----------------  We'll also need the mean function for each variable. 
pm25meanfd = mean(pm25fd); 
hrvmeanfd = mean(hrvfd); 

%  ---------------  centered functions
hrvfd0 = center(hrvfd); 
pm25fd0 = center(pm25fd); 


%   ---------------------------------------------------------------
%   ------------------  Plotting the Bivariate Correlation Function 

%  --------------  First we define a fine mesh of time values. 
% Work
nfiney = 103;

% Leisure
% nfiney = 112; 

nfinSt = nfiney + 1;
nfinEnd = 2*nfiney;

tmesh = linspace(0,T,nfiney)'; 

%   --------    get the discrete data from the functions. 
hrvfdDisc = eval_fd(hrvfd, tmesh); 
pm25fdDisc = eval_fd(pm25fd, tmesh); 

%   -------- plot the discretized version of the functions
subplot(2,1,1) 
plot(tmesh, pm25fdDisc) 
axis([0,T,-5,3]) 
subplot(2,1,2) 
plot(tmesh, hrvfdDisc) 
axis([0,T,0,6])


%   --------- compute the correlations between the measures across curves. 
hrvpm25corr = corrcoef([hrvfdDisc',pm25fdDisc']); 
hrvpm25corr = hrvpm25corr(nfinSt:nfinEnd,1:nfiney); 
for i=2:nfiney, hrvpm25corr(i,1:i-1) = 0; end 

%% JUST for plotting set the variable "tplot"
tplot = tmesh/60 + 8 + 1/3;

%  ------------- display the correlation surface
%  ----------------Use the rotation and zoom features of Matlab's 
%  -----surface display function to examine the surface from various angles. 
%  -- COLOR plotting: a) save as encapsulated color postscript file
%                     b) in GhostView use "PS printing" option
subplot(1,1,1) 
colormap(hot) 
%surf(tmesh, tmesh, hrvpm25corr') 
H=imagesc(tplot, tplot, hrvpm25corr'); 
%surf(tmesh, tmesh, hrvpm25corr') 
xlabel('\fontsize{16} s') 
ylabel('\fontsize{16} t') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 
%axis([0,T,0,T,-1,1]) 
axis([8.3,16.85,8.3,16.85]) 
axis('square') 
title('Correlations between HRV and PM2.5')
colorbar

%  ------------------ Defining the Finite Element Basis 
%  M - number of intervals to split the time range into
%  lambda - width of the time interval
M = 10; 
% M = 103/4 (Leisure = 224/4) 
lambda = T/M; 

% --- Number of time blocks going back
B = 5; 

% -- got the functions "NodeIndexation" and "ParalleloGrid" from J. Ramsay
eleNodes = NodeIndexation(M, B); 

[Si, Ti] = ParalleloGrid(M, T, B); 

%   -- Estimating the Regression Function 
%  - The actual computation requires a discretization of the continuous 
%  - variable t. We define the spacing between these discrete values by 
%  - specifying the number of discrete values within each of the M intervals. 
%  - A value of two or four is usually sufficient to ensure a reasonably 
%  - accurate approximation. 

npts = 2; 
ntpts = M*npts; 
delta = lambda/(2*npts); 
tpts = linspace(delta, T - delta, M*npts)';

%  - Set up the design matrix that will be used in the discrete version 
%  - of the regression analysis. 
%%%%%   ------  This is a fairly length calculation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psiMat = DesignMatrixFD(pm25fd0, npts, M, eleNodes, Si, Ti, B); 

%  - check design matrix is NOT singular

singvals = svd(full(psiMat)); 
condition = max(singvals)/min(singvals); 
disp(['Condition number = ',num2str(condition)]) 

%  - vector of dependent variable values.  
yMat = eval_fd(hrvfd0,tpts)'; 
yVect = reshape(yMat, N*M*npts, 1);

%  - Least squares approximation that gives us our vector of regression 
%  - coefficients multiplying our basis functions. 

bHat = psiMat\yVect; 


%%%       ----------    Plotting the Regression Function 
subplot(1,1,1) 
colormap(hot) 
trisurf(eleNodes, Si, Ti, bHat) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t');   
title(['hrv/pm25 Data ','M = ',num2str(M),', B = ',num2str(B)]) 

%%%%        -----       Triangle plot for Enar
elEnar = size(Si);
elEnar = elEnar(2);
trEnar = diag(zeros(elEnar));
Siplot = Si/60 + 8 + 1/3;
Tiplot = Ti/60 + 8 + 1/3;

subplot(1,1,1) 
%colormap(hot) 
trisurf(eleNodes, Siplot, Tiplot, trEnar) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t'); 
title(['Region of integration',' M = ',num2str(M),', B = ',num2str(B)]) 
axis([8.3,16.85,8.3,16.85]) 
axis('square') 


% - More refined plot (color): use special function BetaEvalFD. 

svec = linspace(0,T,50*4+1)'; 
tvec = linspace(0,T,50*4+1)'; 
betaMat = BetaEvalFD(svec, tvec, bHat, M, T, lambda, ... 
eleNodes, Si, Ti, B); 

subplot(1,1,1) 
colormap(hot) 
H = imagesc(tvec, tvec, betaMat); 
xlabel('\fontsize{16} s (min)') 
ylabel('\fontsize{16} t (min)') 
axis([0,T,0,T]) 
axis('square') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 
hold on
plot([45,45], [0,510], ...
     [70,70], [0,510], ...
     [215, 215], [0,510], ...
     [250, 250], [0,510], ...
     [405, 405], [0,510], ...
     [430, 430], [0,510])
hold off
title(['Coefficients b(s,t) ','M = ',num2str(M),', B = ',num2str(B)]) 
colorbar

%%%%    -------  Computing the Fit to the Data 
%  - Set up a large super­matrix containing the approximated hrv 
%  - using the special function XMatrix. 
%%%         --- This is also a lengthy calculation. 

psiArray = XMatrix(pm25fd, tmesh, M, eleNodes, Si, Ti, B); 

%  - Matrix of approximation values for the hrv acceleration curves 

yHat1 = zeros(nfiney,N); 

for i=1:N 
Xmati = squeeze(psiArray(i,:,:))'; 
yHat1(:,i) = Xmati*bHat; 
end 

%  - Approximation is based only on the estimated regression function b(s,t). 
%  - To complete the approximation, we must get the intercept function a(t). 
%  - This requires using the mean pm25 curve as a model, and subtracting 
%  - the fit that this gives from the mean hrv acceleration. 

psimeanArray = XMatrix(pm25meanfd, tmesh, M, eleNodes, Si, Ti, B); 
yHatMean = squeeze(psimeanArray)'*bHat; 
hrvmeanvec = eval_fd(hrvmeanfd, tmesh); 
alphaHat = hrvmeanvec - yHatMean; 

%  - Plot the intercept function. 
subplot(2,1,1)
plot(tmesh,alphaHat) 
axis([0,T,0,8]) 
title('Estimated intercept function'); 

%  - Plot the fit at the mean pm25 concentrations
subplot(2,1,2)
plot(tmesh,yHatMean) 
axis([0,T,-4,4]) 
title('Estimated predictions at the mean PM2.5 concentrations'); 

%  - Final fit to the data. 
yHat = alphaHat*ones(1,N) + yHat1; 

%  - Plot the data and the fit to the data. 
subplot(2,1,1) 
plot(tmesh, yHat) 
axis([0,T,0,6]) 
title('Fit to the data')
subplot(2,1,2) 
plot(tmesh, hrvfdDisc) 
axis([0,T,0,6])
title('Original data')

%  - Plot the residuals. 
subplot(1,1,1) 
resmat = hrvfdDisc - yHat; 
plot(tmesh, resmat) 
axis([0,T,-2,2]) 


%%%%%       ----- Assessing the Fit 
%  - Error sum of squares function. 
SSE = sum(resmat.^2,2); 

%  - Benchmark against which we can assess this fit, we need to get the 
%  - corresponding error sum of squares function when the model is simply 
%  - the mean hrv acceleration curve. 

hrvmat0 = eval_fd(hrvfd0, tmesh); 
SSY = sum(hrvmat0.^2,2); 

%  - compute a squared multiple correlation function and plot it. 
%  - Don't be suprised, though, to see it go below zero; 
%  - the fit from the mean is not embedded within the fit by the model. 
RSQ = (SSY - SSE)./SSY; 
subplot(1,1,1) 
plot(tmesh,RSQ); 
axis([0,T,-0.2,1]) 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%       MIXED model approximation using function mixedFDA.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Points from the orals paper (Chapter 7: Computational issues)


% 5. Create functional object of time on an interval [0,0.64]

% a. creation of Z matrix
hrtno = 50;
timebasis = linspace(0, 510, hrtno);
timefd = bsplineFDA(tpts, timebasis, 2);

Ztime = full(timefd);
Ztime = Ztime(:,1:(hrtno-1));
ZtimeAll = repmat(Ztime,14,1);


XhMat = DesignMatrixFD(pm25fd0, npts, M, eleNodes, Si, Ti, B); 

ZmixAll = XhMat;

size_tpts = size(tpts);
size_tpts = size_tpts(1);

% b. creation of X matrix
Xmix = diag(eye(size_tpts));
XmixAll = repmat(Xmix,14,1);

% c. creation of id vector
id = linspace(1,size_tpts,size_tpts);
idMat = repmat(id,14,1);
idVect = reshape(idMat, 14*size_tpts, 1);

% d. dimensions of the random effects
dimRand = size(Si);
dimRand = dimRand(2);
s20start = [0.002, 0.25];


% 6. Set up "y" vector
%  - vector of dependent variable values.  
yMat = eval_fd(hrvfd0,tpts)'; 
yVect = reshape(yMat, N*M*npts, 1);

%  - noncentered response
%yMat = eval_fd(lipfd,tpts)'; 
%yVect = reshape(yMat, N*M*npts, 1);

% 7. First try at solving the Henderson's equations

[s2,b,u,Is2,C,H,q,loglik,loops] = mixedFDA(yVect,XmixAll,ZmixAll,idVect,dimRand,s20start,dimRand,1);

%% Fixed variance components
[s2,b,u,Is2,C,H,q,loglik,loops] = mixedFDA(yVect,XmixAll,ZmixAll,idVect,dimRand,s20start,dimRand,0);


norm(bHat - u(1:dimRand))

norm(bHat - diag(zeros(dimRand)))

norm(u - diag(zeros(dimRand)))

plot(bHat,u(1:dimRand), 'o')

subplot(1,1,1)
plot([1:dimRand], u, '-r', 'LineWidth',1.5)
hold on
plot([1:dimRand], bHat, '-.b')
hold off
title('PHFLM (red) vs. HFLM (blue)')
 
% - More refined plot (color): use special function BetaEvalFD. 

svec = linspace(0,T,50*4+1)'; 
tvec = linspace(0,T,50*4+1)'; 
uMat = BetaEvalFD(svec, tvec, u, M, T, lambda, ... 
eleNodes, Si, Ti, B); 

subplot(1,1,1) 
colormap(hsv) 
H = imagesc(tvec, tvec, uMat); 
xlabel('\fontsize{16} s (min)') 
ylabel('\fontsize{16} t (min)') 
axis([0,T,0,T]) 
axis('square') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 
title(['Random coefficients "u" for ',' M = ',num2str(M),', B = ',num2str(B)]) 
colorbar


%% LUnch and coffee breaks
hold on
plot([45,45], [0,510], ...
     [70,70], [0,510], ...
     [215, 215], [0,510], ...
     [250, 250], [0,510], ...
     [405, 405], [0,510], ...
     [430, 430], [0,510])
hold off
 
 
%%%%%%%%%%%%%%%     USE mixed models solutions      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%    -------  Computing the Fit to the Data 
%  - Set up a large super­matrix containing the approximated hrv 
%  - using the special function XMatrix. 
%%%         --- This is also a lengthy calculation. 

%%%%%%%%%        DONE for bHat already
%psiArray = XMatrix(pm25fd, tmesh, M, eleNodes, Si, Ti, B); 

%  - Matrix of approximation values for the hrv acceleration curves 

yHat1u = zeros(nfiney,N); 

for i=1:N 
Xmati = squeeze(psiArray(i,:,:))'; 
yHat1u(:,i) = Xmati*u; 
end 

%  - Approximation is based only on the estimated regression function b(s,t). 
%  - To complete the approximation, we must get the intercept function a(t). 
%  - This requires using the mean pm25 curve as a model, and subtracting 
%  - the fit that this gives from the mean hrv acceleration. 

%%%%%%%%%        DONE for bHat already
%psimeanArray = XMatrix(pm25meanfd, tmesh, M, eleNodes, Si, Ti, B); 

yHatuMean = squeeze(psimeanArray)'*u; 
hrvmeanvec = eval_fd(hrvmeanfd, tmesh); 
alphaHatu = hrvmeanvec - yHatuMean; 

%  - Plot the intercept function. 
subplot(2,1,1)
plot(tmesh,alphaHatu) 
axis([0,T,0,8]) 
title('Estimated intercept function'); 

%  - Plot the fit at the mean pm25 concentrations
subplot(2,1,2)
plot(tmesh,yHatuMean) 
axis([0,T,-4,4]) 
title('Estimated predictions at the mean PM2.5 concentrations'); 

%  - Final fit to the data. 
yHatu = alphaHatu*ones(1,N) + yHat1u; 

%  - Plot the data and the fit to the data. 
subplot(2,1,1) 
plot(tmesh, yHatu) 
axis([0,T,0,6]) 
title('Fit to the data')
subplot(2,1,2) 
plot(tmesh, hrvfdDisc) 
axis([0,T,0,6])
title('Original data')

%  - Plot the residuals. 
subplot(1,1,1) 
resmatu = hrvfdDisc - yHatu; 
plot(tmesh, resmatu) 
axis([0,T,-2,2]) 


%%%%%       ----- Assessing the Fit 
%  - Error sum of squares function. 
SSEu = sum(resmatu.^2,2); 

%  - Benchmark against which we can assess this fit, we need to get the 
%  - corresponding error sum of squares function when the model is simply 
%  - the mean hrv acceleration curve. 

hrvmat0 = eval_fd(hrvfd0, tmesh); 
SSY = sum(hrvmat0.^2,2); 

%  - compute a squared multiple correlation function and plot it. 
%  - Don't be suprised, though, to see it go below zero; 
%  - the fit from the mean is not embedded within the fit by the model. 
RSQu = (SSY - SSEu)./SSY; 
subplot(1,1,1) 
plot(tmesh,RSQu); 
hold on;
  plot(tmesh,RSQ, '.-');
hold off;
axis([0,T,-0.2,1]) 

mean(RSQu)
mean(RSQ)


 
 
 
 
 
subplot(2,1,1) 
colormap(hot) 
H = imagesc(tvec, tvec, betaMat); 
xlabel('\fontsize{16} s (min)') 
ylabel('\fontsize{16} t (min)') 
axis([0,T,0,T]) 
axis('square') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 
%trisurf(eleNodes, Si, Ti, bHat) 
%xlabel('\fontsize{12} s'); 
%ylabel('\fontsize{12} t'); 
title(['HRV vs. PM2.5 - HFLM ','M = ',num2str(M),', B = ',num2str(B)]) 

colormap(hot) 
trisurf(eleNodes, Si, Ti, u) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t'); 
title(['HRV vs. PM2.5 - PHFLM ','M = ',num2str(M),', B = ',num2str(B)]) 




