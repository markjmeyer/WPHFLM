addpath ('c:\matlab6p5\fdaM')
addpath ('c:\matlab6p5\fdaM\examples\lip')

% JH (10Mar03): changed paths above to Nick's computer
%  Last modified 10 March 2003


%  -----------------------------------------------------------------------
%       Lip Movement Data (extended analysis)
%	Modified from annotated analysis from Applied FDA
%			
%  -----------------------------------------------------------------------

%  ----------------  input the data  ------------------------

EMGmat = load ('EMG.dat'); 
lipmat = load ('LipAcc.dat'); 
lipposmat = load ('LipPos.dat'); 

% lipmat    - response
% EMGmat    - predictor

%  ---------------  define number of records and time values in seconds
N = size(EMGmat,2); 
timevec = linspace(0,0.69,501)'; 

%  --------------- set up data with EMG lagged by 50 milliseconds
tfine = (0:0.001:0.69)'; 
nfine = length(tfine); 
% nfine50 = nfine - 50; % JH change to include nfine50 in zeros below
lipmatlag = zeros(nfine-50,N); 
EMGmatlag = zeros(nfine-50,N); 

for i=1:N 
    liptemp = interp1(timevec,lipmat(:,i),tfine); 
    EMGtemp = interp1(timevec,EMGmat(:,i),tfine); 
    lipmatlag(:,i) = liptemp(51:nfine); 
    EMGmatlag(:,i) = EMGtemp(1:(nfine-50)); 
end 


nfine = nfine - 50; 
tfine = tfine(1:nfine); 
T = 0.64; 

%  ----------------------- Now we convert these discrete data 
%  ----------------------- to functional data objects using a B­spline basis. 
nbasis = 93; 
norder = 6; 
basis = create_bspline_basis([0,T], nbasis, norder); 
lipfd = data2fd(lipmatlag, tfine, basis); 
emgfd = data2fd(EMGmatlag, tfine, basis); 

%  ----------------  We'll also need the mean function for each variable. 
emgmeanfd = mean(emgfd); 
lipmeanfd = mean(lipfd); 

%  ---------------  centered functions
lipfd0 = center(lipfd); 
emgfd0 = center(emgfd); 


%   ------------------  Plotting the Bivariate Correlation Function 

%  --------------  First we define a fine mesh of time values. 
nfiney = 101; 
tmesh = linspace(0,T,nfiney)'; 

%   --------    get the discrete data from the lagged functions. 
lipmat = eval_fd(lipfd, tmesh); 
emgmat = eval_fd(emgfd, tmesh); 

%   --------- compute the correlations between the measures across curves. 
lipemgcorr = corrcoef([lipmat',emgmat']); 
lipemgcorr = lipemgcorr(102:202,1:101); 
for i=2:101, lipemgcorr(i,1:i-1) = 0; end 

% JH add: 30Apr2003
subplot(2,1,1)
plot(tmesh,lipmat) 
title('Lip acceleration'); 

subplot(2,1,2)
plot(tmesh,emgmat) 
title('EMG activity'); 
% JH: end addition


%  ------------- display the correlation surface
%  ----------------Use the rotation and zoom features of Matlab's 
%  -----surface display function to examine the surface from various angles. 
%  -- COLOR plotting: a) save as encapsulated color postscript file
%                     b) in GhostView use "PS printing" option
subplot(1,1,1) 
colormap(hot) 
surf(tmesh, tmesh, lipemgcorr') 
xlabel('\fontsize{16} s') 
ylabel('\fontsize{16} t') 
axis([0,.64,0,.64,-1,1]) 
axis('square') 

%  ------------------ Defining the Finite Element Basis 
M = 13; 
lambda = T/M; 

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

npts = 4; 
ntpts = M*npts; 
delta = lambda/(2*npts); 
tpts = linspace(delta, T - delta, M*npts)';

%  - Set up the design matrix that will be used in the discrete version 
%  - of the regression analysis. 
%%%%%   ------  This is a fairly lengthy calculation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psiMat = DesignMatrixFD(emgfd0, npts, M, eleNodes, Si, Ti, B); 

%  - check design matrix is NOT singular

singvals = svd(full(psiMat)); 
condition = max(singvals)/min(singvals); 
disp(['Condition number = ',num2str(condition)]) 

%  - vector of dependent variable values.  
yMat = eval_fd(lipfd0,tpts)'; 
yVect = reshape(yMat, N*M*npts, 1);

%  - Least squares approximation that gives us our vector of regression 
%  - coefficients multiplying our basis functions. 

bHat = psiMat\yVect; 


%%%       ----------    Plotting the Regression Function 

trisurf(eleNodes, Si, Ti, bHat) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t'); 
title(['Lip/EMG Data ','M = ',num2str(M),', B = ',num2str(B)]) 


% - More refined plot (color): use special function BetaEvalFD. 

svec = linspace(0,T,13*6+1)'; 
tvec = linspace(0,T,13*6+1)'; 
betaMat = BetaEvalFD(svec, tvec, bHat, M, T, lambda, ... 
eleNodes, Si, Ti, B); 
subplot(1,1,1) 
colormap(hot) 
H = imagesc(tvec*1000, tvec*1000, betaMat); 
xlabel('\fontsize{16} s (msec)') 
ylabel('\fontsize{16} t (msec)') 
axis([0,640,0,640]) 
axis('square') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 

%% JH: addition 30Apr2003: Plot function \hat{\beta}(t,t)
betatt = diag(betaMat);

subplot(1,1,1)
plot(tvec,betatt);
ylabel('\fontsize{14} Estimated Beta')
xlabel('\fontsize{14} seconds (s=t)')
%% JH: end addition (very different plot from the paper)


%%%%    -------  Computing the Fit to the Data 
%  - Set up a large super­matrix containing the approximated lip 
%  - acceleration curves using the special function XMatrix. 
%%%         --- This is also a lengthy calculation. 

psiArray = XMatrix(emgfd, tmesh, M, eleNodes, Si, Ti, B); 

%  - Matrix of approximation values for the lip acceleration curves 

yHat1 = zeros(nfiney,N); 

for i=1:N 
Xmati = squeeze(psiArray(i,:,:))'; 
yHat1(:,i) = Xmati*bHat; 
end 

%  - Approximation is based only on the estimated regression function b(s,t). 
%  - To complete the approximation, we must get the intercept function a(t). 
%  - This requires using the mean EMG curve as a model, and subtracting 
%  - the fit that this gives from the mean lip acceleration. 

psimeanArray = XMatrix(emgmeanfd, tmesh, M, eleNodes, Si, Ti, B); 
yHatMean = squeeze(psimeanArray)'*bHat; 
lipmeanvec = eval_fd(lipmeanfd, tmesh); 
alphaHat = lipmeanvec - yHatMean; 

%  - Plot the intercept function. 
plot(tmesh,alphaHat) 
title('Estimated intercept function'); 

%  - Final fit to the data. 
yHat = alphaHat*ones(1,N) + yHat1; 

%  - Plot the data and the fit to the data. 
subplot(2,1,1) 
plot(tmesh, yHat) 
axis([0,T,-4,4])
title('Fit to the data');
subplot(2,1,2) 
plot(tmesh, lipmat) 
axis([0,T,-4,4])
title('Original data');

%  - Plot the residuals. 
subplot(1,1,1) 
% JH: used to be strange operation "--" changed to "-"
resmat = lipmat - yHat; 
plot(tmesh, resmat) 
axis([0,T,-4,4]) 


%%%%%       ----- Assessing the Fit 
%  - Error sum of squares function. 
SSE = sum(resmat.^2,2); 

%  - Benchmark against which we can assess this fit, we need to get the 
%  - corresponding error sum of squares function when the model is simply 
%  - the mean lip acceleration curve. 

lipmat0 = eval_fd(lipfd0, tmesh); 
SSY = sum(lipmat0.^2,2); 

%  - compute a squared multiple correlation function and plot it. 
%  - Don't be suprised, though, to see it go below zero; 
%  - the fit from the mean is not embedded within the fit by the model. 
RSQ = (SSY - SSE)./SSY; 
subplot(1,1,1) 
plot(tmesh,RSQ); 
axis([0,0.64,-0.5,1]) 
title('RSQ');
