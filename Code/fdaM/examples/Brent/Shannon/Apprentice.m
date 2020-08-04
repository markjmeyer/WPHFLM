% Analysis to be performed on Shannon Magari's data (based on "dogs.m")
% Prediction of "hrv" - Heart Rate Variability 
% from "pm25.epa" - pm 2.5 level
% -------------------------------------------------------------------
% Created: 20Jun03
% Modified: 26Jun03

addpath ('c:\matlab6p5\fdaM')
addpath ('c:\matlab6p5\fdaM\examples\Brent\Shannon')

%  -----------------------------------------------------------------------
%       Apprentice air pollution data
%  -----------------------------------------------------------------------

%  ----------------  input the data  ------------------------

hrvmat = load ('hrvAppr.dat'); 
pm25mat = load ('pm25Appr.dat'); 

% hrvmat    - response matrix (pefmat)
% pm25mat   - predictor matrix (chamat)

%  ---------------  define number of records and time values in days
N = size(pm25mat,2); 

% Starting time is 9:15 and end time is 13:55
timevec = linspace(0,275,56)'; 

% -- plotting raw data
subplot(1,1,1) 
plot(timevec, hrvmat) 
axis([0,T,0,280])

%  --------------- set up data
tfine = (0:1:275)'; 
nfine = length(tfine); 

 hrvmatlag = zeros(nfine,N); 
 pm25matlag = zeros(nfine,N); 

for i=1:N 
    hrvtemp = interp1(timevec,hrvmat(:,i),tfine); 
    pm25temp = interp1(timevec,pm25mat(:,i),tfine); 
    hrvmatlag(:,i) = hrvtemp(1:nfine);
    pm25matlag(:,i) = pm25temp(1:nfine); 
end 


T = 275; 

%  ----------------------- Now we convert these discrete data 
%  ----------------------- to functional data objects using a B­spline basis. 
nbasis = 20; 
norder = 4; 
basis = create_bspline_basis([0,T], nbasis, norder); 
hrvfd = data2fd(hrvmatlag, tfine, basis); 
pm25fd = data2fd(pm25matlag, tfine, basis); 

%  ----------------  We'll also need the mean function for each variable. 
pm25meanfd = mean(pm25fd); 
hrvmeanfd = mean(hrvfd); 

%  ---------------  centered functions
hrvfd0 = center(hrvfd); 
pm25fd0 = center(pm25fd); 


%   ---------------------------------------------------------------
%   ------------------  Plotting the Bivariate Correlation Function 

%  --------------  First we define a fine mesh of time values. 
nfiney = 99; 
tmesh = linspace(1,T,nfiney)'; 

%   --------    get the discrete data from the lagged functions. 
hrvmat = eval_fd(hrvfd, tmesh); 
pm25mat = eval_fd(pm25fd, tmesh); 

%   --------- compute the correlations between the measures across curves. 
hrvpm25corr = corrcoef([hrvmat',pm25mat']); 
hrvpm25corr = hrvpm25corr(100:198,1:99); 
for i=2:99, hrvpm25corr(i,1:i-1) = 0; end 


%  ------------- display the correlation surface
%  ----------------Use the rotation and zoom features of Matlab's 
%  -----surface display function to examine the surface from various angles. 
%  -- COLOR plotting: a) save as encapsulated color postscript file
%                     b) in GhostView use "PS printing" option
subplot(1,1,1) 
colormap(hot) 
surf(tmesh, tmesh, hrvpm25corr') 
xlabel('\fontsize{16} s') 
ylabel('\fontsize{16} t') 
axis([0,275,0,275,-1,1]) 
axis('square') 


%  ------------------ Defining the Finite Element Basis 
%  M - number of intervals to split the time range into
%  lambda - width of the time interval
M = 25; 
lambda = T/M; 

% --- Number of time blocks going back
B = 6; 

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

trisurf(eleNodes, Si, Ti, bHat) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t'); 
title(['hrv/pm25 Data ','M = ',num2str(M),', B = ',num2str(B)]) 

% - More refined plot (color): use special function BetaEvalFD. 

svec = linspace(0,T,20*4+1)'; 
tvec = linspace(0,T,20*4+1)'; 
betaMat = BetaEvalFD(svec, tvec, bHat, M, T, lambda, ... 
eleNodes, Si, Ti, B); 
subplot(1,1,1) 
colormap(hot) 
H = imagesc(tvec, tvec, betaMat); 
xlabel('\fontsize{16} s (min)') 
ylabel('\fontsize{16} t (min)') 
axis([0,275,0,275]) 
axis('square') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 
title(['Coefficients b(s,t) ','M = ',num2str(M),', B = ',num2str(B)]) 

%%%%    -------  Computing the Fit to the Data 
%  - Set up a large super­matrix containing the approximated hrv 
%  - acceleration curves using the special function XMatrix. 
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
plot(tmesh,alphaHat) 
title('Estimated intercept function'); 

%  - Final fit to the data. 
yHat = alphaHat*ones(1,N) + yHat1; 

%  - Plot the data and the fit to the data. 
subplot(2,1,1) 
plot(tmesh, yHat) 
axis([0,T,0,150]) 
subplot(2,1,2) 
plot(tmesh, hrvmat) 
axis([0,T,0,150])

%  - Plot the residuals. 
subplot(1,1,1) 
resmat = hrvmat - yHat; 
plot(tmesh, resmat) 
axis([0,T,-50,80]) 


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
axis([0,275,-0.05,.4]) 

