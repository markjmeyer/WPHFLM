% Analysis to be performed on dogs data (based on "utah.m"
% -------------------------------------------------------------------
% Created: 5Apr03
% Modified: 16Apr03

addpath ('c:\matlab6p5\fdaM')
addpath ('c:\matlab6p5\fdaM\examples\Brent')

%  -----------------------------------------------------------------------
%       dogs air pollution data
%  -----------------------------------------------------------------------

%  ----------------  input the data  ------------------------

pefmat = load ('dogsPef.dat'); 
% need to transpose pefmat (trouble reading ASCII file)
pefmat = pefmat';

chamat = load ('dogsChamber.dat'); 

% pefmat    - response matrix (lipmat)
% chamat   - predictor matrix (EMGmat)

%  ---------------  define number of records and time values in days
N = size(chamat,2); 
timevec = linspace(0,210,22)'; 

% -- plotting raw data
subplot(1,1,1) 
plot(timevec, pefmat) 
axis([0,T,0,1500])

%  --------------- set up data with EMG lagged by 50 milliseconds
tfine = (0:1:210)'; 
nfine = length(tfine); 

 pefmatlag = zeros(nfine,N); 
 chamatlag = zeros(nfine,N); 

for i=1:N 
    peftemp = interp1(timevec,pefmat(:,i),tfine); 
    chatemp = interp1(timevec,chamat(:,i),tfine); 
    pefmatlag(:,i) = peftemp(1:nfine);
    chamatlag(:,i) = chatemp(1:nfine); 
end 


%nfine = nfine - 50; 
%tfine = tfine(1:nfine); 
T = 210; 

%  ----------------------- Now we convert these discrete data 
%  ----------------------- to functional data objects using a B­spline basis. 
nbasis = 10; 
norder = 4; 
basis = create_bspline_basis([0,T], nbasis, norder); 
peffd = data2fd(pefmatlag, tfine, basis); 
chafd = data2fd(chamatlag, tfine, basis); 

%  ----------------  We'll also need the mean function for each variable. 
chameanfd = mean(chafd); 
pefmeanfd = mean(peffd); 

%  ---------------  centered functions
peffd0 = center(peffd); 
chafd0 = center(chafd); 


%   -------------   CAN't be done for dogs data
%   ---------------------------------------------------------------
%   ------------------  Plotting the Bivariate Correlation Function 

%  --------------  First we define a fine mesh of time values. 
nfiney = 49; 
tmesh = linspace(1,T,nfiney)'; 

%   --------    get the discrete data from the lagged functions. 
pefmat = eval_fd(peffd, tmesh); 
chamat = eval_fd(chafd, tmesh); 

%   --------- compute the correlations between the measures across curves. 
pefchacorr = corrcoef([pefmat',chamat']); 
pefchacorr = pefchacorr(50:98,1:49); 
for i=2:49, pefchacorr(i,1:i-1) = 0; end 


%  ------------- display the correlation surface
%  ----------------Use the rotation and zoom features of Matlab's 
%  -----surface display function to examine the surface from various angles. 
%  -- COLOR plotting: a) save as encapsulated color postscript file
%                     b) in GhostView use "PS printing" option
subplot(1,1,1) 
colormap(hot) 
surf(tmesh, tmesh, pefchacorr') 
xlabel('\fontsize{16} s') 
ylabel('\fontsize{16} t') 
axis([0,99,0,99,-1,1]) 
axis('square') 

%   -----------------------------------------------------------------------
%      NOT done for the dogs data above
%   -----------------------------------------------------------------------


%  ------------------ Defining the Finite Element Basis 
%  M - number of intervals to split the time range into
%  lambda - width of the time interval
M = 21; 
lambda = T/M; 

% --- Number of time blocks going back
B = 2; 

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

psiMat = DesignMatrixFD(chafd0, npts, M, eleNodes, Si, Ti, B); 

%  - check design matrix is NOT singular

singvals = svd(full(psiMat)); 
condition = max(singvals)/min(singvals); 
disp(['Condition number = ',num2str(condition)]) 

%  - vector of dependent variable values.  
yMat = eval_fd(peffd0,tpts)'; 
yVect = reshape(yMat, N*M*npts, 1);

%  - Least squares approximation that gives us our vector of regression 
%  - coefficients multiplying our basis functions. 

bHat = psiMat\yVect; 
%  - Warning: Rank deficient, rank = (depends on M and B)


%%%       ----------    Plotting the Regression Function 

trisurf(eleNodes, Si, Ti, bHat) 
xlabel('\fontsize{12} s'); 
ylabel('\fontsize{12} t'); 
title(['pef/cha Data ','M = ',num2str(M),', B = ',num2str(B)]) 

% - More refined plot (color): use special function BetaEvalFD. 

svec = linspace(0,T,13*4+1)'; 
tvec = linspace(0,T,13*4+1)'; 
betaMat = BetaEvalFD(svec, tvec, bHat, M, T, lambda, ... 
eleNodes, Si, Ti, B); 
subplot(1,1,1) 
colormap(hot) 
H = imagesc(tvec, tvec, betaMat); 
xlabel('\fontsize{16} s (sec)') 
ylabel('\fontsize{16} t (sec)') 
axis([0,210,0,210]) 
axis('square') 
Haxes = gca; 
set(Haxes,'Ydir','normal') 
title(['Coefficients b(s,t) ','M = ',num2str(M),', B = ',num2str(B)]) 

%%%%    -------  Computing the Fit to the Data 
%  - Set up a large super­matrix containing the approximated pef 
%  - acceleration curves using the special function XMatrix. 
%%%         --- This is also a lengthy calculation. 

psiArray = XMatrix(chafd, tmesh, M, eleNodes, Si, Ti, B); 

%  - Matrix of approximation values for the pef acceleration curves 

yHat1 = zeros(nfiney,N); 

for i=1:N 
Xmati = squeeze(psiArray(i,:,:))'; 
yHat1(:,i) = Xmati*bHat; 
end 

%  - Approximation is based only on the estimated regression function b(s,t). 
%  - To complete the approximation, we must get the intercept function a(t). 
%  - This requires using the mean cha curve as a model, and subtracting 
%  - the fit that this gives from the mean pef acceleration. 

psimeanArray = XMatrix(chameanfd, tmesh, M, eleNodes, Si, Ti, B); 
yHatMean = squeeze(psimeanArray)'*bHat; 
pefmeanvec = eval_fd(pefmeanfd, tmesh); 
alphaHat = pefmeanvec - yHatMean; 

%  - Plot the intercept function. 
plot(tmesh,alphaHat) 
title('Estimated intercept function'); 

%  - Final fit to the data. 
yHat = alphaHat*ones(1,N) + yHat1; 

%  - Plot the data and the fit to the data. 
subplot(2,1,1) 
plot(tmesh, yHat) 
axis([0,T,200,350]) 
subplot(2,1,2) 
plot(tmesh, pefmat) 
axis([0,T,0,1500])

%  - Plot the residuals. 
subplot(1,1,1) 
resmat = pefmat - yHat; 
plot(tmesh, resmat) 
axis([0,T,-4,4]) 


%%%%%       ----- Assessing the Fit 
%  - Error sum of squares function. 
SSE = sum(resmat.^2,2); 

%  - Benchmark against which we can assess this fit, we need to get the 
%  - corresponding error sum of squares function when the model is simply 
%  - the mean pef acceleration curve. 

pefmat0 = eval_fd(peffd0, tmesh); 
SSY = sum(pefmat0.^2,2); 

%  - compute a squared multiple correlation function and plot it. 
%  - Don't be suprised, though, to see it go below zero; 
%  - the fit from the mean is not embedded within the fit by the model. 
RSQ = (SSY - SSE)./SSY; 
subplot(1,1,1) 
plot(tmesh,RSQ); 
axis([0,210,-0.05,.1]) 

