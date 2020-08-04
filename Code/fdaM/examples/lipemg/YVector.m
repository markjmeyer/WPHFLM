function yVect = YVector(ytime,yfunct,tpts)
%YVECTOR creates the vector of the points were the 
% functions are evaluated:
%
% Inputs :
%
% ytime : Y's time scale.
% yfunct: Y's functions.
% tpts :  (an array-1 X total number of nodes-of reals) values of
%         used to construct the design matrix.
%
% Outputs :
%
% YVector: vector of N*nbtEle Y-values.

%  Last modified:  7 July 2001

N = size(yfunct,1);
nbPts = length(tpts);
yMat = zeros(N,nbPts);
for r = 1:N
    yMat(r,:) = interp1(ytime,yfunct(r,:),tpts);
end
% Data are regrouped by t values.
yVect = reshape(yMat,N*nbPts,1);
    
    
    
    




