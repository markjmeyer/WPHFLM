function [bHat, psiMat, yVect] = histLMFD(xfd, yMat, lag, M, B, npts)

% construction of the Y's vector:

yVect = reshape(yMat, N*M*npts, 1); 

%  Set up triangulation

eleNodes = NodeIndexation(M,B);
[Si, Ti] = ParalleloGrid(M,B,lag,tn);

%  Set up design matrix

psiMat = DesignMatrixFD(xfd, M, B, lag, eleNodes, Si, Ti, npts);

singvals = svd(psiMat);
condition = max(singvals)/min(singvals);
disp(['Condition number = ',num2str(condition)])

%  solve for b-coefficients using least squares

bHat = psiMat\yVect;

