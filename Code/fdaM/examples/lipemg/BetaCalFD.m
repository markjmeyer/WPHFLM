function psiMat = DesignMatrixFD(xfd, npts, M, eleNodes, Si, Ti, ...
                                 B, lag, DomainType)
%DesignMatrix.m constructs the design matrix used to compute 
%   node coefficients
%
% Inputs :
%
% XFD:       A functional data object for the independent variable
% NPTS :     number of t points per element used to 
%               construct the design matrix.
% M :        number of elements spanning domain [0,TN].
% eleNodes:  matrix of node indices for each element
% SI :       node x-coordinates.
% TI :       node y-coordinates.
% B :        number of elements spanning the domain of integration
% LAG : difference between M and the index of the upper element.
%       By default, this is 0.
% DomainType : Type of domain as follows
%      DomainType = 1:  Triangular,    domain of s = [ 0,T]
%  In these figures, elements are lettered and nodes are numbered.
%  Both elements and nodes are numbered row-wise from left to right
%    and from the bottom to the top.
%
%          4--5--6
%          |b/|d/
%          |/c|/ 
%          2--3
%          |a/
%          |/
%          1
%
%      DomainType = 2:  Parallelogram, domain of s = [-T,T]
%  For M = 2 and B = 2, there are 8 elements and 9 nodes.
%
%          7--8--9
%         /|f/|h/
%        /e|/g|/ 
%       4--5--6
%      /|b/|d/
%     /a|/c|/
%    1--2--3
%
%  The number of elements and nodes are:
%      DomainType = 1:  Elements 2*B(M-1), and nodes (B+1)(M+1-B/2).
%      DomainType = 2:  Elements 2*M*B,    and nodes (M+1)*(B+1).
%
%  However, when B = 0, the two domain types are equivalent, the number of 
%    elements is M and the number of nodes is M+1. 
%  The default value for DomainType is 1.
%
% Outputs :
%
% PSIMAT :  N*M*NPTS by number of nodes design matrix 
%             returned in sparse storage mode.
 
%  Last modified:  3 January 2002

%  check arguments

%tic;

if nargin < 9, DomainType = 1; end

if DomainType ~= 1 & DomainType ~= 2
    error('Argument DomainType has incorrect value');
end

if nargin < 8, lag = 0; end
lag = floor(lag);
if lag < 0 | lag >= M
    error('Lag is not between 0 and M - 1');
end
Mmlag = M - lag;

if nargin < 7, B = Mmlag; end

if B < 0 | B > Mmlag
    error('Argument B has incorrect value');
end
    
Bp1 = B + 1;

if ~isa_fd(xfd)
    error('XFD is not a functional data object');
end

%  get number of curves and upper limit of integration.

xcoef  = getcoef(xfd);
N      = size(xcoef,2);
xbasis = getbasis(xfd);
xrng   = getbasisrange(xbasis);
tn     = xrng(2);

%  set up values of t to be used for the design matrix
    
ntpts  = Mmlag*npts;
lambda = tn/M;         %  width of an element
delta  = lambda/(2*npts);
tpts   = linspace(delta, tn-delta, ntpts);

%  ----------------------------------------------------------
%  Case of B == 0, domain = [tn, tn]:
%  ----------------------------------------------------------

if B == 0
    nNodes = Mmlag + 1;
    psiMat = zeros(N*ntpts,nNodes);
    for it = 1:ntpts
        t = tpts(it);  %  current value of t
        m = ceil(it/npts);
        index    = (1:N) + (it-1)*N;
        lefttver = Si(m);
        rghttver = Si(m+1);
        tdelta   = t - lefttver;
        leftphi  = 1 - tdelta/lambda;
        xvec     = eval(xfd,t)';
        psiMat(index,m)   = psiMat(index,m)   + leftphi.*xvec;
        rghtphi  = tdelta/lambda;
        psiMat(index,m+1) = psiMat(index,m+1) + rghtphi.*xvec;
    end
    return;
end

coef = zeros(2,3);

%  ----------------------------------------------------------
%  Domain [0, tn]
%  ----------------------------------------------------------

if DomainType == 1
    
    %  set up design matrix dimensions
    
    nNodes = (B+1)*(Mmlag - B/2 + 1);  %  number of nodes
    psiMat = zeros(N*ntpts,nNodes);
    
    %  loop through values of t

    mm2   = 0;
    mm3   = 0;
    mm4   = 1;
    m     = 1;
    for it = 1:ntpts
        t   = tpts(it);  %  current value of t
        if m < ceil(it/npts)
            if m <= B-lag
                mm2 = mm2 + 2*(m-1) + 1;
                mm3 = mm3 + m;
                mm4 = mm4 + m + 1;
            else
                mm2 = mm2 + 2*B;
                mm3 = mm3 + B + 1;
                mm4 = mm4 + B + 1;
            end
            m = m + 1;
        end
        % loop through elements spanning interval of integration
        index = (1:N) + (it-1)*N;
        if m <= B-lag
            %  first element in row is right
            ielem = mm2 + 1;
            %  set up basis function values for right element
            rghtNodes = eleNodes(ielem,  :);  % right node indices
            rghtsver  = Si(rghtNodes(2));
            rghttver  = Ti(rghtNodes(2));
            tdelta    = t - rghttver;
            sval      = [0;      tdelta];
            tval      = [tdelta; tdelta];
            rghtPhi   = [(tval-sval)./lambda, ...
                         1-tval./lambda, ...
                           sval./lambda];
            sup       = rghtsver + tdelta;    
            sdn       = rghtsver;
            coef(2,:) = (rghtPhi(2,:) - rghtPhi(1,:))./(sup-sdn);
            coef(1,:) = rghtPhi(1,:) - coef(2,:).*sdn;
            rangeR    = [sdn,sup];
            basisR    = create_polynomial_basis(rangeR, 2);
            linearfdR = fd(coef, basisR);
            rghtIntg = inprod(xfd, linearfdR, 0, 0, rangeR);
            psiMat(index,rghtNodes) = psiMat(index,rghtNodes) + rghtIntg;
            %  loop through subsequent element pairs
            for b = 1:(m-1)
                %  left element
                ielem = ielem + 1;
                leftNodes = eleNodes(ielem,:);  % left  node indices
                leftsver  = Si(leftNodes(2));
                lefttver  = Ti(leftNodes(2));
                tdelta    = t - lefttver;
                sval      = [tdelta; lambda];
                tval      = [tdelta; tdelta];
                leftPhi   = [       tval./lambda, ...
                                  1-sval./lambda, ...
                             (sval-tval)./lambda];
                sup       = leftsver + lambda;    
                sdn       = leftsver + tdelta;
                coef(2,:) = (leftPhi(2,:)-leftPhi(1,:))./(sup-sdn);
                coef(1,:) = leftPhi(1,:) - coef(2,:).*sdn;
                rangeL    = [sdn,sup];
                basisL    = create_polynomial_basis(rangeL, 2);
                linearfdL = fd(coef, basisL);
                %  set up basis function values for right element
                ielem = ielem + 1;
                rghtNodes = eleNodes(ielem,:);  % right node indices
                rghtsver  = Si(rghtNodes(2));
                sval      = [0;      tdelta];
                tval      = [tdelta; tdelta];
                rghtPhi   = [(tval-sval)./lambda, ...
                                  1-tval./lambda, ...
                                    sval./lambda];
                sup       = rghtsver + tdelta;    
                sdn       = rghtsver;
                coef(2,:) = (rghtPhi(2,:) - rghtPhi(1,:))./(sup-sdn);
                coef(1,:) = rghtPhi(1,:) - coef(2,:).*sdn;
                rangeR    = [sdn,sup];
                basisR    = create_polynomial_basis(rangeR, 2);
                linearfdR = fd(coef, basisR);
                %  approximate integrals for left and right elements by
                %    trapezoidal rule
                leftIntg = inprod(xfd, linearfdL, 0, 0, rangeL);
                rghtIntg = inprod(xfd, linearfdR, 0, 0, rangeR);
                psiMat(index,leftNodes) = psiMat(index,leftNodes) + ...
                                          leftIntg;
                psiMat(index,rghtNodes) = psiMat(index,rghtNodes) + ...
                                          rghtIntg;
            end
        else
            for b = 1:B
                %  set up basis for left element
                ielem     = mm2 + 2*(b-1) + 1;
                leftNodes = eleNodes(ielem,:);  % left  node indices
                leftsver  = Si(leftNodes(2));
                lefttver  = Ti(leftNodes(2));
                tdelta    = t - lefttver;
                sval      = [tdelta; lambda];
                tval      = [tdelta; tdelta];
                leftPhi   = [       tval./lambda, ...
                                  1-sval./lambda, ...
                             (sval-tval)./lambda];
                sup       = leftsver + lambda;    
                sdn       = leftsver + tdelta;
                coef(2,:) = (leftPhi(2,:)-leftPhi(1,:))./(sup-sdn);
                coef(1,:) = leftPhi(1,:) - coef(2,:).*sdn;
                rangeL    = [sdn,sup];
                basisL    = create_polynomial_basis(rangeL, 2);
                linearfdL = fd(coef, basisL);
                %  set up basis function values for right element
                ielem = ielem + 1;
                rghtNodes = eleNodes(ielem,:);  % right node indices
                rghtsver  = Si(rghtNodes(2));
                sval      = [0;      tdelta];
                tval      = [tdelta; tdelta];
                rghtPhi   = [(tval-sval)./lambda, ...
                             1-tval./lambda, ...
                               sval./lambda];
                sup       = rghtsver + tdelta;    
                sdn       = rghtsver;
                coef(2,:) = (rghtPhi(2,:) - rghtPhi(1,:))./(sup-sdn);
                coef(1,:) = rghtPhi(1,:) - coef(2,:).*sdn;
                rangeR    = [sdn,sup];
                basisR    = create_polynomial_basis(rangeR, 2);
                linearfdR = fd(coef, basisR);
                %  approximate integrals for left and right elements by
                %    trapezoidal rule
                leftIntg = inprod(xfd, linearfdL, 0, 0, rangeL);
                rghtIntg = inprod(xfd, linearfdR, 0, 0, rangeR);
                psiMat(index,leftNodes) = psiMat(index,leftNodes) + ...
                                          leftIntg;
                psiMat(index,rghtNodes) = psiMat(index,rghtNodes) + ...
                                          rghtIntg;
            end
        end
    end    
end

%  ----------------------------------------------------------
%  Domain [-tn,tn]:
%  ----------------------------------------------------------

if DomainType == 2
    
    %  set up design matrix dimensions
    
    nNodes = (Mmlag+1)*(B+1);  %  number of nodes
    psiMat = zeros(N*ntpts,nNodes);
    
    %  loop through values of t

    for it = 1:ntpts
        t   = tpts(it);  %  current value of t
        m   = ceil(it/npts);
        % loop through elements spanning interval of integration
        index = (1:N) + (it-1)*N;
        for b = 1:B
            %  set up node indices for left and right elements
            rghtElem = 2*(b + (m-1)*B);          % index of right element
            %  set up basis for left element
            leftNodes = eleNodes(rghtElem-1,:);  % left  node indices
            leftsver  = Si(leftNodes(2));
            lefttver  = Ti(leftNodes(2));
            tdelta    = t - lefttver;
            sup       = leftsver + lambda;    
            sdn       = leftsver + tdelta;
            rangeL    = [sdn,sup];
            basisL    = create_polynomial_basis(rangeL, 2);
            sval      = [tdelta; lambda];
            tval      = [tdelta; tdelta];
            leftPhi   = [       tval./lambda, ...
                              1-sval./lambda, ...
                         (sval-tval)./lambda];
            coef(2,:) = (leftPhi(2,:)-leftPhi(1,:))./(sup-sdn);
            coef(1,:) = leftPhi(1,:) - coef(2,:).*sdn;
            linearfdL = fd(coef, basisL);
            %  set up basis function values for right element
            rghtNodes = eleNodes(rghtElem,  :);  % right node indices
            rghtsver  = Si(rghtNodes(2));
            sup       = rghtsver + tdelta;    
            sdn       = rghtsver;
            rangeR    = [sdn,sup];
            basisR    = create_polynomial_basis(rangeR, 2);
            sval      = [0;      tdelta];
            tval      = [tdelta; tdelta];
            rghtPhi   = [(tval-sval)./lambda, ...
                             1-tval./lambda, ...
                               sval./lambda];
            coef(2,:) = (rghtPhi(2,:) - rghtPhi(1,:))./(sup-sdn);
            coef(1,:) = rghtPhi(1,:) - coef(2,:).*sdn;
            linearfdR = fd(coef, basisR);
            %  approximate integrals for left and right elements by
            %    trapezoidal rule
            leftIntg = inprod(xfd, linearfdL, 0, 0, rangeL);
            rghtIntg = inprod(xfd, linearfdR, 0, 0, rangeR);
            psiMat(index,leftNodes) = psiMat(index,leftNodes) + leftIntg;
            psiMat(index,rghtNodes) = psiMat(index,rghtNodes) + rghtIntg;
        end
    end
end

psiMat = sparse(psiMat);
%toc;