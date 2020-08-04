function psiArray = XMatrix(xfd, argvals, M, eleNodes, Si, Ti, B, lag, DomainType)
%XMatrix.m constructs the three-way design matrix used to compute 
%   the approximation of y_i(t_j), i=1,...,N; j=1,...,J
%  The design matrix PSIARRAY is N by K by J where K is the number of nodes.
%
% Inputs :
%
% XFD :      functional data object for the independent variable functions.
% ARGVALS:   J order values of t.
% M :        number of elements spanning domain [0,TN].
% eleNodes:  matrix of node indices for each element
% SI :       node x-coordinates.
% TI :       node y-coordinates.
% B :        number of elements spanning the domain of integration
% LAG : difference between M and the index of the upper element.
%       By default, this is 0.
% DomainType : Type of domain as follows
%      DomainType = 1:  Triangular,    domain of s = [ 0,T]
%  In the following figures, elements are lettered and nodes are numbered.
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
% PSIARRAY :  N by number of nodes by number of t-values design matrix
 
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

%  check XFD

if ~isa_fd(xfd)
    error('XFD is not a functional data object');
end

%  check argvals

if any(diff(argvals) <= 0) 
    error('ARGVALS not strictly increasing');
end

%  get some dimensions and the range of t-values

xbasis = getbasis(xfd);
xrange = getbasisrange(xbasis);
tn     = xrange(2);
if any(argvals < 0) | any(argvals > tn)
    error('ARGVALS not within correct range');
end
lambda = tn/M;              %  width of an element

%  get sample size

xcoef = getcoef(xfd);
N     = size(xcoef,2);
onesN = ones(1,N);

%  number of t values

J = length(argvals);

%  Case of B == 0, domain = [tn, tn]:

if B == 0
    nNodes = Mmlag + 1;
    psiArray = zeros(N,nNodes,J);
    for m = 1:Mmlag
        %  select values of t within this element band
        if m < Mmlag
            indext = find(argvals >= (m-1)*lambda & argvals <  m*lambda);
        else
            indext = find(argvals >= (m-1)*lambda & argvals <= m*lambda);
        end
        if ~isempty(indext)
            lefttver = Si(m);
            rghttver = Si(m+1);
            for it = 1:length(indext)
                t      = argvals(indext(it));
                xvec = eval(xfd,t)';
                tdelta   = t - lefttver;
                leftphi  = 1 - tdelta/lambda;
                psiArray(:,m,  indext(it)) = ...
                psiArray(:,m,  indext(it)) + ...
                                             leftphi.*xvec;
                rghtphi  = tdelta/lambda;
                psiArray(:,m+1,indext(it)) = ...
                psiArray(:,m+1,indext(it)) + ...
                                             rghtphi.*xvec;
            end
        end
    end
    return;
end

coef = zeros(2,3);

%  Domain [0,tn]:

if DomainType == 1
    
%  set up design matrix dimensions
    
nNodes   = (B+1)*(Mmlag - B/2 + 1);  %  number of nodes
psiArray = zeros(N,nNodes,J);
    
%  loop through values of elements

ielem = 0;
mm3   = 0;
mm4   = 1;
for m = 1:B
    %  select values of t within this element band
    if m < Mmlag
        indext = find(argvals >= (m-1)*lambda & argvals <  m*lambda);
    else
        indext = find(argvals >= (m-1)*lambda & argvals <= m*lambda);
    end
    if ~isempty(indext)
        %  first element in row is right
        ielem     = ielem + 1;
        rghtNodes = eleNodes(ielem,:);  % right node indices
        rghtsver  = Si(rghtNodes(2));   % right element lower s bound
        rghttver  = Ti(rghtNodes(2));   % right element lower t bound
        for it = 1:length(indext)
            t      = argvals(indext(it));
            tdelta = t - rghttver;
            tval   = [tdelta; tdelta];
            slo    = rghtsver;
            shi    = rghtsver + tdelta;
            srng   = shi - slo;
            if abs(srng) > 1e-7*lambda
                sval     = [0;  tdelta];
                rghtPhi  = [(tval-sval)./lambda, ...
                                 1-tval./lambda, ...
                                   sval./lambda];
                coef(2,:) = (rghtPhi(2,:)-rghtPhi(1,:))./srng;
                coef(1,:) = rghtPhi(1,:) - coef(2,:).*slo;
                rangeR    = [slo, shi];
                basisR    = create_power_basis(rangeR, 2);
                linearfdR = fd(coef, basisR);
                psiArray(:,rghtNodes,indext(it)) = ...
                psiArray(:,rghtNodes,indext(it)) + ...
                                   inprod(xfd, linearfdR, 0, 0, rangeR);
            end
        end
        % loop through elements spanning interval of integration
        for b = 1:(m-1)
            ielem = ielem + 2;
            leftNodes = eleNodes(ielem-1,:);   % left  node indices
            rghtNodes = eleNodes(ielem,  :);   % right node indices
            leftsver  = Si(leftNodes(2));   % left  element lower s bound
            lefttver  = Ti(leftNodes(2));   % left  element lower t bound
            rghtsver  = Si(rghtNodes(2));   % right element lower s bound
            rghttver  = Ti(rghtNodes(2));   % right element lower t bound
            for it = 1:length(indext)
                t      = argvals(indext(it));
                tdelta = t - lefttver;
                tval   = [tdelta; tdelta];
                %  increment psiArray for left element
                slo    = leftsver + tdelta;
                shi    = leftsver + lambda;
                srng   = shi - slo;
                if abs(srng) > 1e-7*lambda
                    sval     = [tdelta; lambda];
                    leftPhi  = [       tval./lambda, ...
                                     1-sval./lambda, ...
                                (sval-tval)./lambda];
                    coef(2,:) = (leftPhi(2,:) - leftPhi(1,:))./srng;
                    coef(1,:) = leftPhi(1,:) - coef(2,:).*slo;
                    rangeL    = [slo, shi];
                    basisL    = create_power_basis(rangeL, 2);
                    linearfdL = fd(coef, basisL);
                    psiArray(:,leftNodes,indext(it)) = ...
                    psiArray(:,leftNodes,indext(it)) + ...
                                   inprod(xfd, linearfdL, 0, 0, rangeL);
                end
                %  increment psiArray for right element
                slo  = rghtsver;
                shi  = rghtsver + tdelta;
                srng = shi - slo;
                if abs(srng) > 1e-7*lambda
                    sval     = [0;  tdelta];
                    rghtPhi  = [(tval-sval)./lambda, ...
                                     1-tval./lambda, ...
                                       sval./lambda];
                    coef(2,:) = (rghtPhi(2,:)-rghtPhi(1,:))./srng;
                    coef(1,:) = rghtPhi(1,:) - coef(2,:).*slo;
                    rangeR    = [slo, shi];
                    basisR    = create_power_basis(rangeR, 2);
                    linearfdR = fd(coef, basisR);
                    psiArray(:,rghtNodes,indext(it)) = ...
                    psiArray(:,rghtNodes,indext(it)) + ...
                                   inprod(xfd, linearfdR, 0, 0, rangeR);
                end
            end
        end
    end
    mm3 = mm3 + m;
    mm4 = mm4 + m + 1;
end
for m = (B+1):Mmlag
    %  select values of t within this element band
    if m < Mmlag
        indext = find(argvals >= (m-1)*lambda & argvals <  m*lambda);
    else
        indext = find(argvals >= (m-1)*lambda & argvals <= m*lambda);
    end
    if ~isempty(indext)
        for b = 1:B
            ielem = ielem + 2;
            leftNodes = eleNodes(ielem-1,:);   % left  node indices
            rghtNodes = eleNodes(ielem,  :);   % right node indices
            leftsver  = Si(leftNodes(2));   % left  element lower s bound
            lefttver  = Ti(leftNodes(2));   % left  element lower t bound
            rghtsver  = Si(rghtNodes(2));   % right element lower s bound
            rghttver  = Ti(rghtNodes(2));   % right element lower t bound
            for it = 1:length(indext)
                t      = argvals(indext(it));
                tdelta = t - lefttver;
                tval   = [tdelta; tdelta];
                %  increment psiArray for left element
                slo    = leftsver + tdelta;
                shi    = leftsver + lambda;
                srng   = shi - slo;
                if abs(srng) > 1e-7*lambda
                    sval     = [tdelta; lambda];
                    leftPhi  = [       tval./lambda, ...
                                     1-sval./lambda, ...
                                (sval-tval)./lambda];
                    coef(2,:) = (leftPhi(2,:) - leftPhi(1,:))./srng;
                    coef(1,:) = leftPhi(1,:) - coef(2,:).*slo;
                    rangeL    = [slo, shi];
                    basisL    = create_power_basis(rangeL, 2);
                    linearfdL = fd(coef, basisL);
                    psiArray(:,leftNodes,indext(it)) = ...
                    psiArray(:,leftNodes,indext(it)) + ...
                                   inprod(xfd, linearfdL, 0, 0, rangeL);
                end
                %  increment psiArray for right element
                slo  = rghtsver;
                shi  = rghtsver + tdelta;
                srng = shi - slo;
                if abs(srng) > 1e-7*lambda
                    sval     = [0;  tdelta];
                    rghtPhi  = [(tval-sval)./lambda, ...
                                     1-tval./lambda, ...
                                       sval./lambda];
                    coef(2,:) = (rghtPhi(2,:)-rghtPhi(1,:))./srng;
                    coef(1,:) = rghtPhi(1,:) - coef(2,:).*slo;
                    rangeR    = [slo, shi];
                    basisR    = create_power_basis(rangeR, 2);
                    linearfdR = fd(coef, basisR);
                    psiArray(:,rghtNodes,indext(it)) = ...
                    psiArray(:,rghtNodes,indext(it)) + ...
                                   inprod(xfd, linearfdR, 0, 0, rangeR);
                end
            end
        end
    end
    mm3 = mm3 + B + 1;
    mm4 = mm4 + B + 1;
end

end

%  Domain [-tn,tn]:

if DomainType == 2
    
%  set up design matrix dimensions
    
nNodes   = (Mmlag+1)*(B+1);  %  number of nodes
psiArray = zeros(N,nNodes,J);
    
%  loop through values of elements

for m = 1:Mmlag
    %  select values of t within this element band
    if m < Mmlag
        indext = find(argvals >= (m-1)*lambda & argvals <  m*lambda);
    else
        indext = find(argvals >= (m-1)*lambda & argvals <= m*lambda);
    end
    if ~isempty(indext)
        % loop through elements spanning interval of integration
        for b = 1:B
            %  At each left boundary there are two elements, 
            %    the left and right
            %  set up X values for left and right elements
            %  set up node indices for left and right elements
            rghtElem  = 2*(b + (m-1)*B);          % index of right element
            leftNodes = eleNodes(rghtElem-1,:);   % left  node indices
            rghtNodes = eleNodes(rghtElem,  :);   % right node indices
            leftsver  = Si(leftNodes(2));   % left  element lower s bound
            lefttver  = Ti(leftNodes(2));   % left  element lower t bound
            rghtsver  = Si(rghtNodes(2));   % right element lower s bound
            rghttver  = Ti(rghtNodes(2));   % right element lower t bound
            for it = 1:length(indext)
                t      = argvals(indext(it));
                tdelta = t - lefttver;
                tval   = [tdelta; tdelta];
                %  increment psiArray for left element
                slo  = leftsver + tdelta;
                shi  = leftsver + lambda;
                srng = shi - slo;
                if abs(srng) > 1e-7*lambda
                    sval     = [tdelta; lambda];
                    leftPhi  = [       tval./lambda, ...
                                     1-sval./lambda, ...
                                (sval-tval)./lambda];
                    coef(2,:) = (leftPhi(2,:) - leftPhi(1,:))./srng;
                    coef(1,:) = leftPhi(1,:) - coef(2,:).*slo;
                    rangeL    = [slo, shi];
                    basisL    = create_power_basis(rangeL, 2);
                    linearfdL = fd(coef, basisL);
                    psiArray(:,leftNodes,indext(it)) = ...
                    psiArray(:,leftNodes,indext(it)) + ...
                                   inprod(xfd, linearfdL, 0, 0, rangeL);
                end
                %  increment psiArray for right element
                slo  = rghtsver;
                shi  = rghtsver + tdelta;
                srng = shi - slo;
                if abs(srng) > 1e-7*lambda
                    sval     = [0;  tdelta];
                    rghtPhi  = [(tval-sval)./lambda, ...
                                     1-tval./lambda, ...
                                       sval./lambda];
                    coef(2,:) = (rghtPhi(2,:)-rghtPhi(1,:))./srng;
                    coef(1,:) = rghtPhi(1,:) - coef(2,:).*slo;
                    rangeR    = [slo, shi];
                    basisR    = create_power_basis(rangeR, 2);
                    linearfdR = fd(coef, basisR);
                    psiArray(:,rghtNodes,indext(it)) = ...
                    psiArray(:,rghtNodes,indext(it)) + ...
                                   inprod(xfd, linearfdR, 0, 0, rangeR);
                end
            end
        end
    end
end

end

%toc;
