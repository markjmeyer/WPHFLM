function BetaMat = BetaEvalFD(svec, tvec, bHat, M, tn, lambda, eleNodes, Si, Ti, ...
                              B, lag, DomainType)
%BetaCalFD.m evaluates a regression function at all pairs of the values
%  in SVEC and TVEC.  
%  Pairs not within the domain of definition are set to NaN.  
%
% Inputs :
%
% SVEC :     A strictly increasing vector of values of argument s.
% TVEC :     A strictly increasing vector of values of argument t.
% bHat :  (an array-total number of nodes X 1-of reals) estimated
%        parameters values.
% M :        number of elements spanning domain [0,TN].
% TN :       time of the last observation of the Y's.
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

if nargin < 12, DomainType = 1; end

if DomainType ~= 1 & DomainType ~= 2
    error('Argument DomainType has incorrect value');
end

if nargin < 11, lag = 0; end
lag = floor(lag);
if lag < 0 | lag >= M
    error('Lag is not between 0 and M - 1');
end
Mmlag = M - lag;

if nargin < 10, B = Mmlag; end

if B < 0 | B > Mmlag
    error('Argument B has incorrect value');
end
    
Bp1 = B + 1;

%  check argument values

index = svec >= 0.0 & svec <= tn;
svec  = svec(index);
index = tvec >= 0.0 & tvec <= tn;
tvec  = tvec(index);

Ns = length(svec);
Nt = length(tvec);

if Ns == 0
    error('No legitimate values in SVEC');
end
if Nt == 0
    error('No legitimate values in TVEC');
end

if Ns > 1
    if min(diff(svec)) <= 0.0
        error('SVEC not strictly increasing');
    end
end
if Nt > 1
    if min(diff(tvec)) <= 0.0
        error('TVEC not strictly increasing');
    end
end

BetaMat = zeros(Ns, Nt);
BetaMat = NaN;

%  ----------------------------------------------------------
%  Case of B == 0, domain = [tn, tn]:
%  ----------------------------------------------------------

if B == 0
    nNodes = Mmlag + 1;
    for it = 1:Nt
        t = tvec(is);
        sindex = find(svec == t);
        if ~isempty(sindex)
            if s == t
                m = ceil(it/Nt);
                lefttver = Si(m);
                rghttver = Si(m+1);
                tdelta   = t - lefttver;
                leftphi  = 1 - tdelta/lambda;
                rghtphi  = tdelta/lambda;
                BetaMat(sindex,it) = bHat(m)*leftPhi + bHat(m)*rghtphi
            end
        end
    end
    return;
end

%  ----------------------------------------------------------
%  Domain [0, tn]
%  ----------------------------------------------------------

if DomainType == 1
    
    %  loop through values of t

    mm2   = 0;  %  index of last left element
    m     = 1;  %  number of interval on t axis
    for it = 1:Nt
        t   = tvec(it);  %  current value of t
        if m < ceil(M*t/tn) 
            if m <= B-lag
                mm2 = mm2 + 2*(m-1) + 1;
            else
                mm2 = mm2 + 2*B;
            end
            m = m + 1;
        end
        % loop through elements spanning interval of integration
        if m <= B-lag
            %  first element in row is right
            ielem = mm2 + 1;
            %  set up basis function values for right element
            rghtNodes = eleNodes(ielem,  :);  % right node indices
            rghtsver  = Si(rghtNodes(2));
            rghttver  = Ti(rghtNodes(2));
            tdelta    = t - rghttver;
            index     = find(svec >= rghtsver & ...
                             svec <= rghtsver + tdelta);
            if ~isempty(index)
                sval      = svec(index) - rghtsver;
                tval      = (t - rghttver).*ones(length(index),1);
            % fixed the next line by adding one dot: JH 3/13/03
                rghtPhi   = bHat(rghtNodes(1)).*(tval-sval)./lambda  + ...
                           bHat(rghtNodes(2)).*    (1-tval./lambda) + ...
                           bHat(rghtNodes(3)).*       sval./lambda;
                BetaMat(it,index) = rghtPhi';
            end
            %  loop through subsequent element pairs
            for b = 1:(m-1)
                %  left element
                ielem = ielem + 1;
                leftNodes = eleNodes(ielem,:);  % left  node indices
                leftsver  = Si(leftNodes(2));
                lefttver  = Ti(leftNodes(2));
                tdelta    = t - lefttver;
                index     = find(svec >= leftsver + tdelta & ...
                                 svec <= leftsver + lambda);
                if ~isempty(index)
                    sval      = svec(index) - leftsver;
                    tval      = (t - lefttver).*ones(length(index),1);
                    leftPhi   = bHat(leftNodes(1)).*       tval./lambda  + ...
                                bHat(leftNodes(2)).*    (1-sval./lambda) + ...
                                bHat(leftNodes(3)).*(sval-tval)./lambda;
                    BetaMat(it,index) = leftPhi';
                end
                %  set up basis function values for right element
                ielem = ielem + 1;
                rghtNodes = eleNodes(ielem,:);  % right node indices
                rghtsver  = Si(rghtNodes(2));
                rghttver  = Ti(rghtNodes(2));
                tdelta    = t - rghttver;
                index     = find(svec >= rghtsver & ...
                                 svec <= rghtsver + tdelta);
                if ~isempty(index)
                    sval      = svec(index) - rghtsver;
                    tval      = (t - rghttver).*ones(length(index),1);
                    rghtPhi   = bHat(rghtNodes(1)).*(tval-sval)./lambda  + ...
                                bHat(rghtNodes(2)).*    (1-tval./lambda) + ...
                                bHat(rghtNodes(3)).*       sval./lambda;
                    BetaMat(it,index) = rghtPhi';
                end
            end
        else
            for b = 1:B
                %  set up basis for left element
                ielem     = mm2 + 2*(b-1) + 1;
                leftNodes = eleNodes(ielem,:);  % left  node indices
                leftsver  = Si(leftNodes(2));
                lefttver  = Ti(leftNodes(2));
                tdelta    = t - lefttver;
                index     = find(svec >= leftsver + tdelta & ...
                                 svec <= leftsver + lambda);
                if ~isempty(index)
                    sval      = svec(index) - leftsver;
                    tval      = (t - lefttver).*ones(length(index),1);
                    leftPhi   = bHat(leftNodes(1)).*       tval./lambda  + ...
                                bHat(leftNodes(2)).*    (1-sval./lambda) + ...
                                bHat(leftNodes(3)).*(sval-tval)./lambda;
                    BetaMat(it,index) = leftPhi';
                end
                %  set up basis function values for right element
                ielem = ielem + 1;
                rghtNodes = eleNodes(ielem,:);  % right node indices
                rghtsver  = Si(rghtNodes(2));
                rghttver  = Ti(rghtNodes(2));
                tdelta    = t - rghttver;
                index     = find(svec >= rghtsver & ...
                                 svec <= rghtsver + tdelta);
                if ~isempty(index)
                    sval      = svec(index) - rghtsver;
                    tval      = (t - rghttver).*ones(length(index),1);
                    rghtPhi   = bHat(rghtNodes(1)).*(tval-sval)./lambda  + ...
                                bHat(rghtNodes(2)).*    (1-tval./lambda) + ...
                                bHat(rghtNodes(3)).*       sval./lambda;
                    BetaMat(it,index) = rghtPhi';
                end
            end
        end
    end    
end

%  ----------------------------------------------------------
%  Domain [-tn,tn]:
%  ----------------------------------------------------------

if DomainType == 2
    
    %  loop through values of t

    for it = 1:Nt
        t   = tvec(it);  %  current value of t
        m   = ceil(M*t/tn);
        % loop through elements spanning interval of integration
        for b = 1:B
            %  set up node indices for left and right elements
            rghtElem = 2*(b + (m-1)*B);          % index of right element
            %  set up basis for left element
            leftNodes = eleNodes(rghtElem-1,:);  % left  node indices
            leftsver  = Si(leftNodes(2));
            lefttver  = Ti(leftNodes(2));
            tdelta    = t - lefttver;
            index     = find(svec >= leftsver + tdelta & ...
                             svec <= leftsver + lambda);
            if ~isempty(index)
                sval      = svec(index) - leftsver;
                tval      = (t - lefttver).*ones(length(index),1);
            % fixed the next line by adding one dot: JH 3/13/03
                leftPhi   = bHat(leftNodes(1)).*       tval./lambda  + ...
                            bHat(leftNodes(2)).*    (1-sval./lambda) + ...
                            bHat(leftNodes(3)).*(sval-tval)./lambda;
                BetaMat(it,index) = leftPhi';
            end
            %  set up basis function values for right element
            rghtNodes = eleNodes(rghtElem,  :);  % right node indices
            rghtsver  = Si(rghtNodes(2));
            rghttver  = Ti(rghtNodes(2));
            tdelta    = t - rghttver;
            index     = find(svec >= rghtsver & ...
                             svec <= rghtsver + tdelta);
            if ~isempty(index)
                sval      = svec(index);
                tval      = t*ones(length(index),1);
            % fixed the next line by adding one dot: JH 3/13/03
                rghtPhi   = bHat(rghtNodes(1)).*(tval-sval)./lambda  + ...
                            bHat(rghtNodes(2)).*    (1-tval./lambda) + ...
                            bHat(rghtNodes(3)).*       sval./lambda;
                BetaMat(it,index) = rghtPhi';
            end
        end
    end
end

