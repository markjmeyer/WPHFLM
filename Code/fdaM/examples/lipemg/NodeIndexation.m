function eleNodes = NodeIndexation(M, B, lag, DomainType)
%NodeIndexation.m computes indices of the three nodes 
%   defining each triangular element.
%
% Inputs :
%
% M :   number of intervals into which the domain 
%         boundaries are divided.
% B :   width of the domain of integration in terms
%         of number elements.  Default value is M.
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
%  The total number of elements and nodes are:
%      DomainType = 1:  Elements ... 2*B(M-1), and nodes ... (B+1)(M-1+B/2).
%      DomainType = 2:  Elements ... 2*M*B,    and nodes ... (M+1)*(B+1).
%
%  However, when B = 0 and lag = 0, the two domain types are equivalent, 
%    the number of elements is M and the number of nodes is M+1. 
%  The default value for DomainType is 1.
%
% Outputs :
%
% eleNodes:  a matrix with 
%          number rows = total number of elements
%          number cols = 3
% Each row contains the node numbers defining the element, numbered starting
%    with the upper left node and going anti-clockwise.
%
%  However, when B = 0, the first column contains the left node, 
%                       the second the right node, and the third column is 0.
%
% For the above Domaintype = 2 diagram, eleNodes returns
%   
%  4  1  2
%  4  2  5
%  5  2  3
%  5  3  6
%  7  4  5
%  7  5  8
%  8  5  6
%  8  6  9

%  Last modified:  3 January 2002

if nargin < 4, DomainType = 1; end

if DomainType ~= 1 & DomainType ~= 2
    error('Argument DomainType has incorrect value');
end

if nargin < 3, lag = 0; end
lag = floor(lag);
if lag < 0 | lag >= M
    error('Lag is not between 0 and M - 1');
end
Mmlag = M - lag;

if nargin < 2, B = Mmlag; end

if B < 0 | B > Mmlag
    error('Argument B has incorrect value');
end
    
Bp1 = B + 1;

%  Domain [tn, tn]:

if B == 0
    nbEle = Mmlag;
    eleNodes = zeros(nbEle,3);
    eleNodes(:,1) = (1:(Mmlag)  )';
    eleNodes(:,2) = (2:(Mmlag+1))';
    return;
end

%  Domain [0, tn]

if DomainType == 1
    nbEle    = Mmlag^2 - (Mmlag-B)^2;
    eleNodes = zeros(nbEle,3);
    ielem = 0;
    mm3   = 0;
    mm4   = 1;
    for m = 1:B
        %  first element in row is right
        ielem = ielem + 1;
        eleNodes(ielem,:) = [mm4+1,mm3+1,mm4+2];
        for b = 1:(m-1)
            bp1 = b + 1;
            %  left element
            ielem = ielem + 1;
            eleNodes(ielem,:) = [mm4+bp1,mm3+b,  mm3+bp1];
            %  right element
            ielem = ielem + 1;
            eleNodes(ielem,:) = [mm4+bp1,mm3+bp1,mm4+b+2];
        end
        mm3 = mm3 + m;
        mm4 = mm4 + m + 1;
    end
    for m = (B+1):Mmlag
        for b = 1:B
            bp1 = b + 1;
            %  left element
            ielem = ielem + 1;
            eleNodes(ielem,:) = [mm4+b,mm3+b,  mm3+bp1];
            %  right element
            ielem = ielem + 1;
            eleNodes(ielem,:) = [mm4+b,mm3+bp1,mm4+bp1];
        end
        mm3 = mm3 + B + 1;
        mm4 = mm4 + B + 1;
    end
end

%  Domain [-tn,tn]:

if DomainType == 2
    nbEle = 2*Mmlag*B;
    eleNodes = zeros(nbEle,3);
    for m = 1:Mmlag
        mm1 = m - 1;
        m2  = 2*mm1*B;
        m3  = mm1*Bp1;
        m4  = m*Bp1;
        for b = 1:B
            bp1 = b + 1;
            %  left element
            b2  = 2*b + m2 - 1;
            eleNodes(b2,:) = [m4+b,m3+b  ,m3+bp1];
            %  right element
            b2  = b2 + 1;
            eleNodes(b2,:) = [m4+b,m3+bp1,m4+bp1];
       end
    end
end

