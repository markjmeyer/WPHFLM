function [Si, Ti] = ParalleloGrid(M, tn, B, lag, DomainType)
%ParalleloGrid creates the grid for a parallel domain:
%  See function NodeIndexation for numbering system for elements and nodes.
%
% Inputs :
% M :   number of intervals into which the domain [0, TN] is divided.
% TN :  time of the last observation of the Y's.
% B :   number of elements spanning the interval of integration (1 <= B <= M)
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
%  When lag = 0, the total number of elements and nodes are:
%      DomainType = 1:  Elements ... 2*B(M-1), and nodes ... (B+1)(M-1+B/2).
%      DomainType = 2:  Elements ... 2*M*B,    and nodes ... (M+1)*(B+1).
%
%  However, when B = 0, the two domain types are equivalent, the number of 
%    elements is M and the number of nodes is M+1. 
%  The default value for DomainType is 1.
%
% Outputs :
% Si: abscissa values for nodes
% Ti: ordinate values for nodes

%  Last modified:  3 January 2002

if nargin < 5, DomainType = 1; end

if DomainType ~= 1 & DomainType ~= 2
    error('Argument DomainType has incorrect value');
end

if nargin < 4, lag = 0; end
lag = floor(lag);
if lag < 0 | lag >= M
    error('Lag is not between 0 and M - 1');
end
Mmlag = M - lag;

if nargin < 3, B = Mmlag; end

if B < 0 | B > Mmlag
    error('Argument B has incorrect value');
end

lambda = tn/M;  %  width of an element

%  Domain [tn, tn]:

if B == 0
    Si = linspace(0,tn,Mmlag+1)';
    Ti = Si;
    return;
end

%  Domain [0, tn]

if DomainType == 1
    Si = 0;
    Ti = 0;
    for m=1:B
        Si = [Si, linspace(0, m*lambda, m+1)];
        Ti = [Ti, (m*lambda).*ones(1,m+1)];
    end
    for m = (B+1):Mmlag
        Si = [Si, linspace((m-B)*lambda, m*lambda, B+1)];
        Ti = [Ti, (m*lambda).*ones(1,B+1)];
    end
end

%  Domain [-tn,tn]:

if DomainType == 2
    s0 = (-lag - B)*lambda;  %  beginning of interval of integration
    is = linspace(s0,-lag*lambda,B+1);  %  abscissa boundaries of elements
    Si = is;
    Ti = zeros(1,B+1);
    for m = 1:M
        Si = [Si,is+m*lambda];
        Ti = [Ti,ones(1,B+1)*m*lambda];
    end
end
    
    
    
    
