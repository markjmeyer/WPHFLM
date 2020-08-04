% 'DomainTriangulation.m' discretizes the domain of integration.
% ------------------------------------------------------------ %

% Inputs :

% M :   (an integer) number of intervalles in which the domain 
%       boundaries are divided.
% B :   (an integer) width of the domain of integration in terms
%       of number elements.
% lag : (a real) lag.
% tn :  (a real) time of the last observation of the Y's.

% Outputs :

% eleNodes:(an array-3 X total number of nodes-of integers) nodes
%        indexation.
% [Si, Ti]: (an array-2 X total number of nodes-of reals) nodes
%        cartisians coordinates.

% ------------------------------------------------------------ %
function [eleNodes, Si, Ti] = DomainTriangulation(M,B,lag,tn)
    eleNodes = NodeIndexation(M,B);
    [Si, Ti] = ParalleloGrid(M,B,lag,tn);
