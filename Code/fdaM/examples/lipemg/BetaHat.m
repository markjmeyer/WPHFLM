function betaHat = BetaHat(s, t, s0, B, eleNodes, Si, Ti, lambda, bHat)
% BetaHat.m computes the estimate of the regression function 
% for a given point (s,t).
% Inputs :
% S :        s-coordinate.
% T :        t-coordinate.
% S0 :       the first s-value.
% B :        number of elements spanning the domain of integration
% eleNodes:  matrix of node indices for each element
% SI :       node x-coordinates.
% TI :       node y-coordinates.
% lambda : element width.
% bHat :  (an array-total number of nodes X 1-of reals) estimated
%        parameters values.
% Outputs :
% betaHat : estimated value of Beta for the given point (s,t).

%  Last modified 25 July 2001

    %eps = 0.000001;
    eps = 0;
    if t == 0
        m = 1;
    elseif abs(rem(t,lambda)) < eps
        m = round(t/lambda);
    else 
        m = fix(t/lambda)+1;
    end
    
    if     abs((s-s0)-t)             < eps
        b = 1;
    elseif abs(rem((s-s0-t),lambda)) < eps
        b = round((s-s0-t)/lambda);
    else 
        b = fix((s-s0-t)/lambda)+1; 
    end
    
    b2 = 2*b+(m-1)*2*B;  % index of left element
    v = Si(eleNodes(b2,1));

    if ((s-s0) == t | v > s)
        nodes = eleNodes(b2-1,:);  % right element
    else
        nodes = eleNodes(b2,:);    % left  element
    end
    M  = [1,1,1;Si(nodes)-s0;Ti(nodes)];
    M1 = M;
    M2 = M;
    M1(2:3,1) = [s-s0; t];
    M2(2:3,2) = [s-s0; t];
    A  = det(M);
    A1 = det(M1);
    A2 = det(M2);
    %disp('nodes');
    %disp(nodes);
    L1 = A1/A;
    L2 = A2/A;
    %disp('L1');
    %disp(L1);
    %disp('L2');
    %disp(L2);
    betaHat = [ L1 L2 (1-L1-L2)]*bHat(nodes);
                        
    
