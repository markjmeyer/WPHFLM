function dy = derivs(tnow, y) 
% DERIVS sets up the 1st order system corresponding to   
%   linear differential operator defined by wfd.

%  last modified 7 January 2003

global wfd;   
w  = eval_fd(wfd, tnow);
m  = length(w);
wmat = zeros(m, m);
wmat(1:(m-1),2:m) = eye(m-1);
wmat(m,:) = -w;
dy = wmat * y;

