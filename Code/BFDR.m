function [ psi, pst ] = BFDR( Beta, delt, alf )
%% BFDR For a specified delta (delt) and alpha (alf) computes Bayesian False
%       Discovery Rate using MCMC samples Beta
%
%   Created: 3/11/2014
%   By: Mark John Meyer
    B       = size(Beta,1); % number of samples
    R       = size(Beta,2); % number of coefficients
    pst     = NaN(R,1);
    lFDR    = NaN(R,1);
    for r = 1:R
        Betar       = abs(Beta(:,r));
        pst(r)      = (1/B)*size(find(Betar > delt),1);
        if pst(r) == 1;
            pst(r)  = 1 - (2*B)^(-1);
        end;
        lFDR(r)     = 1 - pst(r);
    end
    if sum(pst) == 0
        psi     = zeros(size(pst));
        fprintf('\n No p(s,t) > delta, psi set to zero \n');
    else
        pr      = sort(pst,'descend');
        rstar   = cumsum(1-pr)./linspace(1,R,R)';
        gam     = find(rstar <= alf, 1, 'last' );
        if isempty(gam)
            phi = 1;
        else
            phi = pr(gam);
        end
        psi     = pst >= phi;
    end
end

