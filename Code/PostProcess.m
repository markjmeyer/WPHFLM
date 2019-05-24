function postout = ...
            PostProcess(MCMC_beta,MCMC_zeta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wpspecs)
    %% Outputs:
    %  thetahat,theta_025CI,theta_975CI,tauhat,pihat,Sigma,uhat,u_025CI,u_975CI,g_Ut,g_Ut025,g_Ut975,ghatns
    %% Eventual Inputs
    %  MCMC_beta,MCMC_zeta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wpspecs
    
    %% Test:
%     MCMC_beta           = res.MCMC_beta;
%     MCMC_zeta           = res.MCMC_zeta;
%     model               = res.model;
%     wpspecs             = res.wpspecs;
%     MCMC_tau            = res.MCMC_tau;
%     MCMC_pi             = res.MCMC_pi;
%     MCMC_flag_theta     = res.MCMC_flag_theta;
%     MCMC_theta          = res.MCMC_theta;
%     wpspecs.V           = res.wpspecs.wpspecsy.trackerJacker{1}.T;
%     wpspecs.T           = res.wpspecs.wpspecsx.trackerJacker{1}.T;

    %% function parameters
    p           = model.p;
    K           = wpspecs.K;
    V           = wpspecs.V; %% use V for fsize, may need to change if we explore data reduction in X
    T           = wpspecs.T;
    fsize       = size(model.Dx,2);
    w           = model.p - fsize; %% wsize
    wpspecsy    = wpspecs.wpspecsy;
    wpspecsx    = wpspecs.wpspecsx;
    keep        = model.keep; %% threshold size

    %%
    delt        = model.delt;
    alf         = model.alf;

    %% 
    postout.accept_rate_theta   = mean(MCMC_flag_theta);
    postout.tauhat              = reshape(mean(MCMC_tau),p,wpspecs.J);
    postout.pihat               = reshape(mean(MCMC_pi)',p,wpspecs.J);
    postout.alpha               = reshape(mean(MCMC_alpha)',p,K)';
    
    %% 
    B       = size(MCMC_beta,1);
    Beta    = NaN(B,V*T);
    Zeta    = NaN(B,w*T);

    %%
    for i = 1:B
        %% extract coefficients of interest %%
        Betam       = reshape(MCMC_beta(i,:),fsize,K);
        Zm          = reshape(MCMC_zeta(i,:)',K,w)';
        
        %% add back thresholded coefs %%
        Betat           = zeros(K,K);
        Betat(keep,:)   = Betam;
        
        %% project Beta and Zeta back into Y space %%
        BZ          = [Betat; Zm];
        BZidwpt     = idwpt_rows(BZ,wpspecsy);
        
        yidwpt      = BZidwpt(1:K,:);
        Zidwpt      = BZidwpt((K+1):end,:);
                                    
        %% project Beta back into X space %%
        Bidwpt      = idwpt_rows(yidwpt',wpspecsx)';
        
        %% update Beta and Zeta %%
        Beta(i,:)   = reshape(Bidwpt,1,V*T);
        Zeta(i,:)   = reshape(Zidwpt,1,w*T);

        if mod(i, 10) == 0
            fprintf('.')
        end
        if mod(i,100) == 0
            fprintf('\n Done processing %d MCMC samples \n',i);
        end
    end
    
    %% post-process bhat %%
    bhatpre     = reshape(mean(Beta),V,T);
    bhat        = zeros(V,T);
    for i = 1:T
        for j = i:V
            if j < i+round(T*wpspecs.perlagback)
                bhat(i,j) = bhatpre(i,j);
            end
        end
    end
    
    %% post-process Q025 bhat %%
    Q025            = reshape(quantile(Beta,0.025),V,T);
    Q025_bhat       = zeros(V,T);
    for i = 1:T
        for j = i:V
            if j < i+round(T*wpspecs.perlagback)
                Q025_bhat(i,j) = Q025(i,j);
            end
        end
    end
    
    %% post-process Q975 bhat %%
    Q975            = reshape(quantile(Beta,0.975),V,T);
    Q975_bhat       = zeros(V,T);
    for i = 1:T
        for j = i:V
            if j < i+round(T*wpspecs.perlagback)
                Q975_bhat(i,j) = Q975(i,j);
            end
        end
    end
    
    %% average and find credible intervals %%
    postout.beta        = Beta;
    postout.zeta        = Zeta;
    postout.Q025_bhat   = Q025_bhat;
    postout.Q975_bhat   = Q975_bhat;
    postout.bhat        = bhat;
    postout.bsd         = reshape(std(Beta),V,T);
    postout.Q025_zhat   = reshape(quantile(Zeta,0.025),T,w)';
    postout.Q975_zhat   = reshape(quantile(Zeta,0.975),T,w)';
    postout.zhat        = reshape(mean(Zeta),T,w)';
    fprintf('\n Done averaging and finding CIs.\n \n');
    
    %% coefficients for inference %%
    bh      = zeros(V,T);
    for i = 1:T
        for j = i:V
            if j < i+round(T*wpspecs.perlagback)
                bh(i,j) = 1;
            end
        end
    end
    bFlag           = reshape(bh,1,V*T);
    bINF            = Beta(:,bFlag == 1);
    postout.bINF    = bINF;
    
    %% implement FDR %%
    R       = sum(bFlag);
    pst     = NaN(R,1);
    lFDR    = NaN(R,1);
    for r = 1:R
        Betar       = abs(bINF(:,r));
        pst(r)      = (1/B)*size(find(Betar > delt),1);
        if pst(r) == 1
            pst(r)  = 1 - (2*B)^(-1);
        end
        lFDR(r)     = 1 - pst(r);
    end
    if sum(pst) == 0
        error('No coef > delta');
    end
    pr      = sort(pst,'descend');
    rstar   = cumsum(1-pr)./linspace(1,R,R)';
    gam     = find(rstar <= alf, 1, 'last' );
    if isempty(gam)
        phi = 1;
    else
        phi = pr(gam);
    end
    psi     = pst >= phi;
    
    %% insert psi back in to full vector %%
    psif                = zeros(size(bFlag'));
    pstf                = zeros(size(bFlag'));
    psif(bFlag == 1)    = psi;
    pstf(bFlag == 1)    = pst;

    postout.psi     = reshape(psif,V,T);
    postout.pst     = reshape(pstf,V,T);
    fprintf('\n Done calculating FDR.\n \n');

    %% implement SimBaS %%
    [SBSb, upper_CIb, lower_CIb]	= jointband_sbs(bINF,alf);
    SBS                             = zeros(size(bFlag));
    SBS(bFlag == 1)                 = SBSb;
    upper_CI                        = zeros(size(bFlag));
    upper_CI(bFlag == 1)            = upper_CIb;
    lower_CI                        = zeros(size(bFlag));
    lower_CI(bFlag == 1)            = lower_CIb;

    postout.MAPs                    = reshape(SBS,V,T);
    postout.UMAPs                   = reshape(upper_CI,V,T);
    postout.LMAPs                   = reshape(lower_CI,V,T);
    fprintf('\n Done calculating MAPs.\n \n');

    
    %% output thetas %%
    thetahat                = reshape(mean(MCMC_theta),size(theta,1),size(theta,2));
    postout.thetahat        = thetahat;
    postout.theta_025CI     = reshape(quantile(MCMC_theta,0.05),size(theta,1),size(theta,2));
    postout.theta_975CI     = reshape(quantile(MCMC_theta,0.95),size(theta,1),size(theta,2));
    fprintf('\n Done with regularization parameters.\n \n');
    
