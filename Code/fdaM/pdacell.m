function [bfdcell, resfdcell, afdcell] = ...
    pdacell(xfdcell, ufdcell, awtcell, bwtcell, norder, nfine)
%PDA computes the basis function expansions of the
%  estimates of the coefficient functions a_k(t) and b_j(t) 
%  in the possibly nonhomogeneous linear differential operator
%
%    Lx(t) = a_1(t)u_1(t) + ... + a_k(t)u_K(t) + 
%       b_0(t)x(t) + b_1(t)Dx(t) + ... + b_{m-1}D^{m-1}x(t) + D^m x(t)
%
%  of order m = NORDER that minimizes in a least squares sense the residual
%  functions f(t) = Lx(t).  
%
%  If NORDER = 0, PDACELL fits the varying coefficient or pointwise
%  linear model using the functions x(t) as dependent variables and
%  the forcing functions u(t) as independent variables.  In this case,
%  there must be at least one forcing function.
%
%  The functions x(t) are in functional data object XFDOBJ.
%  The forcing functions u_k(t) are in functional data object UFDOBJ.
%  The coefficient functions for u_k(t) and x(t) are expanded in terms of the
%  basis functions specified in AWTCELL and BWTCELL, respectively.
%
%  The functions u_k(t) and x(t) are assumed to be vector valued 
%    of dimension J. 
%  That is, the differential equation can be a system of J equations rather
%    than a single equation.   
%  Each coefficient function b_j(t) is matrix valued, with a column
%    for each scalar function in the system.  

%  Arguments:
%  XFDCELL   ...  cell array of functional data objects for the functions
%                 whose derivatives define the DIFE
%                 dimensions are J and 1
%  UFDCELL   ...  cell array of independent variables or u-variables
%                 dimensions are J and K
%  AWTCELL   ...  cell array of weight function specs for u-variables
%                 dimensions are J and K
%  BWTCELL   ...  cell array of weight function specs for functions x
%                 dimensions are J, J and NORDER
%  NORDER    ...  order of the linear differential operator, that is, 
%                 the order of the highest derivative.
%  NFINE     ...  number of sampling points for numerical integration

%  The value in each cell of XFDCELL, UFDCELL, AFDCELL and BFDCELL is a
%      scalar FD object

%  The value in each cell of AWTCELL and BWTCELL is a struct containing:
%       FD object containing to be used as the fixed value if not estimated
%       ESTIMATE: 1 if weight is to be estimated, 0 if not.
%  We are not bothering with smoothing at this point.

%  Returns:
%  BFDCELL   ...  cell array of weights for x functions
%                 dimensions are J, J and NORDER
%  RESFDCELL ...  FD object for residual functions.
%  AFDCELL   ...  cell array of weights for u-variables
%                 dimension J and K

%  last modified 28 January 2003

if nargin < 6, nfine   = 101; end
if nargin < 5, norder  = 1;   end
if nargin < 4, bwtcell = {};  end
if nargin < 3, awtcell = {};  end
if nargin < 2, ufdcell = {};  end

if nargin < 2
    error('There are less than four arguments.');
end

norder = floor(norder);
if norder < 0, error('NORDER is negative.');  end
nordp1 = norder + 1;

%  check dimensions of cells

nvar = size(xfdcell,1);
if size(xfdcell,2) ~= 1
    error('XFDCELL has more than one column.');
end

%  check the dimensions of UFDCELL and AWTCELL if there are
%  
if isempty(ufdcell) | isempty(awtcell)
    nu = 0;
    afdcell = {};
else
    nu = size(ufdcell,2);
    if size(ufdcell,1) ~= nvar
        error(['The number of rows of UFDCELL]', ...
               ' does not match that of XFDCELL.']);
    end
    if any(size(awtcell) ~= [nvar, nu])
        error('The dimensions of AWTCELL are incorrect.');
    end
end

%  check to see if there is anything to estimate

if norder == 0 & nu == 0
    error('There are no coefficient functions to estimate.');
end

%  check the dimensions of BWTCELL

if norder == 0
    bfdcell = {};
end
if norder == 1
    if any(size(bwtcell) ~= [nvar, nvar])
        error('The dimensions of BWTCELL are incorrect.');
    end
end
if norder > 1
    if any(size(bwtcell) ~= [nvar, nvar, norder])
        error('The dimensions of BWTCELL are incorrect.');
    end
end

%  check XFDCELL and extract NCURVE and XRANGE

for ivar=1:nvar
    if ~isa_fd(xfdcell{ivar})
        error(['XFDCELL{',num2str(ivar), ...
                '} is not a functional data object.']);
    end 
    if ivar == 1
        xrange     = getbasisrange(getbasis(xfdcell{ivar}));
        ncurve     = size(getcoef(xfdcell{ivar}),2);
        bfdnames   = getnames(xfdcell{ivar});
        resfdnames = bfdnames;
    else
        if any(getbasisrange(getbasis(xfdcell{ivar})) ~= xrange)
            error('Ranges are incompatible for XFDCELL.');
        end
        if size(getcoef(xfdcell{ivar}),2) ~= ncurve
            error('Number of curves is incompatible for XFDCELL.');
        end
    end
end

nbasmax = 0;  %  This will be the maximum number of basis functions

%  check UFDCELL and extract URANGE

if nu > 0
    for ivar=1:nvar
        for iu=1:nu
            if ~isa_fd(ufdcell{ivar,iu})
                error(['UFDCELL{',num2str(ivar),',',num2str(iu), ...
                        '} is not a functional data object.']);
            end  
            if ivar == 1 & iu == 1
                urange = getbasisrange(getbasis(ufdcell{ivar,iu}));
                afdnames = getnames(ufdcell{ivar,iu});
            else
                if any(getbasisrange(getbasis(ufdcell{ivar,iu})) ~= urange)
                    error('Ranges are incompatible for UFDCELL.');
                end
            end
        end
    end

    %  check AWTCELL and extract the max. no. basis fns.
    
    for ivar=1:nvar
        for iu=1:nu
            awtstructi = awtcell{ivar,iu};
            afdi       = awtstructi.fd;
            if ~isa_fd(afdi)
                error(['AFDI is not a functional data object.']);
            end  
            basisi = getbasis(afdi);
            if any(getbasisrange(basisi) ~= urange)
                error('Ranges are incompatible for AWTCELL.');
            end
            nbasi   = getnbasis(basisi);
            nbasmax = max([nbasmax,nbasi]);
        end
    end

end

%  check BWTCELL

if norder > 0
    for ivar1=1:nvar
        for ivar2=1:nvar
            for j=1:norder
                if norder == 1
                    bstruct12 = bwtcell{ivar1,ivar2};
                else
                    bstruct12 = bwtcell{ivar1,ivar2,j};
                end
                if ~isa_fd(bstruct12.fd)
                    error(['BWTCELL{',num2str(ivar1), ', ', ...
                            num2str(ivar2), ', ', ...
                            num2str(iu),      ...
                            '} is not a functional data object.']);
                end
                basisi = getbasis(bstruct12.fd);
                if any(getbasisrange(basisi) ~= xrange)
                    error('Ranges are incompatible for BWTCELL.');
                end
                nbasi   = getnbasis(basisi);
                nbasmax = max([nbasmax,nbasi]);
            end  
        end
    end
end

%  At this point we assume that the ranges for XFDCELL and UFDCELL 
%  are the same, but this will be changed later to allow for lags.

if nu > 0
    if any(xrange ~= urange)
        error('Ranges for XFDCELL and UFDCELL are not compatible.');
    end
end

%  set up sampling values to be used in numerical integration
%    and set up matrix of basis values.  The number of sampling
%  NFINE is here set to a usually workable value if too small.

if nfine < 5*nbasmax, nfine = 5*nbasmax;  end

deltax = (xrange(2)-xrange(1))/(nfine-1);
tx     = (xrange(1):deltax:xrange(2))';

if nu > 0
    deltau = (urange(2)-urange(1))/(nfine-1);
    tu     = urange(1):deltau:urange(2);
end

%  set up  YARRAY to hold values of x functions and their derivatives

yarray = zeros([nfine,ncurve,nvar,nordp1]);
for ivar=1:nvar
    for j=1:nordp1
        yarray(:,:,ivar,j) = eval_fd(xfdcell{ivar}, tx, j-1);
    end
end

%  set up  UARRAY to hold values of u functions

if nu > 0
    uarray = zeros([nfine,nu]);
    for iu=1:nu
        uarray(:,iu) = eval_fd(ufdcell{ivar,iu}, tu);
    end
end

%  set up array YPROD to hold mean of products of values in YARRAY

mmat  = m2ij(nvar,nordp1);
yprod = zeros([nfine,nvar,nordp1,nvar,nordp1]);
for m1=1:nvar*nordp1
    i1 = mmat(m1,1);
    j1 = mmat(m1,2);
    for m2=1:m1;
        i2 = mmat(m2,1);
        j2 = mmat(m2,2);
        if ncurve == 1
            yprodval = squeeze(yarray(:,1,i1,j1)).*squeeze(yarray(:,1,i2,j2));
        else
            yprodval = mean(squeeze(yarray(:,:,i1,j1)).* ...
                            squeeze(yarray(:,:,i2,j2)),2);
        end
        yprod(:,i1,j1,i2,j2) = yprodval;
        yprod(:,i2,j2,i1,j1) = yprodval;
    end
end

%  set up array YUPROD to hold mean of u-variables u times 
%    x functions and their derivatives

onesncurve = ones(1,ncurve);
if nu > 0
    yuprod = zeros(nfine, nvar, nu, nordp1);
    for iu=1:nu
        for i1=1:nvar
            for j1=1:nordp1
                if ncurve == 1
                    yuprodval = yarray(:,1,i1,j1).*uarray(:,iu);
                else
                    yuprodval = ...
                        mean(squeeze(yarray(:,:,i1,j1).* ...
                                    (uarray(:,iu)*onesncurve)),2);
                end
                yuprod(:,i1,iu,j1) = yuprodval;
            end
        end
    end
end

clear yarray

%  set up array UPROD to hold mean of products of u-variables u 

if nu > 0
    uprod = zeros(nfine, nu, nu);
    for iu=1:nu
        for ju=1:iu
            uprodval = uarray(:,iu).*uarray(:,ju);
            uprod(:,iu,ju) = uprodval;
            uprod(:,ju,iu) = uprodval;
        end
    end
    clear uarray
end

%  set up an index array and some arrays of 1's

mmat  = m2ij(nvar,norder);
onesn = ones(nfine,1);

%  set up array to hold coefficients for basis expansions

if nu > 0
    aarray = zeros(nfine,nu);  
else
    aarray = [];
end
if norder > 0
    barray = zeros(nfine,nvar,norder);  
else
    barray = [];
end


%  --------------  beginning of loop through variables  -------------------

for ivar=1:nvar
    
    %  get number of coefficients to be estimated for this equation
    
    % loop through u-variables
    neqns  = 0;
    for iu = 1:nu
        astructi = awtcell{ivar,iu};
        if astructi.estimate 
            neqns = neqns + getnbasis(getbasis(astructi.fd));
        end
    end
    % loop through x functions and their derivatives
    for m2=1:nvar*norder
        i2 = mmat(m2,1);
        j2 = mmat(m2,2);
        if norder == 1
            bstructij = bwtcell{ivar,i2};
        else
            bstructij = bwtcell{ivar,i2,j2};
        end
        if bstructij.estimate
            neqns = neqns + getnbasis(getbasis(bstructij.fd));
        end
    end
    if neqns < 1
        error('Number of equations to solve is not positive.');
    end

    %  set up coefficient array and right side array for linear equation
    
    cmat   = zeros(neqns, neqns);
    dmat   = zeros(neqns, 1);
    
    %  evaluate default weight functions for this variable
    
    for iu=1:nu
        astructi = awtcell{ivar,iu};
        aarray(:,iu) = eval_fd(astructi.fd, tu);
    end
    for i=1:nvar
        for j=1:norder
            if norder == 1
                bstructij = bwtcell{ivar,i};
            else
                bstructij = bwtcell{ivar,i,j};
            end
            barray(:,i,j) = eval_fd(bstructij.fd, tx);
        end
    end
    
    %  loop through equations, 
    %    corresponding to rows for CMAT and DMAT
    
    %  loop through equations for u-variables
    
    mi12 = 0;
    for iu1 = 1:nu
        astructi1   = awtcell{ivar,iu1};
        if astructi1.estimate
            abasisi1    = getbasis(astructi1.fd);
            abasismati1 = getbasismatrix(tu, abasisi1);
            mi11 = mi12 + 1;
            mi12 = mi12 + getnbasis(abasisi1);
            indexi1 = mi11:mi12;
            %  DMAT entry for u-variable
            weighti1 = yuprod(:,ivar,iu1,nordp1);
            dmat(indexi1) = ...
                trapzmat(abasismati1,onesn,deltax,weighti1);
            %  loop through weight functions to be estimated,
            %    corresponding to columns for CMAT
            %  begin with u-variables
            mi22 = 0;
            for iu2=1:nu
                astructi2   = awtcell{ivar,iu2};
                if astructi2.estimate
                    abasisi2    = getbasis(astructi2.fd);
                    abasismati2 = getbasismatrix(tu, abasisi2);
                    weighti2    = uprod(:,iu1,iu2);
                    Cprod  = ...
                        trapzmat(abasismati1, abasismati2, deltau, weighti2);
                    if astructi2.estimate
                        mi21 = mi22 + 1;
                        mi22 = mi22 + getnbasis(abasisi2);
                        indexi2 = mi21:mi22;
                        %  coefficient matrix CMAT entry
                        cmat(indexi1,indexi2) = Cprod;
                    else
                        dmat(indexi1) = dmat(indexi1) + ...
                            Cprod*aarray(:,iu2);
                    end
                end
            end
            %  remaining columns: 
            %    loop through u-variable -- x-derivative pairs
            mij22 = mi22;
            for m2=1:nvar*norder;
                i2 = mmat(m2,1);
                j2 = mmat(m2,2);
                if norder == 1
                    bstructij2   = bwtcell{ivar,i2};
                else
                    bstructij2   = bwtcell{ivar,i2,j2};
                end
                bbasisij2    = getbasis(bstructij2.fd);
                bbasismatij2 = getbasismatrix(tx, bbasisij2);
                weightij12   = yuprod(:,ivar,iu1,j2);
                Cprod = ...
                    trapzmat(abasismati1,bbasismatij2,deltax,weightij12);
                if bstructij2.estimate
                    mij21 = mij22 + 1;
                    mij22 = mij22 + getnbasis(bbasisij2);
                    indexij2  = mij21:mij22;
                    cmat(indexi1,indexij2) = Cprod;
                else
                    dmat(indexj) = dmat(indexj) + ...
                        Cprod*barray(:,i2,j2);
                end
            end
        end
    end
    
    %  loop through equations for x-derivatives
    
    mij12 = mi12;
    for m1=1:nvar*norder
        i1 = mmat(m1,1);
        j1 = mmat(m1,2);
        if norder == 1
            bstructij1 = bwtcell{ivar,i1};
        else
            bstructij1 = bwtcell{ivar,i1,j1};
        end
        if bstructij1.estimate
            bbasisij1    = getbasis(bstructij1.fd);
            bbasismatij1 = getbasismatrix(tx,bbasisij1);
            mij11 = mij12 + 1;
            mij12 = mij12 + getnbasis(bbasisij1);
            indexij1 = mij11:mij12;
            %  DMAT entry for u-variable -- x-derivative pair
            weightij1 = yprod(:,i1,j1,ivar,nordp1);
            dmat(indexij1) = ...
                trapzmat(bbasismatij1,onesn,deltax,weightij1);
            %  first columns of CMAT: u-variable entries
            mi22 = 0;
            for iu2=1:nu
                astructi2   = awtcell{ivar,iu2};
                if astructi2.estimate
                    abasisi2    = getbasis(astructi2.fd);
                    abasismati2 = getbasismatrix(tx, abasisi2);
                    weighti2    = yuprod(:,i1,iu2,j1);
                    Cprod = ...
                        trapzmat(bbasismatij1,abasismati2,deltax,weighti2);
                    if astructi2.estimate
                        mi21 = mi22 + 1;
                        mi22 = mi22 + getnbasis(abasisi2);
                        indexi2 = mi21:mi22;
                        cmat(indexij1,indexi2) = Cprod;
                    else
                        %dmat(indexij1) = dmat(indexij1) + ...
                        %    Cprod*aarray(:,iu2);
                    end
                end
            end
            %  remaining columns: x-derivative pairs
            mij22 = mi22;
            for m2=1:nvar*norder;
                i2 = mmat(m2,1);
                j2 = mmat(m2,2);
                if norder == 1
                    bstructij2   = bwtcell{ivar,i2};
                else
                    bstructij2   = bwtcell{ivar,i2,j2};
                end
                bbasisij2    = getbasis(bstructij2.fd);
                bbasismatij2 = getbasismatrix(tx, bbasisij2);
                weightij22   = yprod(:,i1,j1,i2,j2);
                Cprod = ...
                    trapzmat(bbasismatij1,bbasismatij2,deltax,weightij22);
                if bstructij2.estimate
                    mij21 = mij22 + 1;
                    mij22 = mij22 + getnbasis(bbasisij2);
                    indexij2 = mij21:mij22;
                    cmat(indexij1,indexij2) = Cprod;
                else
                    %dmat(indexij1) = dmat(indexij1) + ...
                    %    Cprod*barray(:,i2,j2);
                end
            end
        end
    end
    
  
    dvec = -cmat\dmat;
  
    %  set up u-function weight functions
    
    mi2 = 0;
    for iu=1:nu
        astructi = awtcell{ivar,iu};
        if astructi.estimate
            mi1 = mi2 + 1;
            mi2 = mi2 + getnbasis(getbasis(astructi.fd));
            indexi = mi1:mi2;
            afdcell{ivar,iu} = putcoef(astructi.fd, dvec(indexi));
        else
            afdcell{ivar,iu} = astructi.fd;
        end
    end
    
    %  set up X-function derivative weight functions
    
    mij2 = mi2;
    for m1=1:nvar*norder;
        i1 = mmat(m1,1);
        j1 = mmat(m1,2);
        bstructij = bwtcell{ivar,i1,j1};
        if bstructij.estimate
            mij1 = mij2 + 1;
            mij2 = mij2 + getnbasis(getbasis(bstructij.fd));
            indexij = mij1:mij2;
            bfdcell{ivar,i1,j1} = putcoef(bstructij.fd, dvec(indexij));
        else
            bfdcell{ivar,i1,j1} = bstructij.fd;
        end
    end
    
end

%  --------------  end of loop through variables  -------------------

%  set up residual cell RESFDCELL

resfdcell     = cell(nvar);
resfdnames{2} = 'Residual function';
resfdnames{3} = 'Residual function value';

resbasis = getbasis(xfdcell{1});
for ivar=1:nvar
    %  initialize with highest order derivative for this variable
    resmat  = eval_fd(xfdcell{ivar}, tx, norder);
    %  add contributions from weighted u-functions
    for iu=1:nu
        amati    = eval_fd(afdcell{ivar,iu}, tu);
        umati    = eval_fd(ufdcell{ivar,iu}, tu);
        resmat   = resmat + (amati.*umati)*onesncurve;
    end
    %  add contributions from weighted x-function derivatives
    for m1=1:nvar*norder;
        i1 = mmat(m1,1);
        j1 = mmat(m1,2);
        bmatij = eval_fd(bfdcell{ivar,i1,j1}, tx)*onesncurve;
        xmatij = eval_fd(xfdcell{i1},         tx, j1-1);
        resmat = resmat + bmatij.*xmatij;
    end
    %  set up the functional data object
    resfdi = data2fd(resmat, tx, resbasis);
    resfdi = putnames(resfdi, resfdnames);
    resfdcell{ivar} = resfdi;
end

%  --------------------------------------------------------------

function mmat = m2ij(nrow,ncol)
%M2IJ sets up a NROW*NCOL by 2 matrix of row-col indices associated
%  with a number of matrix entries row-wise
%  Example:  m2ij(2,3) produces
%     1     1
%     1     2
%     1     3
%     2     1
%     2     2
%     2     3
nval = nrow*ncol;
if nval > 0
    mmat = [reshape(ones(ncol,1)*(1:nrow), nval,1), ...
            reshape((1:ncol)'*ones(1,nrow),nval,1)];
else
    mmat = [];
end
