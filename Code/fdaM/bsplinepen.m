function penaltymat = bsplinepen(basisobj, Lfdobj, sparsewrd)
%BSPLINEPEN computes the bspline penalty matrix for penalty LFDOBJ.
%  Arguments:
%  BASISOBJ  ... a basis object
%  LFDOBJ    ... a linear differential operator object.  
%                Default int2Lfd(2)
%  SPARSEWRD ... if 1, return penaltymatrix in sparse storage mode.

%  Last modified:  20 January 2003

%  check BASISOBJ

if ~isa_basis(basisobj)
    error('BASISOBJ is not a basis object.');
end

%  check basis type

type = getbasistype(basisobj);
if ~strcmp(type, 'bspline')
    error('basisobj not of type bspline');
end

%  set up default value for SPARSEWRD

if nargin < 3, sparsewrd = 0; end

%  set up default linear differential operator

if nargin < 2, 
    range = getbasisrange(basisobj);
    Lfdobj = int2Lfd(2); 
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);
  
%  get basis information

nbasis    = getnbasis(basisobj);
params    = getbasispar(basisobj);
norder    = nbasis - length(params);
range     = getbasisrange(basisobj);
breaks    = [range(1),params,range(2)];  %  break points
nbreaks   = length(breaks);
ninterval = nbreaks - 1;    %  number of intervals defined by breaks
if length(breaks) < 2
    error('The length of argument breaks is less than 2.');
end
if min(diff(breaks)) <= 0
    error('Break points are not strictly increasing.');
end

%  get highest order of derivative and check

nderiv = getnderiv(Lfdobj);
if nderiv < 0, error('NDERIV is negative'); end
if nderiv >= norder
    error(['Derivative of order', num2str(nderiv),                   ...
           'cannot be taken for B-spline of order', num2str(norder), ...
           'Probable cause is a value of the NBASIS argument in',    ...
           ' function FD that is too small.']);
end

%  special case where LFD is D^NDERIV and NDERIV = NORDER - 1

if isinteger(Lfdobj) & nderiv == norder - 1
    %  special case of nderiv = norder - 1
    halfseq    = (breaks(2:nbreaks) + breaks(1:(nbreaks-1)))./2;
    halfmat    = bsplineM(halfseq, breaks, norder, nderiv);
    brwidth    = diff(breaks);
    penaltymat = sparse(halfmat' * diag(brwidth) * halfmat);
end

%  if LFD is D^NDERIV, use exact computation

if isinteger(Lfdobj)
  
    %  Set up the knot sequence

    knots = [range(1)*ones(1,norder), ...
             breaks(2:(nbreaks-1)),          ...
             range(2)*ones(1,norder)]; 

    % Construct  the piecewise polynomial representation
    
    polyorder = norder - nderiv;
    ndegree   = polyorder - 1;
    prodorder = 2*ndegree + 1;   % order of product
    polycoef  = zeros(ninterval, polyorder, norder); 
    indxdown  = norder:-1:nderiv+1;
    for i = 1:nbasis 
        %  compute polynomial representation of B(i,norder,t)(x)
        [Coeff,index] = ppBspline(knots(i:i+norder));
        nrowcoef = size(Coeff,1);
        onescoef = ones(nrowcoef,1);
        % convert the index of the breaks in t to the index in the
        % variable 'breaks'
        index = index + i - norder;
        CoeffD = Coeff(:,1:polyorder);
        if nderiv > 0
            for ideriv=1:nderiv
                fac = indxdown - ideriv;
                CoeffD = (onescoef*fac).*CoeffD;
            end
        end
        % add the polynomial representation of B(i,norder,t)(x) to f
        if i >= norder, k = norder;  else k = i;       end
        if i <= norder, m = i;       else m = norder;  end
        for j=1:nrowcoef
            polycoef(i-k+j,:,m-j+1) = CoeffD(j,:);
        end
    end
    
    % Compute the scalar products 
    
    prodmat = zeros(nbasis);
    convmat = zeros(norder, norder, prodorder);
    for in = 1:ninterval
        %  get the coefficients for the polynomials for this interval
        Coeff = squeeze(polycoef(in,:,:));
        %  compute the coefficients for the products
        for i=0:ndegree-1
            ind = (0:i) + 1;
            convmat(:,:,i+1        ) = ...
                Coeff(ind,          :)'*Coeff(i-ind+2,      :);
            convmat(:,:,prodorder-i) = ...
                Coeff(ndegree-ind+2,:)'*Coeff(ndegree-i+ind,:);
        end
        ind = (0:ndegree)+1;
        convmat(:,:,ndegree+1) = Coeff(ind,:)'*Coeff(ndegree-ind+2,:);
%         convmat
%         convmat = polyprod(Coeff, Coeff);
%         convmat
        %  compute the coefficients of the integral
        delta    = breaks(in+1) - breaks(in);
        power    = delta;
        prodmati = zeros(norder);
        for i=1:prodorder
            prodmati = prodmati + ...
                power.*squeeze(convmat(:,:,prodorder-i+1))./i;
            power = power*delta;
        end
        % add the integral to s
        index = in:in+norder-1;
        prodmat(index,index) = prodmat(index,index) + prodmati; 
    end
    
    if sparsewrd
        penaltymat = sparse(prodmat);
    else
        penaltymat = prodmat;
    end
    
else
    
    %  LFDOBJ is not D^NDERIV, use approximate integration by calling
    %  function INPROD().
    
    prodmat = inprod(basisobj, basisobj, Lfdobj, Lfdobj);
    
    if sparsewrd
        penaltymat = sparse(prodmat);
    else
        penaltymat = prodmat;
    end
    
end
