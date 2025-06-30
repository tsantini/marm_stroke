function bfuncOut = basisFuncHarm(x,y,z,n,subInds,weights)
    % Outputs basis function for index n of spherical harmonics
    % n = index of phs_spha in sampHighOrder.m
    if (nargin<6)
        weights = [];
    end
    if (nargin > 4) && ~isempty(subInds)
        x = x(subInds(1):subInds(2));
        y = y(subInds(1):subInds(2));
        z = z(subInds(1):subInds(2));
    end
    if ~isempty(weights)
        % Apply scalings. Assume that weights were determined in units of
        % decimeter. This makes all relative k-coefficient values similar
        % to relative phase values at a distance of 1 dm from isocenter.
        x = x*10;
        y = y*10;
        z = z*10;
    end
    szx = size(x);
    bfuncOut = zeros([numel(x),numel(n)],'like',x);
    for n_ind = 1:numel(n)
        n_a = n(n_ind);
        switch n_a
            case 1
            bfunc = ones(size(x), 'like', x);
        case 2
            bfunc = x;
        case 3                
            bfunc = y;
        case 4                
            bfunc = z;
        case 5                
            bfunc = (x.*y);
        case 6               
            bfunc = (z.*y);
        case 7                
            bfunc = ((3.*z.^2 - (x.^2 + y.^2 + z.^2)));
        case 8                
            bfunc = (x.*z);
        case 9                
            bfunc = (x.^2 - y.^2);
        case 10               
            bfunc = (3.*y.*x.^2 - y.^3);
        case 11                
            bfunc = (x.*z.*y);
        case 12               
            bfunc = ((5.*z.^2 - (x.^2 + y.^2 + z.^2)).*y);
        case 13               
            bfunc = (5*z.^3 - 3.*z.*(x.^2 + y.^2 + z.^2));
        case 14                
            bfunc = ((5.*z.^2 - (x.^2 + y.^2 + z.^2)) .* x);
        case 15                
            bfunc = (x.^2.*z - y.^2.*z);
        case 16                
            bfunc = (x.^3 - 3.*x.*y.^2);
        case 17                
            bfunc = (x.*y).*(x.^2 - y.^2);
        case 18                
            bfunc = (y.*z).*(3.*x.^2 - y.^2);
        case 19               
            bfunc = (x.*y).*(7.*z.^2 - (x.^2 + y.^2 + z.^2));
        case 20                
            bfunc = y.*(7*z.^3 - 3.*z.*(x.^2 + y.^2 + z.^2));
        case 21               
            bfunc = 35.*z.^4 - 30.*z.^2.*(x.^2 + y.^2 + z.^2) + 3.*((x.^2 + y.^2 + z.^2).^2);
        case 22               
            bfunc = x.*(7*z.^3 - 3.*z.*(x.^2 + y.^2 + z.^2));
        case 23                
            bfunc = (x.^2 - y.^2).*(7.*z.^2 - (x.^2 + y.^2 + z.^2));
        case 24                
            bfunc = x.*(x.^2 - 3.*y.^2).*z;
        case 25                
            bfunc = x.^2.*(x.^2 - 3.*y.^2) - y.^2.*(3.*x.^2 - y.^2);
        case 26                
            bfunc = y.*(5.*x.^4 - 10.*x.^2.*y.^2 + y.^4);
        case 27                
            bfunc = x.*y.*z.*(x.^2 - y.^2);
        case 28                
            bfunc = y.*((x.^2 + y.^2 + z.^2) - 9.*z.^2).*(y.^2 - 3.*x.^2);
        case 29                
            bfunc = x.*y.*z.*(3.*z.^2 - (x.^2 + y.^2 + z.^2));
        case 30                
            bfunc = y.*((x.^2 + y.^2 + z.^2).^2 - 14.*(x.^2 + y.^2 + z.^2).*z.^2 + 21.*z.^4);
        case 31               
            bfunc = 63.*z.^5 - 70.*z.^3.*(x.^2 + y.^2 + z.^2) + 15.*z.*(x.^2 + y.^2 + z.^2).^2;
        case 32                
            bfunc = x.*((x.^2 + y.^2 + z.^2).^2 - 14.*(x.^2 + y.^2 + z.^2).*z.^2 + 21.*z.^4);
        case 33               
            bfunc = z.*(3.*z.^2 - (x.^2 + y.^2 + z.^2)).*(x.^2 - y.^2);
        case 34               
            bfunc = x.*((x.^2 + y.^2 + z.^2) - 9.*z.^2).*(3.*y.^2 - x.^2);
        case 35                
            bfunc = z.*(x.^4 - 6.*x.^2.*y.^2 + y.^4);
        case 36                
            bfunc = x.^5 - 10.*x.^3.*y.^2 + 5.*x.*y.^4;  
        otherwise
            error('Basis function undefined.')
        end
        bfuncOut(:,n_ind) = bfunc(:);
    end
    if ~isempty(weights)
        bfuncOut = sum(bfuncOut.*weights(:).',2);
    else
        bfuncOut = squeeze(reshape(bfuncOut,[szx,size(n)]));
    end
end