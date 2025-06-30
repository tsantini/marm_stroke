function [x, resSqAll, mseAll, xnormAll, xdiffAll, stopThresh, status] = cgne(Ain,bin,x0,NitMax,opt,gam,Rin,W)
    % Use conjugate gradient method to solve Ax = b in a least squares sense. 
    %    i.e., "conjugate gradient on the normal equations", cgne
    %
    % x = cgne(A,b,x0,maxIterations,options,gamma,R,W)
    %
    %   A must have a transpose that can be evaluated using A'*b, or as a function with A(b,'transp')
    %   A operates on x via A*x or A(x,'notransp')
    %
    %   gamma, R (optional): specify to perform regularization using matrix R and tuning factor gamma
    %           solves argmin_x ||sqrt(W)(Ax - b)||^2_2 + gamma||Rx||^2_2
    %       Can provide these as cell arrays to have multiple
    %       regularization fcns.
    %
    %   W (optional): allows weighting data consistency cost.
    %       e.g., density compensation based weighting, Magn Reson Med. 2018 May; 79(5): 2685–2692
    %
    % (c) Corey Baron 2020
    %
    
    % Set options
    if nargin<4 || isempty(NitMax)
        % Maximum number of iterations allowed
        NitMax = 1000;
    end
    if nargin<5 || ~isfield(opt,'resThresh')
        % Iterations stop when squared residual is less than this
        opt.resThresh = 1e-10; 
    end
    if nargin<5 || ~isfield(opt,'resChangeThresh')
        % Iterations stop when fractional difference between residuals is less than this
        opt.resChangeThresh = []; 
    end
    if nargin<5 || ~isfield(opt,'xChangeThresh')
        % Iterations stop when square of l2 norm of difference of x
        % between iterations is less than this.
        opt.xChangeThresh = 1e-10; 
    end
    if nargin<5 || ~isfield(opt,'expResN')
        % Set to non-zero to explicitely compute residual every expResN iterations. 
        %   also restarts conj grad iterations by setting beta to 0
        % Can reduce accumulation of roundoff errors
        opt.expResN = 100; 
    end
    if nargin<5 || ~isfield(opt,'nseThreshFact')
        % Set factor to increase estimated noise by for stopping iterations
        % based on the discrepancy principle. 
        % For stopping criterion, see a great review here: 
        % https://doi-org.proxy1.lib.uwo.ca/10.3846/1392-6292.2007.12.61-70
        opt.nseThreshFact = 1.0; 
    end
    if nargin<5 || ~isfield(opt,'nseExtraIt')
        % Stopping on noise threshold seems to often stop a little too
        % early, so we add some extra iterations.
        opt.nseExtraIt = 4; 
    end
    if nargin<5 || ~isfield(opt,'stopOnResInc')
        % Stop on nth iteration where ||res|| increases, where n is the
        % value that opt.stopOnResInc is set to. This works well for single
        % shot spiral MRI. 
        % Similar to methods used for CT, https://doi.org/10.1016/S0096-3003(98)10007-3
        opt.stopOnResInc = 0; 
    end
    if nargin<5 || ~isfield(opt,'resIncThresh')
        % fractional threshold for stopOnResInc. Should be close to 1
        opt.resIncThresh = 1.0; 
    end
    if nargin<5 || ~isfield(opt,'stopOnXdifInc')
        % Stop on nth iteration where ||x_{k+1}-x_k|| increases, where n is the
        % value that opt.stopOnXdifInc is set to. 
        % Similar to methods used for CT, https://doi.org/10.1016/S0096-3003(98)10007-3
        opt.stopOnXdifInc = 0; 
    end
    if nargin<5 || ~isfield(opt,'noiseVar')
        % Variance of the noise expected in the input bin
        opt.noiseVar = []; 
    end
    if nargin<5 || ~isfield(opt,'plotting')
        % Show plots of progress
        opt.plotting = 0; 
    end
    if nargin<5 || ~isfield(opt,'plottingReshape')
        % Reshape x for plotting images. Only 2D matrix allowed.
        opt.plottingReshape = []; 
    end
    if nargin<5 || ~isfield(opt,'gtruth')
        % Useful for simulations. Allows computation of mse per iteration
        opt.gtruth = []; 
    end
    
    if nargin<6 || (~isempty(gam) && ~iscell(gam) && all(gam == 0))
        gam = [];
    end
    if nargin<7 || isempty(gam)
        Rin = [];
    end
    if nargin<8 || isempty(W)
        W = [];
    end

    findxnorm = 0;
    findxdiff = 0;
    if (nargout > 3) || ~isempty(opt.xChangeThresh)
        findxnorm = 1;
    end
    if (nargout > 4) || ~isempty(opt.xChangeThresh) || (opt.stopOnXdifInc>0) || opt.plotting
        findxdiff = 1;
    end
    
    % Account for different ways of supplying A and R
    if isa(Ain,'function_handle')
        A = Ain;
    else
        A = @(x,transp) Asub(x,transp,Ain);
    end
    if ~isempty(Rin)
        if iscell(Rin)
            R = cell(size(Rin));
            for nR = 1:length(Rin)
                if isa(Rin{nR},'function_handle')
                    R(nR) = Rin(nR);
                else
                    R{nR} = @(x,transp) Asub(x,transp,Rin{nR});
                end
            end
        else
            if isa(Rin,'function_handle')
                R = Rin;
            else
                R = @(x,transp) Asub(x,transp,Rin);
            end
        end
    end
    
    % Get RHS for normal equations
    if ~isempty(W)
        b = A(W.*bin,'transp');
    else
        b = A(bin,'transp');
    end

    % Set default starting guess
    if nargin<3 || isempty(x0)
        x0 = b;
    end

    % Prep for stopping based on measured noise
    stopThresh = [];
    if (~isempty(opt.noiseVar) && opt.noiseVar>0)
        % Create a synthetic noise vector
        if ~isreal(bin)
            % The scaling by 0.5 presumes the noise was estimated from
            % complex data
            nseVecIn = sqrt(opt.noiseVar/2)*(randn(size(bin)) + 1i*randn(size(bin)));
        else
            nseVecIn = sqrt(opt.noiseVar)*randn(size(bin));
        end
        % Pass the noise through A' to find the expected noise in b
        nse = A(nseVecIn,'transp');
        % Define the stopping criteria
        stopThresh = opt.nseThreshFact * nse(:)' * nse(:);
    end
    
    % Initialize first iteration
    xold = x0;
    if ~isempty(W)
        res = b - A(W.*A(xold, 'notransp'), 'transp');
    else
        res = b - A(A(xold, 'notransp'), 'transp');
    end
    if ~isempty(gam)
        if iscell(R)
            for nR = 1:length(R)
                res = res - gam{nR}.*R{nR}(R{nR}(xold, 'notransp'), 'transp');
            end
        else
            res = res - gam.*R(R(xold, 'notransp'), 'transp');
        end
    end
    p = res;
    res_sq_old = res(:)' * res(:);
    nit = 0;
    if opt.plotting
        nf = figure;
        Ncol = ceil(sqrt(NitMax));
        Nrow = ceil(NitMax/Ncol);
        if ~isempty(opt.gtruth)
            nf2 = figure;
        end
    end
    
    % Run iterations
    resSqAll = zeros(NitMax+1,1);
    resSqAll(1) = res_sq_old;
    xnormAll = zeros(NitMax+1,1);
    if findxnorm
        xnormAll(1) = x0(:)'*x0(:);
    end
    xdiffAll = zeros(NitMax+1,1);
    if findxdiff
        xdiffAll(1) = xnormAll(1);
    end
    mseAll = zeros(NitMax+1,1);
    if ~isempty(opt.gtruth)
        tmp = opt.gtruth(:)-x0(:);
        mseAll(1) = tmp'*tmp;
    end
    finished = 0;
    resIncsTotal = 0;
    xdifIncsTotal = 0;
    lastxdifInc = 0; % consider sequential climbing xdiff as a single increase
    status = 0;
    while ~finished
        if ~isempty(W)
            Ap = A(W.*A(p, 'notransp'), 'transp');
        else
            Ap = A(A(p, 'notransp'), 'transp');
        end
        if ~isempty(gam)
            if iscell(R)
                for nR = 1:length(R)
                    Ap = Ap + gam{nR}.*R{nR}(R{nR}(p, 'notransp'), 'transp');
                end
            else
                Ap = Ap + gam.*R(R(p, 'notransp'), 'transp');
            end
        end
        alph = res_sq_old / (p(:)' * Ap(:));
        x = xold + alph*p;
        if mod(nit,opt.expResN) == 0
            % Avoid accumulation of rounding errors when doing many iterations
            if ~isempty(W)
                res = b - A(W.*A(x, 'notransp'), 'transp');
            else
                res = b - A(A(x, 'notransp'), 'transp');
            end
            if ~isempty(gam)
                if iscell(R)
                    for nR = 1:length(R)
                        res = res - gam{nR}.*R{nR}(R{nR}(x, 'notransp'), 'transp');
                    end
                else
                    res = res - gam.*R(R(x, 'notransp'), 'transp');
                end
            end
            betaFact = 0;
        else
            res = res - alph*Ap;
            betaFact = 1;
        end
        % Tracking progress
        res_sq_new = res(:)' * res(:);
        resSqAll(nit+2) = res_sq_new; 
        if findxnorm
            xnormAll(nit+2) = x(:)'*x(:);
        end
        if findxdiff
            xdiff = x-xold;
            xdiffAll(nit+2) = xdiff(:)'*xdiff(:);
        end
        if ~isempty(opt.gtruth)
            tmp = opt.gtruth(:)-x(:);
            mseAll(nit+2) = tmp'*tmp;
        end
        % Plotting
        if opt.plotting
            figure(nf);
            subplot(Nrow,Ncol,nit+1);
            tmp = abs(xdiff);
            if ~isempty(opt.plottingReshape)
                tmp = reshape(tmp, opt.plottingReshape);
                imagesc(tmp);
            else
                plot(tmp);
            end
            hold('all')
            if ~isempty(opt.gtruth) && ~isempty(opt.plottingReshape)
                figure(nf2)
                subplot(Nrow,Ncol,nit+1);
                imagesc(reshape(abs(tmp), opt.plottingReshape));
                colorbar;
                hold('all')
            end
        end
        % Check if converged 
        if ~isempty(opt.xChangeThresh) && (xdiffAll(nit+2) < opt.xChangeThresh)
            % Here the change in x is tiny
            finished = 1;
            status = 10;
        end
        if ~isempty(opt.resThresh) && (res_sq_new < opt.resThresh)
            % Here the residual is tiny
            finished = 1;
            status = 1;
        end
        if ~isempty(opt.resChangeThresh) && ( abs((resSqAll(nit+2)-resSqAll(nit+1))/resSqAll(nit+1)) < opt.resChangeThresh )
            % Here the change in residual is tiny
            finished = 1;
            status = 2;
        end
        if (~isempty(stopThresh) && (res_sq_new < stopThresh)) 
            % Here the residual dropped below the residual expected from
            % noise.
            NitMax = nit + opt.nseExtraIt;
            stopThresh = 0;
            status = 3;
        end
        if opt.stopOnResInc>0 && (resSqAll(nit+2)/resSqAll(nit+1) > opt.resIncThresh)
            resIncsTotal = resIncsTotal + 1;
            if resIncsTotal >= opt.stopOnResInc
                finished = 1;
                status = 4;
            end
        end
        if opt.stopOnXdifInc>0 && (xdiffAll(nit+2) > xdiffAll(nit+1)) 
            if ~lastxdifInc
                xdifIncsTotal = xdifIncsTotal + 1;
                lastxdifInc = 1;
                if xdifIncsTotal >= opt.stopOnXdifInc
                    finished = 1;
                    status = 5;
                end
            end
        else
            lastxdifInc = 0;
        end
        % Finish up current iteration
        p = res + betaFact * (res_sq_new / res_sq_old) * p;
        res_sq_old = res_sq_new;
        xold = x;
        nit = nit+1;
        if nit >= NitMax
            finished = 1;
            status = 0;
        end
    end
    resSqAll = resSqAll(1:(nit+1));
    xnormAll = xnormAll(1:(nit+1));
    mseAll = mseAll(1:(nit+1));
    xdiffAll = xdiffAll(1:(nit+1));
end

function y = Asub(x,transp,Ain)
    if numel(Ain) == 1
        % Likely a object of some custom class
        if strcmp(transp,'notransp')
            y = Ain*x;
        else
            y = Ain'*x;
        end
    else
        % Likely a mask of some kind
        if strcmp(transp,'notransp')
            y = Ain.*x;
        else
            y = conj(Ain).*x;
        end
    end
end