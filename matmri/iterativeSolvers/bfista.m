function [x, resSqAll, RxAll, mseAll, maxEig] = bfista(Ain,bin,Rin,lam,x0,NitMax,opt)
    % Use balanced FISTA to solve argmin(||Ax-b||^2_2 + lam*||Rx||_1)  
    %
    % x = fista(Ain,bin,Rin,x0,NitMax,options)
    %
    %   A must have a transpose that can be evaluated using A'*b, or as a function with A(b,'transp')
    %   A operates on x via A*x or A(x,'notransp')
    %   x0 and output of Ain are expected to be Nx1 vectors
    %
    %   Can enforce a bounded support in x using options.xMask 
    %   See comments in code for information on other options.
    %
    %   For FISTA fundamentals, see Beck A, Teboulle M. A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems. SIAM J. Imaging Sci. 2009;2:183–202.
    %   For balanced FISTA, see Ting ST, Ahmad R, Jin N, et al. Fast implementation for
    %       compressive recovery of highly accelerated cardiac cine MRI using
    %       the balanced sparse model. Magn Reson Med. 2017; 77: 1505-1515
    %
    % (c) Corey Baron 2021
    %
    
    % Set options
    if nargin<6 || isempty(NitMax)
        % Maximum number of iterations allowed
        NitMax = 200;
    end
    if nargin<7 || ~isfield(opt,'maxEig')
        % Maximum eigenvalue of A'A
        opt.maxEig = []; 
    end
    if nargin<7 || ~isfield(opt,'resThresh')
        % Threhold for residuals (to automatically stop iterations)
        opt.resThresh = 1e-4; 
    end
    if nargin<7 || ~isfield(opt,'gtruth')
        % Useful for simulations. Allows computation of mse per iteration
        opt.gtruth = []; 
    end
    if nargin<7 || ~isfield(opt,'verbose')
        % whether to provide information
        opt.verbose = 1; 
    end
    if nargin<7 || ~isfield(opt,'xMask')
        % whether to enforce bounded x-space support
        opt.xMask = []; 
    end
    if nargin<7 || ~isfield(opt,'kMask')
        % whether to enforce bounded k-space support
        opt.kMask = []; 
    end
    
    % Account for different ways of supplying A
    if isa(Ain,'function_handle')
        A = Ain;
    else
        A = @(x,transp) Asub(x,transp,Ain);
    end

    % Account for different ways of supplying Rin
    if isa(Rin,'function_handle')
        R = Rin;
    else
        R = @(x,transp) Asub(x,transp,Rin);
    end
    
    % Set default starting guess
    if nargin<3 || isempty(x0)
        x0 = A(bin,'transp');
    end
    
    % Define the step size.  This must be less than the inverse of the
    % smallest Lipschitz constant, which in our l1 regularization problem
    % is equal to twice the maximum eigenvalue of A'*A (see example 2.2 in Beck et al). We make it a little
    % smaller to be safe.
    if ~isempty(opt.maxEig)
        maxEig = opt.maxEig;
    else
        warning('No opt.maxEig provided (recommended to precompute for a sample slice and input with opt.maxEig, since value should be similar for all slices of a single acquisition)')
        fprintf('bfista.m: finding maximum eigenvalue of A''A using power method...')
        tic1 = tic;
        maxEig = powermethod(A,x0);
        fprintf('took %d sec\n', round(toc(tic1)));
    end
    stepSz = gather(0.9/(2*abs(maxEig)));

    % Initialize
    y = x0;
    resSqAll = zeros(NitMax+1,2);
    RxAll = zeros(NitMax+1,2);
    mseAll = zeros(NitMax+1,1);
    if ~isempty(opt.gtruth)
        tmp = opt.gtruth(:)-x0(:);
        mseAll(1) = tmp'*tmp;
    end
    
    % Get first points for tracking convergence of residuals
    if nargout > 1
        warning('Residual output is requested. This DRASTICALLY increases computation time, and should only be done for demos/development/debugging')
        residual = A(x0,'notransp')-bin;
        resSqAll(1,1) = residual(:)'*residual(:);
    end
    if nargout > 2
        l1norm = abs(lam*(R(x0,'notransp')));
        RxAll(1,1) = gather(sum(l1norm(:)));
    end
    if (nargout > 3) && ~isempty(opt.gtruth)
        tmp = opt.gtruth(:)-x0(:);
        mseAll(1) = tmp'*tmp;
    end

    % Iterate
    finished = 0;
    nit = 1;
    t_prev = 1;
    x_prev = x0;
    strmsg = [];
    while ~finished
        % Find gradient for ||Ax-b||^2_2
        residual_y = A(y,'notransp')-bin;
        grad = 2*A(residual_y,'transp');
        resSqAll(nit,2) = residual_y(:)'*residual_y(:);
        
        % Perform the step along the gradient (this is just simple gradient decent)
        g = y - stepSz*grad;
        
        % Perform soft thresholding of the transform
        % See Beck et al., eq 2.6 for why lam is multiplied by stepSz. 
        % Basically, it is because of the factor of L in front of the l2
        % norm (our stepSz is equiv to L in eq 2.6). The solution to 2.6 is
        % softthresholding when g(x) is the l1 norm
        % Note that our implementation of the wavelet here uses balanced
        % FISTA (Ting et al).
        %   See eq 10
        %   Note that Rin'*Rin = I must be true for this case. This is
        %   true for both the decimated and undecimated wavelet tranform.
        x = R(g,'notransp');
        tmp = abs(lam*x);
        RxAll(nit,2) = gather(sum(tmp(:)));
        x = softthresh(x,lam*stepSz); 
        x = R(x,'transp');

        % Constrain support x-space using projection
        if ~isempty(opt.xMask)
            x = x.*opt.xMask;
        end

        % Constrain support k-space using projection. Note that the
        % supplied mask should be fftshifted. 
        %   This can help to control noise amplification, since samples
        %   outside the k-space support (e.g., corners of k-space when a
        %   spherical sampling pattern is used) are poorly conditioned.
        %   Still not sure if this is beneficial (CB), since wavelet
        %   regularization controls noise amplification anyway. When I've
        %   tried this in combination with wavelet, image looks comparable
        %   but ringing is worsened due to the sharp edges in k-space introduced.
        if ~isempty(opt.kMask)
            if (ndims(opt.kMask)==2) && (size(opt.kMask,2)>1)
                x = fftnc(x,2,0,0);
                x = x.*opt.kMask;
                x = ifftnc(x,2,0,0);
            else
                error('TODO: code for other dims, or require an fftFunc supplied to bfista')
            end
        end
        
        % Update tracking
        if nargout > 1
            residual = A(x,'notransp')-bin;
            resSqAll(nit,1) = residual(:)'*residual(:);
        end
        if nargout > 2
            % Note that Rin'*Rin = I does NOT ensure that Rin*Rin' = I
            % (e.g., undecimated wavelet xform), which is why we have to
            % re-evalue the transform
            l1norm = abs(lam*(R(x,'notransp')));
            RxAll(nit,1) = gather(sum(l1norm(:)));
        end
        if (nargout > 3) && ~isempty(opt.gtruth)
            tmp = opt.gtruth(:)-x(:);
            mseAll(nit) = tmp'*tmp;
        end
        
        % Perform FISTA step
        t = 0.5*(1 + sqrt(1+4*t_prev^2));
        y = x + (t-1)/t_prev*(x-x_prev);
        
        % Update vars
        t_prev = t;
        x_prev = x;
        
        % Check for completion
        if nit >= NitMax
            strmsg = 'NitMax';
            finished = 1;
        end
        if nit>1
            testVal = abs((resSqAll(nit,2)-resSqAll(nit-1,2))/resSqAll(nit-1,2)) +...
                abs((RxAll(nit,2)-RxAll(nit-1,2))/RxAll(nit-1,2));
        end
        if (nit>1) && testVal<opt.resThresh
            strmsg = 'opt.resThresh';
            finished = 1;
        end
        nit = nit+1;
    end
    if opt.verbose
        fprintf('bfista: stopped on %s after %d iterations\n', strmsg, nit);
    end

    % Clean up
    resSqAll = resSqAll(1:nit,:);
    RxAll = RxAll(1:nit,:);
    mseAll = mseAll(1:nit,:);
    
end

function y = Asub(x,transp,Ain)
    if strcmp(transp,'notransp')
        y = Ain*x;
    else
        y = Ain'*x;
    end
end