function [sigTissue,sigCSF,Dtissue] = estD_powderFWE(sigVals,bVals,D_CSF,opt)
%
%   Estimate signal fractions of CSF and tissue assuming DTI model with 2 compartments. 
%       i.e., free water elimination
%   Assumes a powder average has already been taken at each b-value.    
%
%   sigVals:    [N_bval * Nx * Ny * Nz] matrix
%   bVals:      [N_bval * 1] vector. [s/mm^2]
%   D_CSF:      presumed diffusion coefficient of CSF. [mm^2/s]
%   Nit:        number of iterations to use
%   Dtissue0:   starting guess for diffusivity of tissue. [mm^2/s]
%
%   (c) 2022, Corey Baron

    % Set default options
    if nargin<4 || ~isfield(opt,'Dtissue0') || isempty(opt.Dtissue0)
        % Starting guess for tissue ADC
        opt.Dtissue0 = 7e-4;
    end
    if ~isfield(opt,'nIterLowb') || isempty(opt.nIterLowb)
        % Number of iterations
        opt.nIterLowb = 100;
    end
    if ~isfield(opt,'sigTisThresh') || isempty(opt.sigTisThresh)
        % Only try to compute diffusion coefficients etc if frac of signal in tissue is larger than this
        opt.sigTisThresh = 0.1;
    end
    if ~isfield(opt,'verbose') || isempty(opt.verbose)
        % Print out information
        opt.verbose = 1;
    end

    % Copy variables to gpu 
    Dtissue = opt.Dtissue0;
    if isgpuarray(sigVals) && ~isgpuarray(Dtissue)
        Dtissue = gpuArray(Dtissue);
    end
    if isgpuarray(sigVals) && ~isgpuarray(D_CSF)
        D_CSF = gpuArray(D_CSF);
    end
    
    % Iteratively determine sig fractions and Dtissue
    sz = [size(sigVals), 1, 1, 1, 1, 1];
    sigTissue = zeros(sz(2:end),'like',sigVals);
    sigCSF = zeros(sz(2:end),'like',sigVals);
    fprintf('Low b signal fraction estimation: 0%%');
    lastprint = 0;
    for n=1:opt.nIterLowb
        % Estimate b0 signal values 
        % Solving sigVals = sigTissue * exp(-b*Dtissue) + sigCSF * exp(-b*D_CSF) using least squares
        A = [exp(-bVals(:).*reshape(Dtissue(:),1,1,[])), exp(-bVals(:).*repmat(D_CSF(:),1,1,numel(Dtissue)))];
        if n > 1
            if isgpuarray(A)
                % pagefun is perfect for gpu
                y = reshape(sigVals, size(sigVals,1),1,[]);
                x = pagefun(@mldivide,A,y);
                x = squeeze(x);
            else
                % Make into a sparse block diagonal matrix to efficiently solve
                % multiple matrix equations. This seems faster than pagefun
                % for CPU
                j = repmat(1:size(A,2),[size(A,1) 1 size(A,3)]) + size(A,2)*reshape(0:size(A,3)-1,1,1,[]);
                i = 1:size(A,1);
                i = repmat(i(:), [1, size(A,2) size(A,3)]) + size(A,1)*reshape(0:size(A,3)-1,1,1,[]);
                A = sparse(i(:),j(:),A(:));
                x = A\sigVals(:);
                x = reshape(x, [2, size(sigTissue)]);
            end
        else
            % There is only one A matrix because there is only one Dtissue
            % value on the first iteration.
            x = A\sigVals(:,:);
        end
        sigTissue(:) = x(1,:);
        sigCSF(:) = x(2,:);
        
        % Account for negative signals. Move the signal so that the net
        % signal stays constant
        sigTissue(sigCSF < 0) = sigTissue(sigCSF < 0) + sigCSF(sigCSF < 0);
        sigCSF(sigCSF < 0) = 0;
        sigCSF(sigTissue < 0) = sigTissue(sigTissue < 0) + sigCSF(sigTissue < 0);
        sigTissue(sigTissue < 0) = 0;
        
        % Values of sigTissue close to 0 will make solution of Dtissue unstable.
            % We replace those values with the CSF val before solving the system,
            % then replace with actual values again later
        smallSigTissue_inds = abs(sigTissue./sigCSF) < opt.sigTisThresh;
        sigTissue_b = sigTissue(smallSigTissue_inds);
        sigTissue(smallSigTissue_inds) = sigCSF(smallSigTissue_inds);
        
        % Estimate Dtissue
        % solves "ln(sigVals - sigCSF*exp(-b*D_CSF)) = ln(sigTissue) - b*Dtissue" using least squares
        % Note that abs of y_a is taken, because in rare situations it can
        % have a small negative value
        y_a = sigVals(:,:) - sigCSF(:).'.*exp(-bVals(:)*D_CSF);
        y = log(abs(sigTissue(:).')) - log(abs(y_a));
        Dtissue = reshape(bVals(:)\y, size(sigCSF));
        
        % Don't allow negative diffusion coefficients. This will happen
        % where the sigCSF was set too high, so we set it to zero to
        % approach from the other direction (more robust)
        inds = Dtissue<0;
        sigTissue(inds) = sigTissue(inds) + sigCSF(inds);
        sigCSF(inds) = 0;
        y_a = sigVals(:,inds) - sigCSF(inds).'.*exp(-bVals(:)*D_CSF);
        y = log(abs(sigTissue(inds).')) - log(abs(y_a));
        Dtissue(inds) = reshape(bVals(:)\y, sum(inds), 1);
        
        % Remove any remaining negative diffusion coefficients
        Dtissue(Dtissue<0) = 0;

        % Replace the values we switched out earlier. Also set Dtissue
        % there to 0 because the values are meaningless
        sigTissue(smallSigTissue_inds) = sigTissue_b;
        Dtissue(smallSigTissue_inds) = 0;
        
        if opt.verbose && (n/opt.nIterLowb*100 - lastprint > 4.99)
            if lastprint<10
                fprintf('\b\b')
            else
                fprintf('\b\b\b')
            end
            fprintf('%d%%', round(n/opt.nIterLowb*100));
            lastprint = round(n/opt.nIterLowb*100);
        end
        if (n == opt.nIterLowb) && (numel(Dtissue) == 1)
            fprintf('sigTissue = %g; sigCSF = %g; Dtissue = %g; nit = %d\n', sigTissue, sigCSF, Dtissue, n);
        end
    end
    fprintf('\n');
    
    % Certain cases can lead to very small negative signal
    sigTissue = abs(sigTissue);
    sigCSF = abs(sigCSF);

end



