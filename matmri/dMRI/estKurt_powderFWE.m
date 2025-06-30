function [Dtissue,K_lte,K_ste,sigTissue,sigCSF] = estKurt_powderfwe(sigVals_lte,sigVals_ste,bVals_lte,bVals_ste,D_CSF,sigTissue0,sigCSF0,opt)
%
%   Estimate signal fractions of CSF and tissue assuming b-tensor kurtosis model with 2 compartments.
%       i.e., free water elimination
%   Assumes a powder average has already been taken at each b-value.        
%
%   sigVals_lte:    [N_bval_lte * Nx * Ny * Nz] matrix
%   sigVals_ste:    [N_bval_ste * Nx * Ny * Nz] matrix
%   bVals_lte:      [N_bval_lte * 1] vector. [s/mm^2]
%   bVals_ste:      [N_bval_ste * 1] vector. [s/mm^2]
%   D_CSF:          presumed diffusion coefficient of CSF. [mm^2/s]
%   sigTissue0:     starting guess for tissue b0 signal; [Nx * Ny * Nz] matrix
%   sigCSF0;        starting guess for CSF b0 signal; [Nx * Ny * Nz] matrix
%   options:        (optional) structure that contains options for algorithm.
%                   See comments where defaults are set in code for info.
%
%   (c) 2022, Corey Baron

    if nargin<8 || ~isfield(opt,'NitEstKurt') || isempty(opt.NitEstKurt) || isempty(opt)
        opt.NitEstKurt = 100;
    end
    if ~isfield(opt,'sigTisThresh') || isempty(opt.sigTisThresh)
        % Only try to compute diffusion coefficients etc if frac of signal in tissue is larger than this
        opt.sigTisThresh = 0.1;
    end
    if ~isfield(opt,'KsteThresh') || isempty(opt.KsteThresh)
        % Threshold for suspiciously Kste 
        % Small neg vals are allowed because microscopic kurtosis, which
        % contributes to Kste, can theoretically be negative
        opt.KsteThresh = -0.5;
    end
    if ~isfield(opt,'verbose') || isempty(opt.verbose)
        % Print out information
        opt.verbose = 1;
    end
    
    if nargin<5 || nargin<6 || nargin<7
        % Here a tissue only model will be used
        D_CSF = [];
        sigTissue0 = [];
        sigCSF0 = [];
    end
    
    if isempty(D_CSF) || isempty(sigTissue0) || isempty(sigCSF0) 
        % Simply assume only tissue in voxels
        %fprintf('Determining lte and ste kurtosis assuming tissue only (i.e., no CSF)\n')
        % Estimate Dtissue and kurtosis metrics by simultaneously solving 
        % the following two equations using least squares
        % ln(sigVals_lte) = ln(sigTissue) - b*Dtissue + b^2*V_lte  
        % ln(sigVals_ste) = ln(sigTissue) - b*Dtissue + b^2*V_ste
        y_lte = log(sigVals_lte(:,:));
        y_ste = log(sigVals_ste(:,:));
        y = [y_lte; y_ste]; clear y_lte y_ste
        sz_ste = [numel(bVals_ste),1];
        sz_lte = [numel(bVals_lte),1];      
        A_lte = [ones(sz_lte,'like',bVals_lte), -bVals_lte(:), bVals_lte(:).^2, zeros(sz_lte,'like',bVals_lte)];
        A_ste = [ones(sz_ste,'like',bVals_lte), -bVals_ste(:), zeros(sz_ste,'like',bVals_lte), bVals_ste(:).^2];
        A = [A_lte; A_ste]; clear A_lte A_ste
        x = A\y;
        sigTissue(:) = exp(x(1,:));
        Dtissue(:) = x(2,:);
        V_lte(:) = x(3,:);
        V_ste(:) = x(4,:);
        sigCSF = zeros(size(sigTissue), 'like', sigTissue);
    else
        % Free water elimination model
        
        % Copy variables to gpu 
        if isgpuarray(sigVals_lte) && ~isgpuarray(D_CSF)
            D_CSF = gpuArray(D_CSF);
        end
        
        % Iteratively determine sig fractions and Dtissue
        sz = [size(sigVals_lte), 1, 1, 1, 1, 1];
        sigTissue = sigTissue0;
        sigCSF = sigCSF0;
        Dtissue = zeros(sz(2:end),'like',sigVals_lte);
        V_lte = zeros(sz(2:end),'like',sigVals_lte);
        V_ste = zeros(sz(2:end),'like',sigVals_ste);
        fprintf('Free water elimination estimation: 0%%');
        lastprint = 0;
        for n=1:opt.NitEstKurt
            % Estimate Dtissue and kurtosis metrics by simultaneously solving 
            % the following two equations using least squares
            % ln( (sigVals_lte - sigCSF*exp(-b*D_CSF)) / sigTissue ) = -b*Dtissue + b^2*V_lte  
            % ln( (sigVals_ste - sigCSF*exp(-b*D_CSF)) / sigTissue ) = -b*Dtissue + b^2*V_ste
            % values of sigTissue close to 0 will make this unstable.
                % We replace those values with 1 before solving the system,
                % then replace with actual values again later
            sigTissue_b = sigTissue;
            sigTissue(abs(sigTissue./sigCSF) < opt.sigTisThresh) = sigCSF(abs(sigTissue./sigCSF) < opt.sigTisThresh);
            y_lte = (sigVals_lte(:,:) - sigCSF(:).'.*exp(-bVals_lte(:)*D_CSF))./sigTissue(:).';
            y_ste = (sigVals_ste(:,:) - sigCSF(:).'.*exp(-bVals_ste(:)*D_CSF))./sigTissue(:).';
            y_lte = log(abs(y_lte)); % A voxel with only CSF can end up rounding to a small neg num
            y_ste = log(abs(y_ste));
            y = [y_lte; y_ste]; clear y_lte y_ste
            A_lte = [-bVals_lte(:), bVals_lte(:).^2, zeros(numel(bVals_lte),1,'like',bVals_lte)];
            A_ste = [-bVals_ste(:), zeros(numel(bVals_ste),1,'like',bVals_lte), bVals_ste(:).^2];
            A = [A_lte; A_ste]; clear A_lte A_ste
            x = A\y;
            Dtissue(:) = x(1,:);
            V_lte(:) = x(2,:);
            V_ste(:) = x(3,:);

            % Regions that have very small CSF signal will cause instabilities
            % in Dtissue
            sigTissue(abs(sigTissue_b./sigCSF) < opt.sigTisThresh) = sigTissue_b(abs(sigTissue_b./sigCSF) < opt.sigTisThresh);
            Dtissue(abs(sigTissue./sigCSF) < opt.sigTisThresh) = 0;
            V_lte(abs(sigTissue./sigCSF) < opt.sigTisThresh) = 0;
            V_ste(abs(sigTissue./sigCSF) < opt.sigTisThresh) = 0;
            
            % Check for negative D or K, which can happen from a poor
            % estimation of signal ratio (CSF sig estimation too high)
            inds = or(abs(Dtissue) < 0, or(V_lte<0, V_ste < opt.KsteThresh*Dtissue.^2 / 6) );
            if any(inds)
                % Adust the signal ratio there. More stable numerically to
                % come at this from the other side during iterations 
                % (i.e., start with too low CSF sig instead of too high)
                sigTissue(inds) = sigTissue(inds) + sigCSF(inds);
                sigCSF(inds) = 0;
                y_lte = (sigVals_lte(:,inds) - sigCSF(inds).'.*exp(-bVals_lte(:)*D_CSF))./sigTissue(inds).';
                y_ste = (sigVals_ste(:,inds) - sigCSF(inds).'.*exp(-bVals_ste(:)*D_CSF))./sigTissue(inds).';
                y = [log(abs(y_lte)); log(abs(y_ste))]; clear y_lte y_ste
                x = A\y;
                Dtissue(inds) = x(1,:);
                V_lte(inds) = x(2,:);
                V_ste(inds) = x(3,:);
            end
              
            % For any remaining negative D or K values to be 0
            Dtissue(Dtissue < 0) = 0;
            V_lte(V_lte < 0) = 0;
            V_ste(V_ste < opt.KsteThresh*Dtissue.^2/6) = opt.KsteThresh*Dtissue(V_ste < opt.KsteThresh*Dtissue.^2/6).^2/6;

            if opt.NitEstKurt>1
                % Estimate b0 signal values. Only bother if multiple iterations
                % Simultaneously solve the following equations with least
                % squares:
                % sigVals_lte = sigTissue * exp(-b*Dtissue+b^2*V_lte) + sigCSF * exp(-b*D_CSF) 
                % sigVals_ste = sigTissue * exp(-b*Dtissue+b^2*V_ste) + sigCSF * exp(-b*D_CSF) 
                A_lte = [exp(-bVals_lte(:).*reshape(Dtissue(:),1,1,[]) + bVals_lte(:).^2.*reshape(V_lte(:),1,1,[])),...
                    exp(-bVals_lte(:).*repmat(D_CSF(:),1,1,numel(Dtissue)))];
                A_ste = [exp(-bVals_ste(:).*reshape(Dtissue(:),1,1,[]) + bVals_ste(:).^2.*reshape(V_ste(:),1,1,[])),...
                    exp(-bVals_ste(:).*repmat(D_CSF(:),1,1,numel(Dtissue)))];
                A = [A_lte; A_ste]; clear A_lte A_ste

                if isgpuarray(A)
                    % pagefun is perfect for gpu
                    y = [sigVals_lte; sigVals_ste];
                    y = reshape(y, size(y,1),1,[]);
                    x = pagefun(@mldivide,A,y);
                    x = squeeze(x);
                else
                    % Make into a sparse block diagonal matrix to efficiently solve
                    % multiple matrix equations
                    j = repmat(1:size(A,2),[size(A,1) 1 size(A,3)]) + size(A,2)*reshape(0:size(A,3)-1,1,1,[]);
                    i = 1:size(A,1);
                    i = repmat(i(:), [1, size(A,2) size(A,3)]) + size(A,1)*reshape(0:size(A,3)-1,1,1,[]);
                    A = sparse(i(:),j(:),A(:));

                    % Solve the system
                    y = [sigVals_lte; sigVals_ste];
                    x = A\y(:);
                    x = reshape(x, [2, size(sigTissue)]);
                end

                sigTissue(:) = x(1,:);
                sigCSF(:) = x(2,:);

                % Account for negative signals. Move the signal so that the
                % net b0 signal stays constant
                sigTissue(sigCSF < 0) = sigTissue(sigCSF < 0) + sigCSF(sigCSF < 0);
                sigCSF(sigCSF < 0) = 0;
                sigCSF(sigTissue < 0) = sigTissue(sigTissue < 0) + sigCSF(sigTissue < 0);
                sigTissue(sigTissue < 0) = 0;
            end

            if opt.verbose && (n/opt.NitEstKurt*100 - lastprint > 4.99)
                if lastprint<10
                    fprintf('\b\b')
                else
                    fprintf('\b\b\b')
                end
                fprintf('%d%%', round(n/opt.NitEstKurt*100));
                lastprint = round(n/opt.NitEstKurt*100);
            end
            if n == opt.NitEstKurt && numel(Dtissue) == 1
                fprintf('sigTissue = %g; sigCSF = %g; Dtissue = %g; nit = %d\n', sigTissue, sigCSF, Dtissue, n);
            end
        end
        fprintf('\n');
    end
    
    % Convert to kurtosis maps
    K_lte = 6*V_lte./Dtissue.^2;
    K_ste = 6*V_ste./Dtissue.^2;
    K_lte(Dtissue == 0) = 0;
    K_ste(Dtissue == 0) = 0;
    
    % Clip non-physical kurtosis
    K_lte(K_lte > 4) = 4;
    K_ste(K_ste > 4) = 4;
    K_lte(K_lte < opt.KsteThresh) = opt.KsteThresh;
    K_ste(K_ste < opt.KsteThresh) = opt.KsteThresh;

end

