function [estParams, seq, LL, iterTime] = update_gps(currentParams, seq, emMaxIters, varargin)
%
% From myem.m in the livrary DataHigh
% Copyright Benjamin Cowley, Matthew Kaufman, Zachary Butler, Byron Yu, 2012-2013



%   emMaxIters   = 150;
  tol          = 1e-8;
  minVarFrac   = 0.01;
  verbose      = false;
  freqLL       = 10;
  extra_opts   = assignopts(who, varargin);
  
  N            = length(seq(:));
  T            = [seq.T];
  [yDim, xDim] = size(currentParams.C);
  LL           = [];
  LLi          = 0;
  iterTime     = [];
  varFloor     = minVarFrac * diag(cov([seq.y]'));
  
  % Loop once for each iteration of EM algorithm
  

  for i = 1:emMaxIters   
    
  
    rand('state', i);
    randn('state', i);
    if nargout > 3
        tic;
    end
    
    if (rem(i, freqLL) == 0) || (i<=2)
      getLL = true;
    else
      getLL = false;
    end
            
    % ==== E STEP =====
    if ~isnan(LLi)
      LLold = LLi;
    end
    [seq, LLi] = exactInferenceWithLL(seq, currentParams, 'getLL', getLL);    
    LL         = [LL LLi];
    
    %update the GP lenght scale.
    res = learnGPparams(seq, currentParams, 'verbose', verbose,... 
                          extra_opts{:});
    switch currentParams.covType
    case 'rbf'
      currentParams.gamma = res.gamma;
    case 'tri'
      currentParams.a     = res.a;
    case 'logexp'
      currentParams.a     = res.a;
    end
    if currentParams.notes.learnGPNoise  
        currentParams.eps = res.eps;  
    end
    
    tEnd     = toc;
    if nargout > 3        
        iterTime = [iterTime tEnd];
    end

    % Display the most recent likelihood that was evaluated
    if getLL
       fprintf('Iter %d  :  lik %g (%.1f sec)\n',i, LLi, tEnd);
    end
    % Verify that likelihood is growing monotonically
    if i<=2
      LLbase = LLi;
    elseif (LLi < LLold)
      fprintf('\nError: Data likelihood has decreased from %g to %g\n',... 
        LLold, LLi);
      keyboard;
    elseif ((LLi-LLbase) < (1+tol)*(LLold-LLbase))
        fprintf('\nError: Data likelihood has decreased from %g to %g\n',... 
        LLold, LLi);
      break;
    end
  end
  
  if length(LL) < emMaxIters  % the user quit dim reduction, so return nothing
    estParams = [];
    return;
  end
  
  if any(diag(currentParams.R) == varFloor)
    fprintf('Warning: Private variance floor used for one or more observed dimensions in GPFA.\n');
  end
 
  estParams = currentParams;
