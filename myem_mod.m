function [estParams, seq, LL, iterTime] = myem_mod(currentParams, seq, varargin)
%A modified version of MYEM file in the DataHIGH library.
%
%
% [estParams, seq, LL] = em(currentParams, seq, ...)
%
% NOTE: The only difference between this and the original em.m it was copied from
% is that it does not print as much output.
%
% Fits GPFA model parameters using expectation-maximization (EM) algorithm.
%
%   yDim: number of neurons
%   xDim: state dimensionality
%
% INPUTS:
%
% currentParams - GPFA model parameters at which EM algorithm is initialized
%                   covType (string) -- type of GP covariance ('rbf')
%                   gamma (1 x xDim) -- GP timescales in milliseconds are 
%                                       'stepSize ./ sqrt(gamma)'
%                   eps (1 x xDim)   -- GP noise variances
%                   d (yDim x 1)     -- observation mean
%                   C (yDim x xDim)  -- mapping between low- and high-d spaces
%                   R (yDim x yDim)  -- observation noise covariance
% seq           - training data structure, whose nth entry (corresponding to
%                 the nth experimental trial) has fields
%                   trialId      -- unique trial identifier
%                   T (1 x 1)    -- number of timesteps
%                   y (yDim x T) -- neural data
%
% OUTPUTS:
%
% estParams     - learned GPFA model parameters returned by EM algorithm
%                   (same format as currentParams)
% seq           - training data structure with new fields
%                   xsm (xDim x T)        -- posterior mean at each timepoint
%                   Vsm (xDim x xDim x T) -- posterior covariance at each timepoint
%                   VsmGP (T x T x xDim)  -- posterior covariance of each GP
% LL            - data log likelihood after each EM iteration
% iterTime      - computation time for each EM iteration
%               
% OPTIONAL ARGUMENTS:
%
% emMaxIters    - number of EM iterations to run (default: 500)
% tol           - stopping criterion for EM (default: 1e-8)
% minVarFrac    - fraction of overall data variance for each observed dimension
%                 to set as the private variance floor.  This is used to combat
%                 Heywood cases, where ML parameter learning returns one or more
%                 zero private variances. (default: 0.01)
%                 (See Martin & McDonald, Psychometrika, Dec 1975.)
% freqLL        - data likelihood is computed every freqLL EM iterations. 
%                 freqLL = 1 means that data likelihood is computed every 
%                 iteration. (default: 5)
% verbose       - logical that specifies whether to display status messages
%                 (default: false)
% wbar          - keeps track of a waitbar for DimReduce
%
% @ 2009 Byron Yu         byronyu@stanford.edu
%        John Cunningham  jcunnin@stanford.edu
%  Copyright Benjamin Cowley, Matthew Kaufman, Zachary Butler, Byron Yu, 2012-2013



  emMaxIters   = 500;
  tol          = 1e-8;
  minVarFrac   = 0.01;
  verbose      = false;
  freqLL       = 10;
  wbar         = -1;  %waitbar handle for cross-validation
  dimreduce_wb           = -1;  % waitbar handle for dim reduction
  bar_time     = 0;   % current progress of the waitbar
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
      
    % ==== M STEP ====    
    sum_Pauto   = zeros(xDim, xDim);
    for n = 1:N
      sum_Pauto = sum_Pauto + ...
          sum(seq(n).Vsm, 3) + seq(n).xsm * seq(n).xsm'; 
    end
    Y           = [seq.y];
    Xsm         = [seq.xsm];
    sum_yxtrans = Y * Xsm';
    sum_xall    = sum(Xsm, 2);
    sum_yall    = sum(Y, 2);

    term = [sum_Pauto sum_xall; sum_xall' sum(T)]; % (xDim+1) x (xDim+1)
    Cd   = ([sum_yxtrans sum_yall]) / term;   % yDim x (xDim+1)
    
    currentParams.C = Cd(:, 1:xDim);
    currentParams.d = Cd(:, end);
    
    % yCent must be based on the new d
    % yCent = bsxfun(@minus, [seq.y], currentParams.d); 
    % R = (yCent * yCent' - (yCent * [seq.xsm]') * currentParams.C') / sum(T);        
    if currentParams.notes.RforceDiagonal
      sum_yytrans = sum(Y .* Y, 2);
      yd          = sum_yall .* currentParams.d;
      term        = sum((sum_yxtrans - currentParams.d * sum_xall') .*... 
                        currentParams.C, 2);
      r           = currentParams.d.^2 + (sum_yytrans - 2*yd - term) / sum(T);
      
      % Set minimum private variance
      r               = max(varFloor, r);
      currentParams.R = diag(r);
    else      
      sum_yytrans = Y * Y';
      yd          = sum_yall * currentParams.d';
      term        = (sum_yxtrans - currentParams.d * sum_xall') * currentParams.C';
      R           = currentParams.d * currentParams.d' +...
                    (sum_yytrans - yd - yd' - term) / sum(T);
                    
      currentParams.R = (R + R') / 2; % ensure symmetry
    end
    
    if currentParams.notes.learnKernelParams
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
    end
    
    if nargout > 3
        tEnd     = toc;
        iterTime = [iterTime tEnd];
    end

    % Display the most recent likelihood that was evaluated
    if getLL
        fprintf('Iter %d  :  lik %g\n',i, LLi);
    end
    % Verify that likelihood is growing monotonically
    if i<=2
      LLbase = LLi;
    elseif (LLi < LLold)
      fprintf('\nError: Data likelihood has decreased from %g to %g\n',... 
        LLold, LLi);
      keyboard;
    elseif ((LLi-LLbase) < (1+tol)*(LLold-LLbase))
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
