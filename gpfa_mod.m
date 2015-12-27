function [gpfa_params, gpfa_traj, LL] = gpfa_mod(D,dims,varargin)
% GPFA_MOD A modified version of GPFA_ENGINE for the DataHigh program
%
%  Inputs:
%   D -- struct of trajectories, make sure it conforms to gpfa's format
%   dims --- the dimension you choose for the number of latent variables
%   emMaxIters (optional) --- number of iterations GPFA will go to
%   wb (optional) --- updates waitbar
%   binWidth (optional, default:20) --- if you want to change binWidth
%  Copyright Benjamin Cowley, Matthew Kaufman, Zachary Butler, Byron Yu,
%  John Cunningham, 2012-2013

% ---GNU General Public License Copyright---
% This file is a modified version of gpfa_engineHD.m inside DataHigh.
% 
% DataHigh is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 2.
% 
% 

startup_gpfa;  % run the script to set up MEX environment
emMaxIters      = 200;  % shown to work for most cases
bin_width       = 25; % 20ms at 1250 Hz in samples
extra_opts      = assignopts(who,varargin);

% Set maxIters option
if nargin < 3
    emMaxIters = 200;
end
gpfa_params     = [];
gpfa_traj       = [];
LL              = [];

% GPs paramters
startTau      = 100; % in msec
startEps      = 1e-3;
   
  
% ==================================
% Initialize state model parameters
% ==================================
startParams.covType = 'rbf';
% GP timescale
% Assume binWidth is the time step size.
startParams.gamma = (bin_width / startTau)^2 * ones(1, dims);
% GP noise variance
startParams.eps   = startEps * ones(1, dims);

% ========================================
% Initialize observation model parameters
% ========================================
yAll             = [D.y];
try
    [faParams, faLL] = fastfa(yAll, dims);
    startParams.d = mean(yAll, 2);
    startParams.C = faParams.L;
    startParams.R = diag(faParams.Ph);
catch
    d  = mean(yAll, 2);
    Y0 = yAll - repmat(d, 1, size(yAll,2));
    Q = cov(Y0');
    [C, Lambda] = eigs(Q, dims);
    fprintf('Initializing GPFA parameters with PCA\n');
    R = diag(diag(Q - C * Lambda * C'));
    
    startParams.d = d;
    startParams.C = C;
    startParams.R = R;
    disp('Succeded');
end



% Define parameter constraints
startParams.notes.learnKernelParams = true;
startParams.notes.learnGPNoise      = false;
startParams.notes.RforceDiagonal    = true;

currentParams = startParams;

  % =====================
  % Fit model parameters
  % =====================


  [gpfa_params, D, LL] =... 
    myem_mod(currentParams, D,'emMaxIters',emMaxIters);
      
  if(nargout > 1)
    % Extract orthonormalized neural trajectories for original, unsegmented trials
    % using learned parameters

    [gpfa_traj, LL_last] = exactInferenceWithLL(D, gpfa_params,'getLL',1);

    
    % orthogonalize the trajectories
    [Xorth, Corth] = orthogonalize([gpfa_traj.xsm], gpfa_params.C);
    gpfa_traj = segmentByTrial(gpfa_traj, Xorth, 'data');
    gpfa_traj = rmfield(gpfa_traj, {'Vsm', 'VsmGP', 'xsm'});
    gpfa_params.Corth = Corth;
  end
  

end