function Xtats = classGPFA(P, models, varargin)
%CLASSGPFA  Given a data struct P containing spike count vectors and GPFA models, this file computes a
%           binary classification as the argmax P(data|each_model).
%
%           INPUTS:
%           P           : DataHigh type struct of dimension (1 x num events)
%           folds       : number of folds used for crossvalidation during training
%           debug       : shows debugging and verbose output
%           models      : cell containing trained GPFA models with dimension (1 x num models)
%
%           OUTPUT:
%           stats       : a struct with dimensions (1 x folds) including the fields
%                       conf_matrix, class_output, real_label, and posterior, which are
%                       the classification confusion matrix where positive samples correspond
%                       to right alternations whereas negative samples are left alternations;
%                       output of the classifier {1:right, 2:left}, real label, and the
%                       log posterior P(data|model).
%
%
%Version 1.0 Ruben Pinzon@2015
folds        = length(models{1}.params);
useAllTrials = false;
scale = false;
if nargin == 3
   scaleK = varargin{1}; 
   if ~isempty(scaleK)
    scale  = true;   
    fprintf('Scaling the GP Kernel with %2.2f\n',scaleK)
   end
elseif nargin == 4
   useAllTrials = true;
   disp('Warning: Using all the trials for testing')
end

n_laps      = length(P);
v_laps      = [P.trialId];
model_like  = zeros(length(models), n_laps);
    
for m = 1 : length(models)
    likelikehood   = -Inf*ones(folds, n_laps);

    for ifold = 1 : folds
        
        if ~useAllTrials
            usedlaps    = models{m}.trainTrials{ifold};
            unseenP     = ones(1,n_laps);
            for u = 1 : length(usedlaps)
                u_idx = find(v_laps == usedlaps(u));
                unseenP(u_idx) = 0;
            end
            unseenP = find(unseenP ==1);
        else
            unseenP = 1:n_laps;
        end
        
        for p = 1 : length(unseenP) 
        
            %select the model parameters from the fold#1 
            model = models{m}.params{ifold};
            lap   = unseenP(p);
            %rescale time scale of the GP if needed.
            if scale
               model.gamma = model.gamma * scaleK;
            end
            
            [traj, ll] = exactInferenceWithLL(P(lap), model,'getLL',1);        
            likelikehood(ifold,lap) = ll;                  

        end       
        %remove trials used during training
    end
    
    model_like(m,:) = max(likelikehood);
end

[~, max_mod]    = max(model_like);

type            = [P.type]; %{P(proto|model) , realtag}


TP            = sum(max_mod == 1 & type == 1)/(sum(type == 1));
FN            = sum(max_mod == 2 & type == 2)/(sum(type == 2));
FP            = sum(max_mod == 1 & type == 2)/(sum(type == 2));
TN            = sum(max_mod == 2 & type == 1)/(sum(type == 1));

Xtats.conf_matrix    = [TP, FP; TN, FN];
Xtats.class_output   = max_mod;
Xtats.real_label     = type;
Xtats.likelihood     = model_like;

%fprintf('Fold %d done\n',ifold)
