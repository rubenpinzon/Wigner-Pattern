function stats = classGPFA(P, debug, models, varargin)
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
%Rev 2.0 Dec 27: only process testing samples.
%Version 2.0 Ruben Pinzon@2015

folds = length(models{1}.params);

scale = false;
if nargin>4
   scaleK = varargin{1}; 
   scale  = true;
   
   fprintf('Scaling the GP Kernel with %2.2f\n',scaleK)
end

colors = hsv(length(P));

likelikehood   = zeros(length(models), length(P));
proto_traj     = cell(length(models), length(P));
    
%trials in the testing set of each model    
for m = 1 : length(models)
    for ifold = 1 : folds
        %select the model parameters from the fold#1 
        model       = models{m}.params{ifold};
        test_p      = models{m}.testTrials{ifold};

        %rescale time scale of the GP if needed.
        if scale
           model.gamma = model.gamma * scaleK;
        end
        for t = 1 : length(test_p)
            [traj, ll] = exactInferenceWithLL(P(test_p(t)), model,'getLL',1);        
            likelikehood(m,test_p(t)) = ll;
            proto_traj{m,test_p(t)} = traj;
        end

    end
end

for p = 1 : length(P)
    [~, max_mod]  = max(likelikehood(:,p));
    loglike_p(p)  = likelikehood(max_mod, p);        

    if debug    
        %get joint neural space
        model = models{m}.params{1};
        traj  = exactInferenceWithLL(P(p), model,'getLL',0);
        figure(ifold)
        Xorth = orthogonalize([traj.xsm], model.C);
        plot3(Xorth(1,:),Xorth(2,:),Xorth(3,:),'color',colors(max_mod,:)), hold on
        plot3(Xorth(1,1),Xorth(2,1),Xorth(3,1),'color',...
            colors(max_mod,:),'markerfacecolor',colors(max_mod,:),'marker','s')
        plot3(Xorth(1,end),Xorth(2,end),Xorth(3,end),'color',...
            colors(max_mod,:),'markerfacecolor',colors(max_mod,:),'marker','o')
        drawnow
        if p == 1
            fprintf('Observation %d..', p); 
        else
            fprintf('%d..', p);
        end

    end
end
disp('done')

%max log like P(event|model)
[val, best_mod]  = max(likelikehood);
type            = [P.type]; %{P(proto|model) , realtag}


TP            = sum(best_mod == 1 & type == 1)/(sum(type == 1));
FN            = sum(best_mod == 2 & type == 2)/(sum(type == 2));
FP            = sum(best_mod == 1 & type == 2)/(sum(type == 2));
TN            = sum(best_mod == 2 & type == 1)/(sum(type == 1));

stats.conf_matrix   = [TP, FP; TN, FN];
stats.class_output  = best_mod;
stats.real_label    = type;
stats.likelihood    = likelikehood;
stats.traj          = proto_traj;

