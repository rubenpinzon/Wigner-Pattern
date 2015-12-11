function stats = classGPFA(P, folds, debug, models, varargin)
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


scale = false;
if nargin>4
   scaleK = varargin{1}; 
   scale  = true;
   
   disp('Scaling the GP Kernel')
end

for ifold = 1 : folds
    posterior   = zeros(length(models), length(P));
    proto_traj  = cell(length(models), length(P));


    for p = 1 : length(P) 
        for m = 1 : length(models)
            %select the model parameters from the fold#1 
            model = models{m}.params{ifold};
            
            %rescale time scale of the GP if needed.
            if scale
               model.gamma = model.gamma * scaleK;  
            end
            
            [traj, ll] = exactInferenceWithLL(P(p), model,'getLL',1);        
            posterior(m,p) = ll;
            proto_traj{m,p} = traj;      

        end       

        [~, max_mod]  = max(posterior(:,p));
        loglike_p(p)  = posterior(max_mod, p);        

        if debug    
            %get joint neural space
            model = result_D.params{1};
            [traj, ll] = exactInferenceWithLL(P_seg(p), model,'getLL',1);
            figure(ifold)
            Xorth = orthogonalize([traj.xsm], model.C);
            plot3(Xorth(1,:),Xorth(2,:),Xorth(3,:),'color',colors(max_mod,:)), hold on
            plot3(Xorth(1,1),Xorth(2,1),Xorth(3,1),'color',...
                colors(max_mod,:),'markerfacecolor',colors(max_mod,:),'marker','s')
            drawnow
            fprintf('Proto event %d\n',p);

        end
    end

    %max log like P(event|model)
    [val, best_mod]  = max(posterior);
    type            = [P.trialType]; %{P(proto|model) , realtag}


    TP            = sum(best_mod == 1 & type == 1)/(sum(type == 1));
    FN            = sum(best_mod == 2 & type == 2)/(sum(type == 2));
    FP            = sum(best_mod == 1 & type == 2)/(sum(type == 2));
    TN            = sum(best_mod == 2 & type == 1)/(sum(type == 1));
    
    stats(ifold).conf_matrix    = [TP, FP; TN, FN];
    stats(ifold).class_output   = best_mod;
    stats(ifold).real_label     = type;
    stats(ifold).posterior      = posterior;
    fprintf('Fold %d done\n',ifold)
end