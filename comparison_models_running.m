%Scrip to compare the GPFA models trained with  and without the middle arm in
%the maze running section
%
%
%Ruben Pinzon@2015
clear all, clc

% with middle arm
name = {'_branch2_results40ms.mat', '_branch2_noMidArm.mat'};
arm = {'w', 'wth'};
conditions      = {'_left', '_right', ''};

% D = wth.D;
% %measuring the 
% for i = 1 : length(D)
%    m_count(i, :) = sum(D(i).data,2);    
% end

% fprintf('Mean+/-std spike count accross laps and cells %f+/-%f\n',mean(mean(m_count)),std(mean(m_count)))

figure(1)
set(gcf, 'position', [0 1 600 300], 'color', 'w')

c = lines(3);
m = {'s', 'o', 'v'};
for a = 1 : 2
    load(['/media/bigdata/i01_maze13.003/i01_maze13.003' name{a}]);  
    [D_left, D_right] = split_trails(D);

    for s = 1 : length(conditions)

        Data            = eval(sprintf('D%s',conditions{s}));
        Result          = eval(sprintf('result_D%s;',conditions{s}));
        mask            = false(1,length(Data)); % for cross validation
        cv_trials       = Result.cv_trials;
        fold_indx       = Result.foldidx;
        fprintf('Condition %s loaded\n',conditions{s})

        ifold           = 1;      
        % prepare masks:
        % test_mask isolates a single fold, train_mask takes the rest
        test_mask       = mask;
        test_mask(cv_trials(fold_indx(ifold):fold_indx(ifold+1)-1)) = true;
        test_data  = Data(test_mask);
        %Posterior of test data given the trained model
        [traj, ll_te] = exactInferenceWithLL(test_data, Result.params{ifold},'getLL',1);
        % orthogonalize the trajectories4
        [Xorth, Corth, TT, EE] = orthogonalize([traj.xsm], Result.params{ifold}.C);

        displayname = [arm{a} conditions{s}];
        exp_var = (EE.^2)/sum((EE.^2));
        plot(cumsum(exp_var), '-o','displayname',displayname, 'color',c(a,:),...
             'marker', m{s}),
        hold on, grid on


    end
end
ylabel('Variance Explained')
xlabel('Eigenvector')
print(gcf,'/media/bigdata/i01_maze13.003/i01_maze13.003_varExp.png','-dpng')

