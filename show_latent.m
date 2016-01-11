function Xorth = show_latent(model, data, colors, labels)
%SHOW_LATENT shows the latent trajectories in the 3D space defined by the 3
%       largest singe values decomposition vectors of A in the model Y = Cx+d+e
%       inputs:
%
%
%       model : {1 x n_models} : parameters of the GPFA model {C, d, R, gamma, eps, covType}
%       data  : {1 x n_models} : including the field y(binned spike bins)
%
%Ruben Pinzon@2015


n_models = length(model);
fprintf('%d models provided\n',n_models);
figure
for m = 1 : n_models    
    Params   = model{m}.params{1}; % fold #1
    traj     = exactInferenceWithLL(data, Params,'getLL',0);
    x        = orthogonalize([traj.xsm], Params.C);     
    Xorth{m} = x;
    
    T              = [0 cumsum([traj.T])];
        
    set(gcf, 'position', [1,1,1424,973], 'color', 'w')
       
    start_traj  = []; end_traj = [];
    for ilap = 1 : length(traj)
       lap_t = T(ilap)+1:T(ilap+1);   
       c        = colors(labels(ilap),:); %this color is model, has to be the trial type

       plot_xorth(x(1,lap_t),x(2,lap_t),x(3,lap_t),[1 2 4 5 7 8],{'X_1','X_2','X_3'},c,num2str(traj(ilap).trialId))           
       plot_xorth(x(1,lap_t),x(2,lap_t),[],3,{'X_1','X_2'},c)
       plot_xorth(x(2,lap_t),x(3,lap_t),[],6,{'X_2','X_3'},c)
       plot_xorth(x(1,lap_t),x(3,lap_t),[],9,{'X_1','X_3'},c)   

      
    end
    %covariance ellipses    
    clear x traj
end