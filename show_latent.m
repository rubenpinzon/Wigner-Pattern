function Xorth = show_latent(model, data)
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


for m = 1 : n_models    
    
    traj     = exactInferenceWithLL(data{m}, model{m},'getLL',0);
    x        = orthogonalize([traj.xsm], model{m}.C);     
    
    Xorth{m} = x;
    
    T              = [0 cumsum([traj.T])];
        
    figure(50+m)
    set(gcf, 'position', [1,1,1424,973], 'color', 'w')
       
    start_traj  = []; end_traj = [];
    for ilap = 1 : length(traj)
       lap_t = T(ilap)+1:T(ilap+1);
       c     = traj(ilap).color; %color of the lap
       plot_xorth(x(1,lap_t),x(2,lap_t),x(3,lap_t),[1 2 4 5 7 8],{'X_1','X_2','X_3'},c,num2str(traj(ilap).trialId))           
       plot_xorth(x(1,lap_t),x(2,lap_t),[],3,{'X_1','X_2'},c)
       plot_xorth(x(2,lap_t),x(3,lap_t),[],6,{'X_2','X_3'},c)
       plot_xorth(x(1,lap_t),x(3,lap_t),[],9,{'X_1','X_3'},c)          

       start_traj(ilap, :) = x(1:3,lap_t(1));
       end_traj(ilap, :)   = x(1:3,lap_t(end));
       
    end
    %covariance ellipses
    ellipse_eig(end_traj(:,1:2), 3, [1, 0, 0])
    ellipse_eig(end_traj(:,2:3), 6,[1, 0, 0])
    ellipse_eig(end_traj(:,[1,3]), 9,[1, 0, 0])
    ellipse_eig(start_traj(:,1:2), 3, [0, 0, 1])
    ellipse_eig(start_traj(:,2:3), 6,[0, 0, 1])
    ellipse_eig(start_traj(:,[1,3]), 9,[0, 0, 1])        
    clear x traj

end