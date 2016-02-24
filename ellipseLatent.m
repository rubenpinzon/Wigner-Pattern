function ellipseLatent(varargin)


figure()
set(gcf,'color','w')

colors  = hsv(256);

diff_col = false;
if nargin == 1
    diff_col = true;
end

for v = 1 : nargin
    varX = varargin{v};
    p_prime = [mean([varX(1:2:end,2),varX(2:2:end,2)]) 2];         %to aling the ellipse with the neural trajectory (two neighbours)
    for n_bin = 1 : len(varX)

        data = [varX(1:2:end,n_bin),varX(2:2:end,n_bin)];          % data to compute the ellipse. This are accross trials
        [e_x, e_y] = ellipse_eig(data);                                        % ellipse
        e_z        = n_bin*ones(len(e_x), 1);                                  % Locate ellipse in Z-space     

        E          = [e_x, e_y, e_z];    
        origin     = mean(E);  

        angle      = acos(dot(origin, p_prime)/(norm(origin)*norm(p_prime)));

        T          = rotx(angle);                                              % Rotation matrix     
        rotE       = E;
        
        colorK     = colors(100*v,:);
        if diff_col
            colorK = colors(n_bin,:);
        end
        
        plot3(rotE(:,1),rotE(:,2),rotE(:,3),'color',colorK), hold on
        %plot3(data(:,1),data(:,2),n_bin*ones(1,len(data)),'*','color',colors(n_bin,:))

        if n_bin < len(varX)
            p_prime = [mean([varX(1:2:end,n_bin+1),varX(2:2:end,n_bin+1)]) n_bin+1]; 
        else
            p_prime = [mean([varX(1:2:end,n_bin),varX(2:2:end,n_bin)]) n_bin]; 
        end
    end
    plot3(mean(varX(1:2:end,:)),mean(varX(2:2:end,:)),1:len(varX),'color','k','linewidth',2)
    
end
xlabel('Latent X1','fontsize',14)
ylabel('Latent X2','fontsize',14)
zlabel('Time','fontsize',14)
grid on