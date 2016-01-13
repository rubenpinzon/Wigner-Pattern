function trainLDA(W,P)
%TRAINLDA perfroms linear classification of data P [features X samples] using
%crossvalidation strategy and the labels provided in W.type. len(W) = samples.
%it uses the bayesian linear classifier by James Barrett@2013, BLC lib
%
%
%Ruben Pinzon @2015

type  = [W.type];                                                          %using only correct trials, thus, type 1 and 2
W     = W(type ==1 | type == 2);                                           %filter correct trials
P     = P(:,type ==1 | type == 2)';                                        %remove error trials from data 
sigma = [W.type]';                                                         %labels
sigma(sigma==2) = -1;                                                      %conver type 2 to -1 as reqired by BLC.lib

[w, w0, werror, w0error, logL] = blc_train(P,sigma)


data = [P sigma];
[N, d] = size(data);
Nminus = sum(data(:,d)==-1);
data = sortrows(data,d);

Xminus = data(1:Nminus,1:d-1);
Xplus = data(Nminus+1:N,1:d-1);

plot(Xplus(:,1),Xplus(:,2),'bo')
hold on
plot(Xminus(:,1),Xminus(:,2),'ro')
legend('+1 class','-1 class')

x = -3:1:3;
y=-(w(1)/w(2))*x -(w0/w(2))*ones(1,length(x));
plot(x,y,'-k')


xlabel('x_1')
ylabel('x_2')
grid on