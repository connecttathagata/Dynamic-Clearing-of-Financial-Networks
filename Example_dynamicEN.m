%% Initial Setup
% close all
clear
format long

%Problem Parameters
cutoff=1e-10;
dt0 = 1e-6;
x0 = 100;
rng(465768789)

%System Parameters
n = 4;
x = [1;3;2;5];
L = [0,7,1,1;
    3,0,3,3;
    1,1,0,1;
    1,2,1,0];
x = [x;x0];
L = [L , 3*ones(n,1) ; zeros(1,n+1)];
n = n+1;
pbar = sum(L,2);
tmp = L./repmat(pbar,[1 n]);
tmp(isnan(tmp)) = 1/n;
Pi = tmp - diag(diag(tmp));

[p,D] = EN(x,Pi,pbar);
equity = x + Pi.'*p - pbar;
equity(end) = equity(end) - x0;

T = 1; %terminal time


%% BROWNIAN BRIDGE calculator for dc
mu = @(t,c)(x + sum(L,1).' - sum(L,2) - c)/(1-t);
% sigma = @(t,c)0*eye(n,n);
% sigma = @(t,c)1*eye(n,n);
sigma = @(t,c)5*eye(n,n);


%% System Calculator and Analysis + Graphs
[time,V,A,c,faroff]=continuousAlg(dt0,T,x,mu,sigma,@(t)L);
V(end,:) = V(end,:) - x0;

figure; hold on
plot(time,V(5,:),time,V(1,:),time,V(2,:),time,V(3,:),time,V(4,:),'Linewidth',1)
legend('Society: V_0(t)-100','Bank 1: V_1(t)','Bank 2: V_2(t)','Bank 3: V_3(t)','Bank 4: V_4(t)','Location','southwest');
plot([0,1],[0,0],'k:')
gax = gca; gax.ColorOrderIndex = 1;

gax.ColorOrderIndex = 1;

scatter([T;T],[equity(5);equity(5)],50,'filled')
scatter([T;T],[equity(1);equity(1)],50,'filled')
scatter([T;T],[equity(2);equity(2)],50,'filled')
scatter([T;T],[equity(3);equity(3)],50,'filled')
scatter([T;T],[equity(4);equity(4)],50,'filled')
axis([0 1 -10 10])
xlabel('Time')
ylabel('Wealth')

