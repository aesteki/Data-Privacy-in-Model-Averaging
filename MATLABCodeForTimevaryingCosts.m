clc
close all
clear all

Aj = [0 1 1 1 1;
      1 0 1 0 1;
      1 1 0 1 0;
      1 0 1 0 1;
      1 1 0 1 0];
[N,~] = size(Aj);
D = diag(Aj*ones(N,1));
L = D - Aj;
I = eye(N);
name_alphac = [188.3;592.5;2567.2;1793.3;2567.2];
name_betac = [45.9;45.9;208.2;166.6;208.2];
d = 25000/N*ones(N,1);
initialorder = 10;
Nsim = 500000;
active = ones(N,1);
%Active Initial Conditions
Dt = .001;
v(:,1) = zeros(N,1);
z(:,1) = 4000*zeros(N,1);
w = 0.00001;

for k = 1:Nsim
    
    if k == 200000
        active(5) = 0;
    end
    name_alpha = name_alphac+name_alphac/100.*[sin(w*1*k);sin(w*2*k);sin(w*3*k);sin(w*4*k);sin(w*5*k)];
    name_beta = name_betac+name_alphac/100.*[cos(w*1*k);cos(w*2*k);cos(w*3*k);cos(w*4*k);cos(w*5*k)];
    d = 25000/N*ones(N,1)+25000/100.*[sin(w*k);sin(w*k);sin(w*k);sin(w*k);sin(w*k)];
    for i = 1:N
        if active(i) == 1
            alpha(i,1) = name_alpha(i);
            beta(i,1) = name_beta(i);
        else
            alpha(i,1) = 0;
            beta(i,1) = 0;
        end
    end
    p_star(:,k) = (sum(d)+sum(alpha))/sum(beta)*beta-alpha;
    %Active
    x(:,k) = z(:,k)+d+alpha;
    z(:,k+1) = z(:,k) - Dt*(beta.*x(:,k)-d-alpha+L*x(:,k)+L*v(:,k));
    v(:,k+1) = v(:,k) + Dt*L*x(:,k);
    p(:,k) = beta.*x(:,k)-alpha+(z(:,k+1)-z(:,k))/Dt;
    ep = p(:,k)-p_star(:,k);
    e(:,k) = log(sqrt(ep.*ep)/norm(p_star(:,k)));
end

figure
sim = 0:1:Nsim-1;
l1 = plot(sim,p,'b','LineWidth',2);
hold on
l2 = plot(sim,p_star,'r','LineWidth',2);
ylim([-1 15000])
set(gcf,'position',[100,100,500,400])
set(gca,'FontSize',20)
legend([l1(1,1),l2(1,1)],{'$p^i(k)$','$p^{i\star}(k)$'},'FontSize',20,'interpreter','latex')
xlabel('$k$','fontsize',27,'interpreter','latex')
ylabel('$p^i(k)$','fontsize',27,'interpreter','latex')
ax = gca;
ax.XAxis.Exponent = 3;
ax.YAxis.Exponent = 3;

figure
sim = 0:1:Nsim-1;
l1 = plot(0:1:200000-2,e(:,1:200000-1),'LineWidth',2);
hold on
l2 = plot(200000-1:1:Nsim-1,e(1:4,200000:Nsim),'LineWidth',2);
set(gcf,'position',[100,100,500,400])
set(gca,'FontSize',20)
%legend([l1(1,1),l1(2,1),l1(3,1),l1(4,1),l1(5,1)],{'$\textup{log}(e^1(k))$','$\textup{log}(e^2(k))$','$\textup{log}(e^3(k))$','$\textup{log}(e^4(k))$','$\textup{log}(e^5(k))$'},'FontSize',20,'interpreter','latex')
xlabel('$k$','fontsize',27,'interpreter','latex')
ylabel('$\textup{log}(e^i(k))$','fontsize',27,'interpreter','latex')
ax = gca;
ax.XAxis.Exponent = 3;
ax.ColorOrder = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880]