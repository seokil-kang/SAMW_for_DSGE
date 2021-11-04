%---------------------------------------------------
%                  simulation result
%---------------------------------------------------
clc;

% optimal policy and value function
pi_star = y_0*alpha*beta/(1+alpha*beta);
V_star = (value_fn(x,.1,theta,pi_star)+value_fn(x,.5,theta,pi_star)+value_fn(x,.9,theta,pi_star))/3;

% corner policies yield -inf: so let's convert them with a very low value like -10^5
V_fn(isinf(V_fn)) = -10^3;

% sequence of simulated policy
pi_bar = (Lambda'*PHI)';

% sequence of simulated value function
Psi_star=0;
for j = 1:N
    Psi_star = Psi_star + value_fn(x,rand,theta,pi_bar(end));
end
Psi_star = Psi_star/N;
V_bar = sum(V_fn.*PHI);
Psi_bar = sum(V_fn,2)/N;
indx = find(V_bar < .6 & V_bar > 0);
for j = 1:length(indx)
    V_bar_sort(j) = V_bar(indx(j));
end
N2 = length(indx);

% plot the distribution sequence
lw = 2.5;
mi = 30;


figure
hold on
p0 = line([pi_star pi_star],[0, max(PHI(:,end))],'linestyle',':','color','[.85 .33 .1]','linewidth',lw);
p1 = plot(Lambda,PHI(:,1),':*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.75,.75,.75]');
p3= plot(Lambda,PHI(:,25),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.7,.7,.7]');
p4 =plot(Lambda,PHI(:,250),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.65,.65,.65]');
p2 = plot(Lambda,PHI(:,1000),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.6,.6,.6]');
p5 =plot(Lambda,PHI(:,5000),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.55,.55,.55]');
p6 =plot(Lambda,PHI(:,10000),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.5,.5,.5]');
p3 = plot(Lambda,PHI(:,N),'-*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[0,0,0]');
hold off
axis tight
xlabel('policy space \Lambda')
ylabel('probability of policy \phi')
%text(Lambda(350),PHI(350,25),'i=10,000','fontsize',16)
legend([p0,p3],'optimal policy','\phi^N(stationary distribution)')
%legend([p0,p6],'optimal policy','\phi^{10000}(10,000th distribution)')
legend boxoff
%title('|\Lambda|=500, N = 10^5')



figure
mi = 1;
hold on
p0 = line([pi_star pi_star],[0, max(PHI(:,end))],'linestyle',':','color','[.85 .33 .1]','linewidth',lw);
p1 = plot(Lambda,PHI(:,1),':*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.75,.75,.75]');
plot(Lambda,PHI(:,2),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.7,.7,.7]');
plot(Lambda,PHI(:,5),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.65,.65,.65]');
p2 = plot(Lambda,PHI(:,10),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.6,.6,.6]');
plot(Lambda,PHI(:,25),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.55,.55,.55]');
plot(Lambda,PHI(:,50),'-.*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[.5,.5,.5]');
p3 = plot(Lambda,PHI(:,N),'-*','MarkerIndices',1:mi:length(Lambda),'linewidth',lw,'color','[0,0,0]');hold off
axis tight
xlabel('policy space \Lambda')
ylabel('probability of policy \phi')
%legend([p1,p3],'\phi^1(initial distribution)','\phi^N(stationary distribution)')
legend([p0,p3],'optimal policy','\phi^N(stationary distribution)')
legend boxoff
%title('|\Lambda|=10, N = 10^2')



figure
grid on
hold on
p1 = line([1 N],[pi_star pi_star],'linestyle',':','color','[.85 .33 .1]','linewidth',5);
p2 = plot(pi_bar,'-dk','MarkerIndices',[1:25:N],'markersize',7,'linewidth',2.25);
%p2 = plot(pi_bar,'-dk','MarkerIndices',[1:10:100, 101:1000:N],'markersize',5,'linewidth',2);
%p4 = plot(V_bar,'ok','MarkerIndices',1:50:N,'markersize',7);
%p4 = plot(V_bar,'ok','MarkerIndices',1:200:N,'markersize',7);
%p3 = line([1 N],[V_star V_star],'linestyle','--','color','b','linewidth',lw);
hold off
xlim([1 N])
%axis tight
xlabel('iteration period i')
ylabel('policy estimated at iteration i')
L = legend([p1,p2],'optimal policy','simulated optimal policy $$\sum_{\pi\in\Lambda}\phi^i(\pi)\pi$$');
%L = legend([p1,p2],'$$\pi^*$$','$$\sum_{\pi\in\Lambda}\phi^i(\pi)\pi$$');
set(L,'Interpreter','latex');
set(L,'fontsize',14);
%set(L,'location','southoutside');
legend boxoff


% Euler Equation Error
EEE = zeros(size(pi_bar));
for j = 1:length(pi_bar)
    EEE(j) = EEE_SAMW(x,theta,pi_bar(j));
end
EEE = log10(abs(EEE));

figure
grid on
hold on
%p2 = plot(EEE,'-dk','MarkerIndices',[1:5:50, 51:N/40:N],'markersize',7,'linewidth',2.25);
p2 = plot(EEE,'-k','linewidth',3);
hold off
xlim([1 N])
axis tight
xlabel('iteration period i')
ylabel('Euler Equation Error at iteration i of log_{10}')
%}

% report the main result
fprintf("optimal policy = %.4f\n",pi_star)
fprintf("simulated policy = %.4f\n",pi_bar(end))
fprintf("Euler equation error in log10 = %.4f\n",EEE(end))
fprintf("variance of sample policy = %.4f\n",(Lambda-pi_bar(end)).^2'*phi_s)
fprintf("optimal value function = %.4f\n",V_star)
fprintf("simulated value function = %.4f\n",Psi_star)
fprintf("finite-time upper bound = %.4f\n",(gamma-1)/log(gamma)*sum(V_bar)/N+log(length(Lambda))/N/log(gamma))
%}