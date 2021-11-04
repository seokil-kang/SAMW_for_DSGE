% This code execute the algorithm Simulated Annealing Multiplicative Weights for the simple three-period RBC model with full depreciation
%  Code is based on the textbook Simulation-Based Algorithms for Markov Decision Processes  by Chang et al.
% Last update: 02/01/2019
% Written by Seokil Kang (Indiana U)

%---------------------------------------------------
%                    House keeping
%---------------------------------------------------
clear;close all;clc;

%---------------------------------------------------
%                    model setup
%---------------------------------------------------
% Parameter value
alpha = 1/3;        % capital share of production
beta = .9;             % % discount rate 
% exogenous value
k_0 = 5;        % initial capital stock
A_0 = 1;        % initial TFP
y_0= A_0*k_0^alpha;     % initial production

% vectors for convenience
x = [k_0, A_0];              % state vector
theta = [alpha;beta];   % parameter vector

%---------------------------------------------------
%                simulation setup
%---------------------------------------------------
% Simulation parameter
N = 5*10^2;               % iteration number

% Construct the policy space
k = 100;         % # of threshold of each subspace
Lambda_0 = [0:y_0/(k-1):y_0]';           % HEURISTIC finite policy space for period 0
y_1 = 1.5*(y_0)^alpha;
kk = 100;
Lambda_1 = [0:y_1/(kk-1):y_1]';        % HEURISTIC finite policy space for period 1
Lambda = [];
for j = 1:k
    for jj = 1:kk
        lambda = [Lambda_0(j), Lambda_1(jj)];
        Lambda = [Lambda;lambda];
    end
end

% Construct the distribution on the policy space
phi_1 = ones(k*kk,1)/k/kk;    % initial distriubtion is uniform
PHI = [phi_1];                  % distribution sequence

% updating the distribution
 gamma = 2;        % annealing parameter which converges slower then variable case

% baskets for computation
V_fn = zeros(k,kk);      % value function basket

% convergence of distribution
convergence = 100;
tol = 10^-4;
i = 0;

%---------------------------------------------------
%                         simulation
%---------------------------------------------------

tic;
for i = 2:N
%while convergence > tol
   % i = i+1;
    % draw a random seed for transition function
    w = rand;
    % compute the value function for all policy in Lambda
    for j = 1:k
        for jj = 1:kk
            V_fn(j,jj) = value_fn(x,w,theta,[Lambda_0(j),Lambda_1(jj)]);
        end
    end
    V_fn = reshape(V_fn',k*kk,1);
    
    % update the distribution sequence
    gamma = 1+sqrt(1/(i));                  % decreasing annealing parameter
    Z = PHI(:,end)'*gamma.^V_fn;    %  normalizing factor
    phi_update = PHI(:,end).*gamma.^V_fn/Z;
    PHI = [PHI phi_update];
    V_fn = reshape(V_fn',k,kk);
    if mod(i,50) == 0
        clc
        fprintf('Current simulation iteration = %.0f\n',i)
    end
    %convergence = norm(PHI(:,end-1)-PHI(:,end));
end
clc
fprintf('Total simulation iteration = %.0f\n',i)

%%
%---------------------------------------------------
%                         simulation
%---------------------------------------------------

phi_s = PHI(:,end);

% analytic optimal policy
pi_star_0 = (alpha*beta+(alpha*beta)^2)/(1+alpha*beta+(alpha*beta)^2)*A_0*k_0^alpha;
% A very weird trouble occurs here: the optimal policy at period 1 must  have E[A_1] = .8 in order to match the simulation result!! I cannot understand this point so far
pi_star_1 = alpha*beta/(1+alpha*beta)*1*pi_star_0^alpha;
pi_star = [pi_star_0;pi_star_1];
toc;

%{
figure
scatter3(Lambda(:,1),Lambda(:,2),phi_s,'.','linewidth',25);
xlabel('\pi_0')
ylabel('\pi_1')
%}
figure
hold on
p1 = plot(Lambda(:,1)',phi_s,'.','markersize',12,'linewidth',2.5);
p2 = plot(Lambda(:,2)',phi_s,'.','markersize',12,'linewidth',2.5);
line([pi_star(1) pi_star(1)],[0, max(phi_s)],'linestyle','--','color','c','linewidth',1.5);
line([pi_star(2) pi_star(2)],[0, max(phi_s)],'linestyle','--','color','m','linewidth',1.5);
hold off
legend([p1 p2],'\phi for \pi_0', '\phi for \pi_1')