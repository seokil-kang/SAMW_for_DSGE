% This code execute the algorithm Simulated Annealing Multiplicative Weights for the simple two-period RBC model with full depreciation
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
beta = .95;             % % discount rate 
% exogenous value
k_0 = 20;        % initial capital stock
A_0 = 1;        % initial TFP
y_0= A_0*k_0^alpha;     % initial production

% vectors for convenience
x = [k_0, A_0];              % state vector
theta = [alpha;beta];   % parameter vector

%---------------------------------------------------
%                simulation setup
%---------------------------------------------------
% Simulation parameter
N = 2*10^4;               % iteration number

% Construct the policy space
% If we let Lambda to have y_0 as its element, it will make the value function diverges to negative infinity
% This is innocuous in computing optimal policy and value function,
% however, theorem for finite-time upper bound fails 
% because it requries an assumption 0<V<1 for any policy
k = 100;         % # of threshold of space
%Lambda = [0:y_0/(k+1):y_0]';        % HEURISTIC finite policy space
%Lambda = Lambda(2:end-1);   % exclude the exhausting policy which is not HEURISTIC
Lambda = [0:y_0/(k-1):y_0]';        % HEURISTIC finite policy space

% Construct the distribution on the policy space
phi_1 = ones(k,1)/k;    % initial distriubtion is uniform
PHI = [phi_1];                  % distribution sequence

% updating the distribution
 %gamma = 2;        % annealing parameter which converges faster then variable case

% baskets for computation
V_fn = zeros(k,N);      % value function basket
V_fn(:,1) = value_fn(x,rand,theta,Lambda);

%---------------------------------------------------
%                         simulation
%---------------------------------------------------

tic;
for i = 2:N
    % draw a random seed for transition function
    w = rand;
    
    % compute the value function for all policy in Lambda
    V_fn(:,i) = value_fn(x,w,theta,Lambda);
    
    % update the distribution sequence
    gamma = 1+sqrt(1/(i-1));              % decreasing annealing parameter
    %gamma = 2;                                      % annealing parameter which converges faster then variable case
    Z = PHI(:,end)'*gamma.^V_fn(:,i);    %  normalizing factor
    phi_update = PHI(:,end).*gamma.^V_fn(:,i)/Z;
    PHI = [PHI phi_update];
    if mod(i*100,N)==0
        clc
        fprintf('Simulation iterated %.0f%s\n',i*100/N,'%');
    end
end
toc;
phi_s = PHI(:,end);
%%
result_SAMW