% Compute the value function of two period RBC
function V = value_fn(x,w,theta,pi)
k_0 = x(1);         % initial capital stock
A_0 = x(2);        % initial TFP
% transition function determines the next state A
if w<=1/9
    A_1 = .5; 
    A_2 = .5;
elseif w<=2/9
    A_1 =.5; 
    A_2 = 1;
elseif w<=3/9
    A_1 =.5; 
    A_2 = 1.5;
elseif w<=4/9
    A_1 = 1; 
    A_2 = .5;
elseif w<=5/9
    A_1 = 1; 
    A_2 = 1;
elseif w<=6/9
    A_1 = 1; 
    A_2 = 1.5;
elseif w<=7/9
    A_1 = 1.5; 
    A_2 = .5;
elseif w<=8/9
    A_1 = 1.5; 
    A_2 = 1;
else
    A_1 = 1.5; 
    A_2 = 1.5;
end

alpha = theta(1);       % capital share of production
beta = theta(2);         % discount rate
pi_0 = pi(1);                % policy at period 0
pi_1 = pi(2);                % policy at period 1

% value function conditional on all values above
if pi_1 >= A_1*pi_0^alpha
    V = -inf;
else
V = log(A_0*k_0^alpha-pi_0) + beta*log(A_1*pi_0^alpha-pi_1) + beta^2*log(A_2*pi_1^alpha);
end
end