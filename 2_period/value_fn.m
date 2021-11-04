% Compute the value function of two period RBC
function V = value_fn(x,w,theta,pi)
k_0 = x(1);         % initial capital stock
A_0 = x(2);        % initial TFP
% transition function determines the next state A
if w<1/3
    A_1 = .5;           % bad TFP
elseif w<2/3
    A_1 = 1;            % normal TFP
else
    A_1 = 1.5;         % good TFP
end

alpha = theta(1);       % capital share of production
beta = theta(2);         % discount rate
V = log(A_0*k_0^alpha-pi) + beta*log(A_1) + alpha*beta*log(pi);

%{
% scale the value function to be bounded in [0,1] a.s.
% by modifying the utility function
scle1 = 50;
scle2 = 100;
%V = log(A_0*k_0^alpha-pi)+scle1 + beta*log(A_1) + alpha*beta*log(pi)+beta*scle1;
%V = V/scle2;
%}
end