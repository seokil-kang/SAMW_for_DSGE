function EEE = EEE_SAMW(x,theta,pi)
k_0 = x(1);
A_0 = x(2);
alpha = theta(1);
beta = theta(2);
y_0 = A_0*k_0^alpha;
c_0 = y_0-pi;
A_1 = [.5;1;1.5];
prob_A = ones(size(A_1))/length(A_1);
c_1 = A_1*pi^alpha;
EEE = (c_0- (beta*sum(prob_A.*A_1./c_1*alpha*pi^(alpha-1)))^(-1))/c_0;
end