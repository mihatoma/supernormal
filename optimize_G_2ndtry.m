clear; clc;

% Define the prior
mu = [1/3; 1/3; 1/3];

% Define the neighborhoods for each X_i
neighborhoods = [1 2];

% Define the omegas and gammas
omega = [1 2 3];
gamma = [1 2 3];
%gamma = 1;

% Construct variables to place results of payoffs and costs
res_gross_payoff = zeros(3,1);
res_cost = zeros(3,1);

% Define the actions and corresponding utilities
A = [0 1];
u = [-1; 1; 0.5]; %utilities for a=1
u=transpose(u);

% Define the constants
I = length(neighborhoods);
kappa = [0.5, 0.5];
pi_gamma0 = ones(1, 3) / 9; % initial guess for pi_gamma

% Define the objective function
f = @(pi_gamma) - (gross_payoff(mu, [pi_gamma(1) pi_gamma(2) pi_gamma(3)], u) - cost(mu, [pi_gamma(1) pi_gamma(2) pi_gamma(3)], neighborhoods, kappa));
sol=f(pi_gamma0)

% Define the lower and upper bounds for pi
lb = zeros(1, 3);
ub = ones(1, 3);

% Define the equality constraint
Aeq = ones(1, 3); % ones(1,3)
beq = 1;

% Define the inequality constraints
Aineq = [];
bineq = [];

% Define the optimization options
options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'iter'); %'quasi-newton'

% Perform the optimization
pi_gamma0 = ones(1, 3) / 10;  % Initial guess for posteriors
problemx = createOptimProblem('fmincon','objective',f,'x0',pi_gamma0,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
[pi_gamma, fval] = fmincon(problemx);

% Print the optimal value of the objective function and pi
fprintf('Optimal objective function value: %f\n', -fval);
disp('Optimal pi:');
disp(pi_gamma);


function G = gross_payoff(mu, pi_gamma, u)
    omega = [1 2 3];
    gamma = [1 2 3];
    G_end = zeros(3,1);
    for i = 1:length(gamma)
        for j = 1:length(omega)
            bool=(u' .* mu .* pi_gamma');
            if (bool(j,1) > 0)
                G = (mu .* pi_gamma') .* (u' .* mu .* pi_gamma');
            else
                G = (mu' .* pi_gamma)*0;
            end
            
        end
        G_end = G_end+G;
    end
    res_gross_payoff(i,1)=G_end(i,1);
end


function K = cost(mu, pi_gamma, neighborhoods, kappa)
    %omega = [1 2 3];
    gamma = [1 2 3];
    K = 0;
    K_end = 0;
    for i = 1:length(neighborhoods) % no. of neighbourhoods for omega_1 and omega_2
        X_i = i; 
        K_intermediary = mu .* kappa(X_i); % First summation element
        for k=1:length(gamma) 
            H1 = -(pi_gamma) .* log(pi_gamma);
            H2 = -(mu .* log(mu));
            K = K + (mu .* pi_gamma' .* (H1' - H2)); % Second summation element
        end
        K_end = K_end + K_intermediary .* K; % Together
    end
    res_cost(i,1) = K_end(i,1);
end


