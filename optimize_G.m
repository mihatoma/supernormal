clear; clc;

% Define the prior
mu = [1/3; 1/3; 1/3];

% Define the neighborhoods for each X_i
neighborhoods = [1 2];

% Define the omegas and gammas
omega = [1 2 3];
gamma = 1;
%gamma = [1 2 3];

% Define the actions and corresponding utilities
A = [0 1];
u = [-1; 1; 0.5]; %utilities for a=1
u=transpose(u);

% Define the constants
I = length(neighborhoods);
kappa = [0.5, 0.5];

% Define the objective function
f = @(pi_gamma) - (gross_payoff(mu, pi_gamma, u) - cost(mu, pi_gamma, neighborhoods, kappa));

% Define the lower and upper bounds for pi
lb = zeros(3, 3);
ub = ones(3, 3);

% Define the equality constraints
Aeq = ones(1,9);
%Aeq = ones(1,3);
beq = 1;

% Define the inequality constraints
Aineq = [];
bineq = [];

% Define the optimization options
options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'iter');

% Perform the optimization
pi_gamma0 = ones(3, 3) / 9;  % Initial guess for posteriors
problemx = createOptimProblem('fmincon','objective',f,'x0',pi_gamma0,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
[pi_gamma, fval] = fmincon(problemx);

% Print the optimal value of the objective function and pi
fprintf('Optimal objective function value: %f\n', -fval);
disp('Optimal pi:');
disp(pi_gamma);


function G = gross_payoff(mu, pi_gamma, u)
    omega = [1 2 3];
    gamma = [1 2 3];
    G = 0;
    for i = 1:length(gamma)
        for j = 1:length(omega)
            G = G + (mu(j).* pi_gamma(i, j))*(max(0,u(j).*mu(j).*pi_gamma(i, j)));
        end
    end
end

function K = cost(mu, pi_gamma, neighborhoods, kappa)
    omega = [1 2 3];
    gamma = [1 2 3];
    K = 0;
    for i = 1:length(neighborhoods)
        X_i = i;
        for j = 1:length(gamma)
            K_intermediary = mu(j) * kappa(X_i); % First summation element
            for k=1:length(omega)
                H1 = -(pi_gamma(j,k))*log(pi_gamma(j,k));
                H2 = -(mu(j,1).*log(mu(j,1)));
                K = K + pi_gamma(j, k)' * (H1 - H2); % Second summation element
            end
            K_end = K_intermediary + K; % Together
        end
    end
end
