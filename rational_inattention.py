import numpy as np
from scipy.optimize import minimize

def gross_payoff(mu, pi_gamma, u):
    omega = np.array([1, 2, 3])
    gamma = np.array([1, 2, 3])
    G_end = np.zeros((3, 3))
    for i in range(len(gamma)):
        for j in range(len(omega)):
            bool_val = (u * mu * pi_gamma)[j, 0]
            if (bool_val > 0):
                G = (mu * pi_gamma) * (u * mu * pi_gamma)
            else:
                G = (mu.T * pi_gamma) * 0

        G_end += G
    return G_end


def cost(mu, pi_gamma, neighborhoods, kappa):
    gamma = [1, 2, 3]
    K = 0
    K_end = 0
    for i in range(len(neighborhoods)):  # no. of neighbourhoods for omega_1 and omega_2
        X_i = i + 1
        K_intermediary = mu * kappa[X_i - 1]  # First summation element
        for k in range(len(gamma)):
            H1 = -(pi_gamma) * np.log(pi_gamma)
            H2 = -(mu * np.log(mu))
            K = K + (mu * pi_gamma.T * (H1.T - H2))  # Second summation element
        K_end = K_end + K_intermediary * K  # Together
    return K_end


# Define the prior
mu = np.array([1 / 3, 1 / 3, 1 / 3])

# Define the neighborhoods for each X_i
neighborhoods = np.array([1, 2])

# Define the omegas and gammas
omega = np.array([1, 2, 3])
gamma = np.array([1, 2, 3])
# gamma = 1

# Define the actions and corresponding utilities
A = np.array([0, 1])
u = np.array([-1, 1, 0.5])  # utilities for a=1
u = u.reshape((3, 1))

# Define the constants
I = len(neighborhoods)
kappa = np.array([0.5, 0.5])
pi_gamma0 = np.ones(3) / 9  # initial guess for pi_gamma


# Define the objective function
def f(pi_gamma):
    return - (gross_payoff(mu, np.array([pi_gamma[0], pi_gamma[1], pi_gamma[2]]), u) - cost(mu, np.array(
        [pi_gamma[0], pi_gamma[1], pi_gamma[2]]), neighborhoods, kappa))


sol = f(pi_gamma0)

# Define the bounds for pi
bounds = [(0, 1)] * 3


# Define the equality constraint
def constraint(pi_gamma):
    return np.sum(pi_gamma) - 1


eq_cons = {'type': 'eq', 'fun': constraint}

# Define the optimization options
options = {'disp': True}

# Perform the optimization
pi_gamma0 = np.ones(3) / 10  # Initial guess for posteriors
res = minimize(f, pi_gamma0, method='SLSQP', bounds=bounds, constraints=[eq_cons], options=options)

# Print the optimal value of the objective function and pi
print('Optimal objective function value: %f' % -res.fun)
print('Optimal pi:')
print(res.x)

