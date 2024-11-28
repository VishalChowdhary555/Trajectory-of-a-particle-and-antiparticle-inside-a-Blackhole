

% Constants
G = 6.67430e-11;    % Gravitational constant (m^3 kg^(-1) s^(-2))
c = 299792458;      % Speed of light (m/s)
h_bar = 1.0545718e-34;  % Reduced Planck constant (J s)
eV = 1.60218e-19;   % Electron volt (J)

% Black hole parameters
M = 2;    % Mass of the black hole (arbitrary units)
a = 0.9;    % Spin parameter (angular momentum per unit mass)

% Particle parameters
m0 = 1.0;   % Rest mass of the particle (arbitrary units)
charge = -1;   % Charge of the particle (in electron charge units)
potential = 1e9;  % Electric potential near the black hole (V)
% ... (your constants and functions)

% Constants for numerical integration
t_span = [0, 100];  % Time span for integration
initial_conditions = [1.5, pi/2, 0, 1.5, 1e-10];

% Numerical integration using ode45
[t, y] = ode45(@(t, y) kerr_particle_trajectory_with_hawking(t, y, M, a, m0, G, c, h_bar, charge, potential), t_span, initial_conditions);

% Extracting results
r = y(:, 1);
theta = y(:, 2);
phi = y(:, 3);
p_phi = y(:, 4);
energy = y(:, 5);

% Calculate additional quantities
angular_momentum = p_phi ./ sin(theta);

% Calculate the line element ds
ds = sqrt(-g_t_t(r, a, M, G, c) .* energy.^2 + g_phi_phi(a, r, G, c) .* angular_momentum.^2);

% Integrate ds to obtain proper time
proper_time = cumtrapz(t, ds) ./ (m0 * c^2);


% Display important aspects
disp('Results:');
disp(['Final Energy: ' num2str(energy(end)) ' J']);
disp(['Final Angular Momentum: ' num2str(angular_momentum(end)) ' J s']);
disp(['Proper Time of Fall: ' num2str(proper_time(end)) ' s']);
disp('Energy at different time points:');
disp([t, energy]);

% Plot key quantities during the simulation
figure;

subplot(2, 2, 1);
plot(t, r);
title('Radius vs Time');
xlabel('Time');
ylabel('Radius');

subplot(2, 2, 2);
plot(t, energy);
title('Energy vs Time');
xlabel('Time');
ylabel('Energy');

subplot(2, 2, 3);
plot(t, angular_momentum);
title('Angular Momentum vs Time');
xlabel('Time');
ylabel('Angular Momentum');

subplot(2, 2, 4);
% Ensure t and proper_time have the same length for plotting
plot(t, proper_time(1:length(t)));  % Use only the portion of proper_time corresponding to t
title('Proper Time vs Time');
xlabel('Time');
ylabel('Proper Time');

% Plot the trajectory of the particle and antiparticle pairs
figure;
plot3(r.*sin(theta).*cos(phi), r.*sin(theta).*sin(phi), r.*cos(theta));
title('Trajectory of Particle and Antiparticle Pairs');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

function dydt = kerr_particle_trajectory_with_hawking(t, y, M, a, m0, G, c, h_bar, charge, potential)
    % Kerr Metric parameters
    Delta = y(1)^2 - 2*M*y(1) + a^2;
    rho_sq = y(1)^2 + a^2*cos(y(2))^2;

    % Energy and angular momentum per unit mass at the current position
    epsilon = -g_t_t(y(1), a, M, G, c)*(y(4) + a*sin(y(2))^2);
    lambda = g_phi_phi(a, y(1), G, c)*(y(4) + a*sin(y(2))^2) - 2*a*g_t_phi(a, y(1), M, G, c)*y(4)/rho_sq;

    % Electric potential energy
    electric_potential = charge * potential;

    % Hawking radiation (conceptual)
    % Generate Hawking radiation based on a probability distribution
    hawking_probability = rand(size(t));  % Simplified random probability
    hawking_energy = h_bar * 2 * pi * c / (m0 * c^2);  % Simplified energy of Hawking particle

    % Adjust the energy of the particle based on Hawking radiation
    dydt = zeros(5, 1);
    dydt(5) = -electric_potential + hawking_probability * hawking_energy;

    % Gravitational potential term with consistent use of black hole mass
    dydt(5) = dydt(5) - G * M / y(1)^2;

    % Equations of motion
    dydt(1) = y(4)/(rho_sq*(epsilon^2 - 1));
    dydt(2) = sqrt(1 - epsilon^2)*cos(y(2))/rho_sq;
    dydt(3) = lambda/(rho_sq*sin(y(2))^2);
    dydt(4) = -((epsilon*(epsilon*y(4) - lambda))^2 - (rho_sq*Delta + a^2*epsilon^2)^2 - (a*epsilon)^2*Delta*sin(y(2))^2)/(rho_sq^2*Delta);
end

function g = g_t_t(r, a, M, G, c)
    g = -(1 - 2*M.*r./(r.^2 + a.^2))./(c^2);
end

function g = g_phi_phi(a, r, G, c)
    g = (r.^2 + a.^2 + 2*a.^2.*r./(r.^2 + a.^2).*sin(pi/4)^2).*sin(pi/4)^2;
end

function g = g_t_phi(a, r, M, G, c)
    g = -2*a*M*r/(r^2 + a^2)*sin(pi/4)^2/c;
end
