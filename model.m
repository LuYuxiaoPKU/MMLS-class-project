% Impementation of the ODE Model from Santurio and Barros

size = 120;

% k_?
% Parameters (from supplementary table 1)
p_C = 0.9; % CAR-T cell proliferation rate
g_T = 1e10; % T cell concentration for half-maximal CAR-T cell proliferation
tau_C = 7; % CAR-T cell lifespan
alpha = 1e-11; % Tumor cell inactivation rate
omega_T = 0.012; % Glioblastoma proiliferation rate
k = 8.5e11; % Carrying capacity
psi_T = 2.571e-15; 
gamma_T = 2.5e-10; % Killing efficiency from the CAR-T cells against GBM
omega_G = 0.0068; % Glial cell proliferation rate
psi = 2.8e-12; % Interaction coefficient between tumor cells and glial cells
psi_g = 2.571e-14; % Competition coefficient between tumor cells and glial cells
gamma_g = 2.5e-10; % Killing efficiency from the CAR-T cells against glial cells

% Initial conditions
c = zeros(size,1, 'double'); c(1) = 5e8;
t = zeros(size,1, 'double'); t(1) = 0.1*k;
h = zeros(size,1, 'double'); h(1) = 0.1*t(1);
g = zeros(size,1, 'double'); g(1) = k-h(1);
n = zeros(size,1, 'double'); n(1) = k-t(1);
k_ = g(1); %k*0.4; % Carrying capacity of antigen-positive glial population

% Formulas/Equations (G_ is H)
c_f = @(t_0) (p_C*(t(t_0)+h(t_0)/g_T+(t(t_0)+h(t_0)))) - alpha*(c(t_0))*t(t_0) - 1/tau_C * c(t_0);
t_f = @(t_0) omega_T*t(t_0)* (1-t(t_0)/k) - psi_T * (g(t_0)+h(t_0))*t(t_0) - gamma_T*t(t_0)*c(t_0);
h_f = @(t_0) omega_G*h(t_0)*(1-h(t_0)/k_)-psi_g*h(t_0)*t(t_0) - gamma_g * h(t_0) * c(t_0);
g_f = @(t_0) omega_G*g(t_0)*(1-g(t_0)/k-k_)-psi_g*g(t_0)*t(t_0);
n_f = @(t_0) exp(1)*psi*heaviside(g(t_0)+h(t_0)); % Heaviside requires symbolic maths toolbox

% Calculations
for i=2:size
    c(i) = c_f(i-1);
    t(i) = t_f(i-1);
    h(i) = h_f(i-1);
    g(i) = g_f(i-1);
    n(i) = n_f(i-1);
end


% Plot all populations
ax = tiledlayout(3,2);
xlabel(ax, "Time (days)")
ylabel(ax, "Cell number")
ax1 = nexttile;
plot(c);
title(ax1,"CAR-T Cells")
ax1 = nexttile;
plot(t);
title(ax1,"Glioblastoma Cells")
ax1 = nexttile;
plot(h);
title(ax1,"Glial Cells with Antigen")
ax1 = nexttile;
plot(g);
title(ax1,"Glial Cells without Antigen")
ax1 = nexttile;
plot(n);
title(ax1,"Neurons")












