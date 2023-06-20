% Impementation of the ODE Model from Santurio and Barros

nmax = 120;

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
c0 = 5e8;
t0 = 0.1*k;
h0 = 0.1*t0;
g0 = k-h0;
n0 = k-t0;
k_ = 0.1*g0; %k*0.4; % Carrying capacity of antigen-positive glial population

y0 = [c0 t0 h0 g0 n0];   
t = linspace(0,120, 100*nmax);

% Runge-Kutta 4th-Order Algorithm
y_Kutta = zeros(nmax, 5);
y_Kutta(1, :) = y0;
h = t(2)-t(1);  % Constant time step
modelfcn = @(t,y) (odefcn(t, y, p_C, g_T, tau_C, alpha, omega_T, k, k_, psi_T, gamma_T, omega_G, psi, psi_g, gamma_g));

for i = 2:length(t)
    k0 = modelfcn(t(i-1),       y_Kutta(i-1,:)).';
    k1 = modelfcn(t(i-1) + h/2, y_Kutta(i-1,:)+k0*h/2).';
    k2 = modelfcn(t(i-1) + h/2, y_Kutta(i-1,:)+k1*h/2).';
    k3 = modelfcn(t(i-1) + h,   y_Kutta(i-1,:)+h*k2).';
    y_Kutta(i,:) = y_Kutta(i-1,:) + h*(k0+2*k1+2*k2+k3)/6;
end

aend =length(t);
ax = tiledlayout(3,2);
xlabel(ax, "Time (days)")
ylabel(ax, "Cell number")
ax1 = nexttile;
plot(t(1:aend), y_Kutta(1:aend,1));
title(ax1,"CAR-T Cells")
ax1 = nexttile;
plot(t(1:aend), y_Kutta(1:aend,2));
title(ax1,"Glioblastoma Cells")
ax1 = nexttile;
plot(t(1:aend), y_Kutta(1:aend,3));
title(ax1,"Glial Cells with Antigen")
ax1 = nexttile;
plot(t(1:aend), y_Kutta(1:aend,4));
title(ax1,"Glial Cells without Antigen")
ax1 = nexttile;
plot(t(1:aend), y_Kutta(1:aend,5));
title(ax1,"Neurons")

