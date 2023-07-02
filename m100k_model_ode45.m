


% Impementation of the ODE Model from Santurio and Barros

nmax = 500;
% k_?
% Parameters (from supplementary table 1)
p_C = 0.9; % CAR-T cell proliferation rate
g_T = 1e10; % T cell concentration for half-maximal CAR-T cell proliferation
tau_C = 7; % CAR-T cell lifespan
alpha = 1e-11; % Tumor cell inactivation rate###
omega_T = 0.012; % Glioblastoma proiliferation rate
k = 8.5e11; % Carrying capacity
psi_T = 2.571e-15; 
gamma_T = 2.5e-10; % Killing efficiency from the CAR-T cells against GBM
omega_G = 0.0068; % Glial cell proliferation rate
psi = 2.8e-12; % Interaction coefficient between tumor cells and glial cells
psi_g = 2.571e-14; % Competition coefficient between tumor cells and glial cells
gamma_g = 2.5e-10; % Killing efficiency from the CAR-T cells against glial cells####

% Initial conditions

% t0 = 0.100*k;
t0 = 0.025*k;
h0 = 0.1*t0;
k_ = 0.001*g0;%switch to 0.0008*g0
g0 = k-h0-t0;
n0 = k-t0;
 %k*0.4; % Carrying capacity of antigen-positive glial population
h = 1;


%c0=5e7
c0 = 5e7;
y0 = [c0 t0 h0 g0 n0];   
% Runge-Kutta 4th-Order Algorithm
y = zeros(nmax, 5);
y(1, :) = y0;
modelfcn = @(t,y) (odefcn(t, y, p_C, g_T, tau_C, alpha, omega_T, k, k_, psi_T, gamma_T, omega_G, psi, psi_g, gamma_g, h));

[t,y] = ode45(modelfcn,[0 nmax],y0);

ax = tiledlayout(3,2);
title("Solution with build-in matlab-ODE function")


ax1 = nexttile;
plot(y(:,1),LineWidth=1);

hold on;
xticks(ax1,[0:20:120])
xlim([0 120])

title(ax1,"CAR-T Cells")

ax2 = nexttile;
plot(y(:,2),LineWidth=1);

xticks(ax2,[0:20:120]);
xlim(ax2,[0 120]);
hold on;
title(ax2,"Glioblastoma Cells")

ax3 = nexttile;
plot(y(:,3),LineWidth=1);

xticks(ax3,[0:20:120]);
xlim(ax3,[0 120]);
hold on;
title(ax3,"Glial Cells with Antigen")

ax4 = nexttile;
plot(y(:,4),LineWidth=1);

xticks(ax4,[0:20:120]);
xlim(ax4,[0 120]);
hold on;
title(ax4,"Glial Cells without Antigen")

ax5 = nexttile;
plot(y(:,5),LineWidth=1);

xticks(ax5,[0:20:120])
xlim(ax5,[0 120])
hold on;
title(ax5,"Neurons")


% ##########################################
%c0=1e8
c0 = 1e8;

y0 = [c0 t0 h0 g0 n0];   
% Runge-Kutta 4th-Order Algorithm
y = zeros(nmax, 5);
y(1, :) = y0;
modelfcn = @(t,y) (odefcn(t, y, p_C, g_T, tau_C, alpha, omega_T, k, k_, psi_T, gamma_T, omega_G, psi, psi_g, gamma_g, h));
[t,y] = ode45(modelfcn,[0 nmax],y0);
plot(ax1,y(:,1),LineWidth=1);
plot(ax2,y(:,2),LineWidth=1);
plot(ax3,y(:,3),LineWidth=1);
plot(ax4,y(:,4),LineWidth=1);
plot(ax5,y(:,5),LineWidth=1);
%#########################
%c0=5e8
c0 = 5e8;

y0 = [c0 t0 h0 g0 n0];   
% Runge-Kutta 4th-Order Algorithm
y = zeros(nmax, 5);
y(1, :) = y0;
modelfcn = @(t,y) (odefcn(t, y, p_C, g_T, tau_C, alpha, omega_T, k, k_, psi_T, gamma_T, omega_G, psi, psi_g, gamma_g, h));
[t,y] = ode45(modelfcn,[0 nmax],y0);
plot(ax1,y(:,1),LineWidth=1);
plot(ax2,y(:,2),LineWidth=1);
plot(ax3,y(:,3),LineWidth=1);
plot(ax4,y(:,4),LineWidth=1);
plot(ax5,y(:,5),LineWidth=1);


%##########################
%c0=1e9
c0=1e9;
y0 = [c0 t0 h0 g0 n0];   
% Runge-Kutta 4th-Order Algorithm
y = zeros(nmax, 5);
y(1, :) = y0;
modelfcn = @(t,y) (odefcn(t, y, p_C, g_T, tau_C, alpha, omega_T, k, k_, psi_T, gamma_T, omega_G, psi, psi_g, gamma_g, h));
[t,y] = ode45(modelfcn,[0 nmax],y0);
plot(ax1,y(:,1),LineWidth=1);
plot(ax2,y(:,2),LineWidth=1);
plot(ax3,y(:,3),LineWidth=1);
plot(ax4,y(:,4),LineWidth=1);
plot(ax5,y(:,5),LineWidth=1);
%plot legend
legend(ax1,"C(0)=5e7 Cells","C(0)=1e8 Cells","C(0)=5e8 Cells","C(0)=1e9 Cells")
legend(ax2,"C(0)=5e7 Cells","C(0)=1e8 Cells","C(0)=5e8 Cells","C(0)=1e9 Cells")
legend(ax3,"C(0)=5e7 Cells","C(0)=1e8 Cells","C(0)=5e8 Cells","C(0)=1e9 Cells")
legend(ax4,"C(0)=5e7 Cells","C(0)=1e8 Cells","C(0)=5e8 Cells","C(0)=1e9 Cells")
legend(ax5,"C(0)=5e7 Cells","C(0)=1e8 Cells","C(0)=5e8 Cells","C(0)=1e9 Cells")
%plot axis title
xlabel(ax1, "Time (days)")
xlabel(ax2, "Time (days)")
xlabel(ax3, "Time (days)")
xlabel(ax4, "Time (days)")
xlabel(ax5, "Time (days)")
ylabel(ax1, "Cell number")
ylabel(ax2, "Cell number")
ylabel(ax3, "Cell number")
ylabel(ax4, "Cell number")
ylabel(ax5, "Cell number")