function dydt = odefcn(t, y, p_C, g_T, tau_C, alpha, omega_T, k, k_, psi_T, gamma_T, omega_G, psi, psi_g, gamma_g) 
%ODEFCN ODE for CAR-T cell Therapy
%    Parameters (from supplementary table 1)
%    p_C        CAR-T cell proliferation rate
%    g_T  T     cell concentration for half-maximal CAR-T cell proliferation
%    tau_C      CAR-T cell lifespan
%    alpha      Tumor cell inactivation rate
%    omega_T    Glioblastoma proiliferation rate
%    k          Carrying capacity
%    k_         Carrying capacity of antigen positive cells
%    psi_T      Competition coefficient between glial cells and tumor cells
%    gamma_T    Killing efficiency from the CAR-T cells against GBM
%    omega_G    Glial cell proliferation rate
%    psi        Interaction coefficient between tumor cells and glial cells
%    psi_g      Competition coefficient between tumor cells and glial cells
%    gamma_g    Killing efficiency from the CAR-T cells against glial cells
    
    dydt = zeros(5,1);

    dydt(1) = (p_C*(y(2)+y(3)/g_T+(y(2)+y(3)))) - alpha*(y(1))*y(2) - 1/tau_C * y(1);
    dydt(2) = omega_T*y(2)* (1-y(2)/k) - psi_T * (y(4)+y(3))*y(2) - gamma_T*y(2)*y(1);
    dydt(3) = omega_G*y(3)*(1-y(3)/k_)-psi_g*y(3)*y(2) - gamma_g * y(3) * y(1);
    dydt(4) = omega_G*y(4)*(1-y(4)/k-k_)-psi_g*y(4)*y(2);
    dydt(5) = psi*(y(3)+y(4))*heaviside(-(y(3)+y(4))) * y(5); % Heaviside requires symbolic maths toolbox
end 