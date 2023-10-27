function [T_out_h, T_out_c, A] = Heat_Exchanger_estimation()

    % epsilon-NTU method is used to determin the outlet temperatures if the
    % heat transfer area is known as well as properties of both inlet
    % streams (temperature, flow rate and heat capacity). Later the heat
    % exchnager cost is estimated.

    %% Heat exchnager properties
    D_i = 0.0525;       % m
    D_o = 0.0603;       % m
    D_a = 0.0779;       % m

    %% cold fluid
    F_c     = 5000;     % kg/h
    F_c     = F_c/3600; % kg/s

    T_c_in  = 20;       % C

    CP_c    = 4.2;      % kJ/kg/K
    k_c     = 0.687;    % W/m/K

    rho_c   = 932.53;   % kg/m3
    mu_c    = 0.207e-3; % Pa*s
    Pr_c    = 1.28;

    C_c     = CP_c * F_c;

    %% warm fluid
    F_w     = 1.36;     % kg/s
    T_w_in  = 250;      % C

    CP_w    = 4.2;      % kJ/kg/K
    k_w     = 0.687;    % W/m/K

    rho_w   = 932.53;   % kg/m3
    mu_w    = 0.207e-3; % Pa*s
    Pr_w    = 1.28;

    C_w     = CP_w * F_w;

    %% Flow cross-section
    A_i     = pi*D_i^2/4;
    A_a     = pi*(D_a^2 - D_o^2)/4;

    %% Flow-rates
    V_i     = F_w / (rho_w * A_i);
    V_a     = F_c / (rho_c * A_a);

    %% Annulus Equivalent Diameter
    D_e     = (D_a^2 - D_o^2) / D_o;

    %% Reynolds
    Re_w    = V_i * D_i / mu_w;
    Re_c    = V_a * D_e / mu_c;
    
    %% Nusselt
    Nu_w    = 0.023 * Re_w^(0.8) * Pr_w^(0.3);
    Nu_c    = 0.023 * Re_c^(0.8) * Pr_c^(0.4);

    h_w     = Nu_w * k_w / D_i;
    h_c     = Nu_c * k_c / D_o;

    %% Heat transfer coefficient

    U = 1/(h_w*A_i) + 1/(h_c*A_a);

    C_min   = min(C_c, C_w);
    C_max   = max(C_c, C_w);
    C_ratio = C_min./C_max;

    %A       = 30;                                               % m^2
    %U       = 500;                                              % W/(m^2 K)

    %% HX properties
    % NTU
    %NTU     = U.*A./C_min;

    % parallel flow HX
    %epsilon = (1 - exp(-NTU.*(1+C_ratio)) ) ./ (1 + C_ratio);

    % counter flow HX
    %epsilon = (1 - exp(-NTU .* (1+C_ratio)) ) ./ (1-C_ratio.*exp(-NTU .* (1-C_ratio)) );

    % 1-shell and 2-tube HX
    %N = 1;
    %epsilon_1 = 2./ ( 1 + C_ratio + sqrt(1+C_ratio.^2) .* (1+exp(-NTU.*sqrt(1+C_ratio.^2))) ./ (1-exp(-NTU.*sqrt(1+C_ratio.^2))) );
    %epsilon   = ( ( (1-epsilon_1.*C_ratio)./(1-epsilon_1) ).^N - 1) ./ ( ( (1-epsilon_1.*C_ratio)./(1-epsilon_1) ).^N - C_ratio) ;

    % Heat flow
    %Q       = epsilon .* C_min .* (T_in_h - T_in_c);

    %% Outlet temperatures
    %T_out_h = T_in_h - Q./C_h;
    %T_out_c = T_in_c + Q./C_c;

end