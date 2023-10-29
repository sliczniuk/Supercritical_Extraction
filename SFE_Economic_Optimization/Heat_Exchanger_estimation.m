function [T_out_h, T_out_c, L] = Heat_Exchanger_estimation()

    % epsilon-NTU method is used to determin the outlet temperatures if the
    % heat transfer area, inner daimeter of the inner diamter of tube, 
    % outer daimeter of the inner diamter of tube, and inner daimeter of
    % the outer diamter of shell are known as well as properties of both inlet
    % streams (temperature, flow rate and heat capacity). Later the heat
    % exchnager cost is estimated.

    %% Heat exchnager properties
    D_i     = 4.094e-2;       % m - inner diamter of tube
    D_io    = 4.830e-2;       % m - outer diamter of tube
    D_o     = 7.500e-2;       % m - inner diamter of shell

    A       = 8.7;
    k       = 60;             % W/m/C - conductivity of the tube

    %% geometry of heat exchanger
    % cross-sectional area of inner tube
    S_i     = (pi/4) * D_i^2;
    % mean hydraulic diamter of inner pipe
    Dei     = D_i;

    % cross-sectional area of annulus
    S_o     = (pi/4) * (D_o^2 - D_io^2);
    % mean hydraulic diamter of annlulus
    Deo     = D_o - D_io;

    %% cold fluid - CO2 - outside
    F_c     = 5;        % kg/s
    
    T_c_in  = 15;       % C

    CP_c    = 1622;     % J/kg/C - specific heat
    k_c     = 0.138;    % W/m/C - conductivity

    rho_c   = 1044;    % kg/m3
    mu_c    = 2.7e-3;  % (N*s)/m^2
    Pr_c    = CP_c * mu_c / k_c ;

    V_c     = F_c / (rho_c * S_o);

    Re_c    = V_c * Deo / (mu_c/rho_c);

    Nu_c    = 0.023 * Re_c^0.8 * Pr_c^0.4;
    h_c     = Nu_c * k_c / Deo;     % W/m2/C

    C_c     = CP_c * F_c;

    %% warm fluid - inside
    F_w     = 4.83;    % kg/s
    T_w_in  = 95;      % C

    CP_w    = 4197;     % J/kg/C
    k_w     = 0.676;    % W/m/C

    rho_w   = 969.0;    % kg/m3
    mu_w    = 3.11e-4;    % (N*s)/m^2
    Pr_w    = CP_w * mu_w / k_w;

    V_w     = F_w / (rho_w * S_i);

    Re_w    = V_w * Dei / (mu_w/rho_w);

    Nu_w    = 0.023 * Re_w^0.8 * Pr_w^0.4;
    h_w     = Nu_w * k_w / Dei;     % W/m2/C

    C_w     = CP_w * F_w;

    %% Heat transfer coefficient

    U_inv   = 1/h_c + D_io/(2*k)*log(D_io/D_i) + 1/h_w*(D_io/D_i);
    U       = 1/U_inv;
    
    %% 
    C_min   = min(C_c, C_w);
    C_max   = max(C_c, C_w);
    C_ratio = C_min./C_max;

    %% HX properties
    % NTU
    NTU     = U.*A./C_min;

    % parallel flow HX
    %epsilon = (1 - exp(-NTU.*(1+C_ratio)) ) ./ (1 + C_ratio);

    % counter flow HX
    epsilon = (1 - exp(-NTU .* (1-C_ratio)) ) ./ (1-C_ratio.*exp(-NTU .* (1-C_ratio)) );

    % 1-shell and 2-tube HX
    %N = 1;
    %epsilon_1 = 2./ ( 1 + C_ratio + sqrt(1+C_ratio.^2) .* (1+exp(-NTU.*sqrt(1+C_ratio.^2))) ./ (1-exp(-NTU.*sqrt(1+C_ratio.^2))) );
    %epsilon   = ( ( (1-epsilon_1.*C_ratio)./(1-epsilon_1) ).^N - 1) ./ ( ( (1-epsilon_1.*C_ratio)./(1-epsilon_1) ).^N - C_ratio) ;

    % Heat flow
    Q       = epsilon .* C_min .* (T_w_in - T_c_in);

    %% Outlet temperatures
    T_w_out = T_w_in - Q./C_w;
    T_c_out = T_c_in + Q./C_c;

    %% 
    L       = A/(pi * D_io);
end