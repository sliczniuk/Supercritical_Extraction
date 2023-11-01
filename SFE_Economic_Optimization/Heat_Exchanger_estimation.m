function [T_w_out, T_c_out, Cost] = Heat_Exchanger_estimation(T_out_compressor, P, Parameters)

    % epsilon-NTU method is used to determin the outlet temperatures if the
    % heat transfer area, inner daimeter of the inner diamter of tube, 
    % outer daimeter of the inner diamter of tube, and inner daimeter of
    % the outer diamter of shell are known as well as properties of both inlet
    % streams (temperature, flow rate and heat capacity). Later the heat
    % exchnager cost is estimated.

    Parameters_table        = readtable('Parameters.csv') ;                 % Table with prameters
    Parameters              = num2cell(Parameters_table{:,3});              % Parameters within the model 

    %% Heat exchnager properties
    D_i     = 4.094e-2;       % m - inner diamter of tube
    D_io    = 4.830e-2;       % m - outer diamter of tube
    D_o     = 7.500e-2;       % m - inner diamter of shell

    A       = 2;
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

    %% fluid outside - water
    F_c     = 0.75;        % kg/s
    
    T_c_in  = 50+273;       % C

    CP_c    = 4119;     % J/kg/C - specific heat
    k_c     = 0.607;    % W/m/C - conductivity

    rho_c   = 997;    % kg/m3
    mu_c    = 1e-3;  % (N*s)/m^2
    Pr_c    = CP_c * mu_c / k_c ;

    V_c     = F_c / (rho_c * S_o);

    Re_c    = V_c * Deo / (mu_c/rho_c);

    Nu_c    = 0.023 * Re_c^0.8 * Pr_c^0.4;
    h_c     = Nu_c * k_c / Deo;     % W/m2/C

    C_c     = CP_c * F_c;

    %% fluid inside - co2
    F_w     = 0.6;    % kg/s % TODO: get the F from the simulation
    T_w_in  = T_out_compressor;      % C
    
    % properties of co2 at the highest T
    Z_in     = Compressibility( T_w_in, P,               Parameters );
    rho_in   = rhoPB_Comp(      T_w_in, P, Z_in,         Parameters );
    CP_in    = SpecificHeatComp(T_w_in, P, Z_in, rho_in, Parameters );      % [J/kg/K]
    k_in     = HeatConductivity_Comp(T_w_in, rho_in)*1e-3;                  % [ 10^-3 (W / m * K) ] -> (W / m * K)
    mu_in    = Viscosity(T_w_in, rho_in);                                   % Pa*s

    % properties of co2 at the lowest T
    Z_out    = Compressibility( T_c_in, P,                 Parameters );
    rho_out  = rhoPB_Comp(      T_c_in, P, Z_out,          Parameters );
    CP_out   = SpecificHeatComp(T_c_in, P, Z_out, rho_out, Parameters );    % [J/kg/K]
    k_out    = HeatConductivity_Comp(T_c_in, rho_out)*1e-3;                 % [ 10^-3 (W / m * K) ] -> (W / m * K)
    mu_out   = Viscosity(T_c_in,rho_out);                                   % Pa*s

    % average properties of co2
    rho_w    = (rho_in + rho_out) / 2;
    CP_w     = (CP_in + CP_out) / 2;
    mu_w     = (mu_in + mu_out) / 2;
    k_w      = (k_in + k_out) / 2;
    
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
    T_w_out = T_w_in - Q./C_w;                                              % CO2
    T_c_out = T_c_in + Q./C_c;                                              % water

    %% 
    L       = A/(pi * D_io);

    %% Cost - Product and Process Design Principles : Synthesis, Analysis, and Evaluation; https://sci-hub.ru/10.3390/en13102656
    %  P is the shell-side pressure in MPa
    %P = P/10;

    A_ft = A * 10.764;                                                      % m2 -> ft2
    P_psi = P * 14.5;

    CB   = exp( 7.2718 + 0.16*log(A_ft) );
    FP   = 0.851 + 0.1292 * (P_psi/600) + 0.0198*(P_psi/600)^2;             % Pressure factor
    FM = (600/567);                                                         % a cost index in 2013 (CE = 567) and purchase cost for CE = 600

    Cost = CB * FP;

end
