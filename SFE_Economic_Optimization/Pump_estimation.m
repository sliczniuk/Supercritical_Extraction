function [T_out, W, Cost] = Pump_estimation(T_in, P_in, P_out, F, Parameters)

    % Pump_estimation function estimate the outlet temperature after
    % isentropic compression of gas entering a pump at T_in and P_in. As
    % entropy is assumed to be constnat, the S_in(T_in,P_in) =
    % S_out(T_out,P_out). P_out is design variable and T_out is unknown. An
    % optimization problem is solved to find T_out. Later the power requre
    % to compress the fluid is calculated. The power is used to evaluate
    % the equipment cost in $.
    % T_in and T_out are in [K]
    % P_in and P_out are in [bar]
    % F is in [kg/s]

    %% import casadi and set decision variable T_out
    import casadi.*
    T_out   = MX.sym('T_out',1);
    
    %% Evaluate properties of fluid at the inlet to the pump/compressor
    Z_in    = Compressibility( T_in, P_in,       Parameters );
    rho_in  = rhoPB_Comp(      T_in, P_in, Z_in, Parameters );
    S_in    = SpecificEntropy( T_in, P_in, Z_in, rho_in, Parameters); %kJ/kg
    H_in    = SpecificEnthalpy(T_in, P_in, Z_in, rho_in, Parameters); %kJ/kg

    %% Evaluate properties of fluid at the outlet to the pump/compressor - symbolic
    Z_out   = Compressibility( T_out, P_out,        Parameters );
    rho_out = rhoPB_Comp(      T_out, P_out, Z_out, Parameters );
    S_out   = SpecificEntropy( T_out, P_out, Z_out, rho_out, Parameters); %kJ/kg

    %% Set optimization problem and find S_in = S_out by adjusting T_out
    cost_function = (S_in - S_out).^2;

    g_compressor = Function('g_compressor', {T_out}, {cost_function});
    G_compressor = rootfinder('G_compressor','newton',g_compressor);
    
    %% Evaluate properties of fluid at the outlet to the pump/compressor
    T_out   = G_compressor(T_in);
    Z_out   = Compressibility( T_out, P_out,        Parameters );
    rho_out = rhoPB_Comp(      T_out, P_out, Z_out, Parameters );
    S_out   = SpecificEntropy( T_out, P_out, Z_out, rho_out, Parameters); %kJ/kg
    H_out   = SpecificEnthalpy(T_out, P_out, Z_out, rho_out, Parameters); %kJ/kg

    %% Ideal work calculations
    W = F .* (H_out - H_in); % kg/s * kJ/kg = kJ/s = kW

    %% Equipment cost 
    Cost = 10167.5*W^0.46; % $ https://sci-hub.st/10.3390/en13236454

end