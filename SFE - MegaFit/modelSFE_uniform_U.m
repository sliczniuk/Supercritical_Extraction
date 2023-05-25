function xdot = modelSFE_uniform_U(x, p, mask)
    % (t, x, u, parameters)
    % Model with (F)luid, (S)olid, (T)emperature
    % Di is a function of temperature (T), the function works with numbers,
    % vectors of numbers and vectors of symbolic variables
    % Rho (Peng-Robinson) are constant numbers

    %% Load Paramters
    T_u           =     p{1};
    P_u           =     p{2};
    F_u           =     p{3};

    parameters    =     p(4:end);

    nstages       =     parameters{1};
    C0solid       =     parameters{2};     % Extractor initial concentration of extract
    r             =     parameters{3};     % Extractor length (m)
    epsi          =     parameters{4};     % Void bed fraction
    dp            =     parameters{5};     % Diameter of the particle (m)
    L             =     parameters{6};     % Length of the extractor (m)
    rho_s         =     parameters{7};     %
    km            =     parameters{8};
    mi            =     parameters{9};

    Di            =     parameters{44};      Di = Di * 1e-12;
    Dx            =     parameters{45};      Dx = Dx * 1e-5;
    C_SAT         =     parameters{47};

    nstages_index =     numel(mask);
    
    %% Properties of the bed
    A             =     pi*r^2 ;       % Cross-section of the extractor (m^2)
    rp            =     dp / 2;
    lp2           =     (rp / 3)^2;

    %% States
    FLUID         =     x(0*nstages_index+1:1*nstages_index);
    SOLID         =     x(1*nstages_index+1:2*nstages_index);
    TEMP          =     x(2*nstages_index+1:3*nstages_index);
    PRESSURE      =     P_u;
      
    %Properties of the fluid in the extractor
    Z             =     Compressibility(TEMP, PRESSURE,    parameters);

    RHO           =     rhoPB_Comp(     TEMP, PRESSURE, Z, parameters);
    
    VELOCITY      =     Velocity(F_u, RHO, parameters);
    
    %% Thermal Properties
    CP            =     SpecificHeatComp(TEMP, PRESSURE, Z, RHO,                 parameters);            % [kJ/kg/K]
    CPRHOCP       =     cpRHOcp_Comp(    TEMP, PRESSURE, Z, RHO, CP, epsi.*mask, parameters);
    KRHOCP        =     kRHOcp_Comp(     TEMP, PRESSURE, Z, RHO, CP, epsi.*mask, parameters);

    %% Saturation
    Sat_coe       =     Saturation_Concentration(FLUID, C_SAT); % Inverse logistic is used to control saturation. Close to saturation point, the Sat_coe goes to zero.

    %% BC
    Cf_0          =     if_else(F_u == 0, FLUID(1), 0);
    Cf_B          =     FLUID(nstages_index);

    T_0           =     T_u;
    T_B           =     TEMP(nstages_index);

    Z_0           =     Compressibility(T_u, P_u,     parameters);
    Z_B           =     Compressibility(T_B, P_u,     parameters);

    rho_0         =     rhoPB_Comp(     T_u, P_u, Z_0, parameters);
    rho_B         =     rhoPB_Comp(     T_B, P_u, Z_B, parameters);

    u_0           =     Velocity(F_u, rho_0, parameters);
    u_B           =     Velocity(F_u, rho_B, parameters);
    
    %% Derivatives

    dz            = L/nstages_index;
    
    dCfdz         = backward_diff_1_order(FLUID,Cf_0    , [],   dz);
    d2Cfdz2       = central_diff_2_order(FLUID, FLUID(1), Cf_B, dz);
    
    dTdz          = backward_diff_1_order(TEMP,T_0, [] , dz);
    d2Tdz2        = central_diff_2_order(TEMP, T_0, T_B, dz);
        
    dudz          = backward_diff_1_order(VELOCITY, u_0, [], dz);

%    d_cons_CF_dz  = central_diff_1_order(VELOCITY .* FLUID, u_0 .* Cf_0, u_B .* FLUID(nstages_index), dz);
    d_cons_CF_dz  = backward_diff_1_order(VELOCITY .* FLUID, u_0 .* Cf_0, [], dz);
   
    %%

    xdot = [
    
    %--------------------------------------------------------------------
    % Concentration of extract in fluid phase | 0
   - 1 ./  ( 1 - epsi .* mask )            .* d_cons_CF_dz                                                                     + ...
    Dx            ./  ( 1 - epsi .* mask ) .* d2Cfdz2                                                                          + ...
    (epsi.*mask)  ./  ( 1 - epsi .* mask ) .* (1 ./ mi ./ lp2 .* Di)  .* Sat_coe .* ( SOLID - FLUID .* (rho_s ./ km ./ RHO ) ) ;
    %zeros(nstages_index,1);

    %--------------------------------------------------------------------
    % Concentration of extract in solid phase | 1
    -  mask                                 .* (1 ./ mi ./ lp2 .* Di)  .* Sat_coe .* ( SOLID - FLUID .* (rho_s ./ km ./ RHO ) );
    %zeros(nstages_index,1);
    
    %--------------------------------------------------------------------
    % Temperature | 2
    - VELOCITY      ./ ( 1 - epsi .* mask )  .* CPRHOCP   .*   dTdz  +   0  .* d2Tdz2;
    
    %--------------------------------------------------------------------
    % 5*nstage+1 = output equation
    %VELOCITY(nstages_index) * A * FLUID(nstages_index) * 1e3 ;   %kg/s - > g/s
    F_u ./ RHO(nstages_index) .* FLUID(nstages_index) * 1e3;
    
    ];

end