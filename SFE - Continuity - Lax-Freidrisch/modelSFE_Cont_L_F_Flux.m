function xdot = modelSFE_Cont_L_F_Flux(x, p, mask, u_in, u_out, dt, flag_cons)
    % (t, x, u, parameters)
    % Model with (F)luid, (S)olid, (T)emperature
    % Di is a function of temperature (T), the function works with numbers,
    % vectors of numbers and vectors of symbolic variables
    % Rho (Peng-Robinson) are constant numbers

    %% Load Paramters
    T_u           =    p{1};
    P_u           =    p{2};
    F_u           =    p{3};

    parameters    =    p(4:end);

    nstages       =    parameters{1};
    C0solid       =    parameters{2};     % Extractor initial concentration of extract
    r             =    parameters{3};     % Extractor length (m)
    epsi          =    parameters{4};     % Void bed fraction
    dp            =    parameters{5};     % Diameter of the particle (m)
    L             =    parameters{6};     % Length of the extractor (m)
    rho_s         =    parameters{7};     %
    km            =    parameters{8};
    mi            =    parameters{9};

    Di            =    parameters{44};      Di = Di * 1e-12;
    Dx            =    parameters{45};      Dx = Dx * 1e-6;

    nstages_index =    numel(mask);
    
    %% Properties of the bed
    A             =     pi*r^2 ;       % Cross-section of the extractor (m^2)
    rp            =     dp / 2;
    lp2           =     (rp / 3)^2;
   

    %% States
    FLUID_NS         =    x(0*nstages_index+1:1*nstages_index);
    SOLID_NS         =    x(1*nstages_index+1:2*nstages_index);
    TEMP_NS          =    x(2*nstages_index+1:3*nstages_index);
    RHO_NS           =    x(3*nstages_index+1:4*nstages_index);
    %VELOCITY_NS      =    x(4*nstages_index+1:5*nstages_index);
    
    %E_Inv         = (1 - epsi.*mask).^(-1);

    FLUID       = FLUID_NS  ./ ( 1 - epsi .* mask );
    SOLID       = SOLID_NS;
    TEMP        = TEMP_NS;
    RHO         = RHO_NS    ./ ( 1 - epsi .* mask );
    %VELOCITY    = VELOCITY_NS ./ RHO;      
    
    %%

    %PRESSURE      =     P_u .* linspace( u_in, u_out, nstages_index)';
    PRESSURE      =     Pressure_PR(TEMP,RHO,parameters);
      
    %Properties of the fluid in the extractor
    Z             =     Compressibility(TEMP, PRESSURE,   parameters);

    %RHO           =     rhoPB_Comp(     TEMP, PRESSURE, Z,parameters);
    
    VELOCITY      =  (F_u / A) .* linspace( u_in, u_out, nstages_index)';
    %VELOCITY      =     Velocity(F_u, RHO, parameters);
    %VELOCITY      =     VELOCITY_NS;
    
    %% Thermal Properties
    CP            =     SpecificHeatComp(TEMP, PRESSURE, Z, RHO,                 parameters);            % [kJ/kg/K]
    CPRHOCP       =     cpRHOcp_Comp(    TEMP, PRESSURE, Z, RHO, CP, epsi.*mask, parameters);
    KRHOCP        =     kRHOcp_Comp(     TEMP, PRESSURE, Z, RHO, CP, epsi.*mask, parameters);

    %alpha         =     ThermalExpansion(TEMP, PRESSURE, Z, RHO,                 parameters);
    %beta          =     compressibility_beta(TEMP, PRESSURE, RHO, parameters);
    
    MU            =     Viscosity(TEMP,RHO);

    %% BC
    Cf_0   = 0;
    Cf_B   = FLUID(nstages_index);

    T_0    = T_u;
    T_B    = TEMP(nstages_index);

    Z_0    = Compressibility(T_u, P_u,      parameters);
    rho_0  = rhoPB_Comp(     T_u, P_u, Z_0, parameters);
    %rho_0  = RHO(1);
    rho_B  = RHO(nstages_index);

    %VELOCITY      =     Velocity(F_u, rho_0, parameters).*linspace(1,1,nstages_index)';

    %u_0    = Velocity(F_u, rho_0, parameters);
    u_0     = F_u / A;
    %u_0    = VELOCITY(1);
    u_B    = VELOCITY(nstages_index);

    %epsi_0 = ( 1-0 ) .^(-1) ;
    %epsi_B = ( 1-0 ) .^(-1) ;

    %P_0    = PRESSURE(1);
    %P_0    = P_u;

    %% Continuity Flux: variables
    epsi_rho   = [(1-0)*rho_0; RHO_NS; (1-0)*rho_B];

    f_cont     = [u_0 .* rho_0; VELOCITY .* RHO; u_B .* rho_B];

    %% Concentration Flux: variables
    epsi_Cf    = [(1-0)*Cf_0; FLUID_NS;  (1-0)*Cf_B];

    f_Cf       = [u_0 .* Cf_0; VELOCITY .* FLUID; u_B .* Cf_B];

    %% Concentration Flux: variables
    c   = SpeedOfSound(TEMP, PRESSURE, Z, RHO, parameters);
    c_0 = SpeedOfSound(T_0, P_u, Z_0, rho_0, parameters);
    c_B = SpeedOfSound(T_B, PRESSURE(nstages_index), Z(nstages_index), rho_B, parameters);

    %% Derivatives

    dz        = L/nstages_index;

    Flux_cont = flux_LF(epsi_rho, f_cont, dt, dz);
    %Flux_cont = flux_Ru(epsi_rho, f_cont, [u_0; VELOCITY; u_B], [c_0; c; c_B], [0; epsi .* mask; 0], dz);

    Flux_Cf   = flux_LF(epsi_Cf, f_Cf, dt, dz);
    
    dCfdz     = backward_diff_1_order(FLUID,Cf_0, [], dz);
    %dCfdz     = central_diff_1_order(FLUID, Cf_0, Cf_B, dz);
    d2Cfdz2   = central_diff_2_order(FLUID, FLUID(1), Cf_B, dz);
    
    dTdz      = backward_diff_1_order(TEMP,T_0, [], dz);
    %dTdz      = central_diff_1_order(TEMP,T_0, T_B, dz);
    d2Tdz2    = central_diff_2_order(TEMP, T_0, T_B, dz);
    
    dRhodz    = backward_diff_1_order(RHO, rho_0, [], dz);
    %dRhodz    = central_diff_1_order(RHO, rho_0, rho_B, dz);
    %d2Rhodz2   = central_diff_2_order(RHO, rho_0, rho_B, dz);
    
    dudz      = backward_diff_1_order(VELOCITY, u_0, [], dz);
    %dudz      = central_diff_1_order(VELOCITY, u_0, u_B, dz);
    %d2udz2  = central_diff_2_order(VELOCITY, u_0, u_B, dz);

    d_cons_rho_dz      = backward_diff_1_order(VELOCITY .* RHO, u_0 .* rho_0, [], dz);
    %d_cons_rho_dz     = central_diff_1_order(VELOCITY .* RHO, u_0 .* rho_0, u_B .* rho_B, dz);

    d_cons_CF_dz      = backward_diff_1_order(VELOCITY .* FLUID, u_0 .* Cf_0, [], dz);
    %d_cons_CF_dz      = central_diff_1_order(VELOCITY .* FLUID, u_0 .* Cf_0, u_B .* Cf_B, dz);

    %dEdz    = backward_diff_1_order((1 - epsi.*mask),(1 - 0),[],dz);
    %d2Edz2   = central_diff_2_order(E_Inv,epsi_0, epsi_B, dz);

    %dPdz    = forward_diff_1_order(PRESSURE,[],PRESSURE(end),dz);
    %dPdz    = backward_diff_1_order(PRESSURE,P_u,[],dz);
    %dPdz    = central_diff_1_order(PRESSURE,P_u,PRESSURE(end),dz);
    %dPdz    = dPdz .* 1e5;                                                                                                % bar = 1e5 * Pa; 
   
    %%

    %--------------------------------------------------------------------
    % Concentration of extract in fluid phase | 0
    if isequal(flag_cons,'conservative')
        C_f = - d_cons_CF_dz                                                                                                    + ...
                0            ./  ( 1 - epsi .* mask ) .* d2Cfdz2                                                                + ...
               (epsi.*mask)   .* (1 ./ mi ./ lp2 .* Di)  .* ( SOLID - FLUID .* (rho_s ./ km ./ RHO ) )                          ;
    else
        C_f = - Flux_Cf                                                                                                         + ...
                0            ./  ( 1 - epsi .* mask ) .* d2Cfdz2                                                                + ...
               (epsi.*mask)   .* (1 ./ mi ./ lp2 .* Di)  .* ( SOLID - FLUID .* (rho_s ./ km ./ RHO ) )                          ;
    end

    %--------------------------------------------------------------------
    % Concentration of extract in solid phase | 1
    C_s  = -  mask             .* (1 ./ mi ./ lp2 .* Di)  .* ( SOLID - FLUID .* (rho_s ./ km ./ RHO ) );
    %zeros(nstages_index,1);

    %--------------------------------------------------------------------
    % Temperature | 2
    Heat = - VELOCITY      ./ ( 1 - epsi .* mask )  .* CPRHOCP   .*   dTdz  +  KRHOCP  .* d2Tdz2;
    %zeros(nstages_index,1);

    %--------------------------------------------------------------------
    % Continuity - 3 
    if isequal(flag_cons,'conservative')
        Continuity = - d_cons_rho_dz     ;
    else
        Continuity = - Flux_cont    ;
    end
    %zeros(nstages_index,1);

    %--------------------------------------------------------------------
    % Momentum - 4
    %{
    - VELOCITY ./ (1-epsi.*mask)                           .* dudz              +...
    - VELOCITY .* VELOCITY                                 .* dEdz              +...
    - (1-epsi.*mask) ./ RHO                                .* dPdz              +...
    +     4/3 .* MU .* (1-epsi.*mask) .*VELOCITY ./ RHO    .* d2Edz2             +...
    + 2 * 4/3 .* MU .* (1-epsi.*mask)            ./ RHO    .* dudz .* dEdz  +...
    +     4/3 .* MU                              ./ RHO    .* d2udz2;
    %}
    if isequal(flag_cons,'conservative')
        Momentum = zeros(nstages_index,1);
    else
        %Momentum = -VELOCITY .* ( 1 - epsi .* mask ) .* dudz - VELOCITY.^2 .* dEdz - ( 1 - epsi .* mask ) ./ RHO .* dPdz ;
        Momentum = zeros(nstages_index,1);
    end
    
    %- VELOCITY .* dudz - 1./RHO .* dPdz + 4/3 .* MU ./ RHO .* d2udz2;

    %--------------------------------------------------------------------
    % 5*nstage+1 = output equation
    %yield = F_u / RHO(nstages_index) * FLUID(nstages_index) * 1e3 ;   %kg/s - > g/s
    %yield = F_u * FLUID(nstages_index) * 1e3 ;   %m3/s - > g/s
    yield = VELOCITY(nstages_index) * A * FLUID(nstages_index) * 1e3 ;   %m/s - > g/s
    

    %%
    xdot = [C_f; C_s; Heat; Continuity; Momentum; yield ];

end