function xdot = modelSFE(x, p, mask)
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

    FLUID         =    x(0*nstages_index+1:1*nstages_index);
    SOLID         =    x(1*nstages_index+1:2*nstages_index);
    TEMP          =    x(2*nstages_index+1:3*nstages_index);
    RHO_NS        =    x(3*nstages_index+1:4*nstages_index);
    VELOCITY_NS   =    x(4*nstages_index+1:5*nstages_index);
    
    %PRESSURE      =    P_u * ones(nstages_index,1);
    PRESSURE      = Pressure_PR(TEMP,RHO_NS,parameters);
    P_u           = PRESSURE(1);
    
    %% Properties 
    A             =     pi*r^2 ;       % Cross-section of the extractor (m^2)
    rp            =     dp / 2;
    lp2           =     (rp / 3)^2;

    % Properties of the fluid at the inlet
    Z_0           =     Compressibility(T_u, P_u,     parameters);
    rho_0         =     rhoPB_Comp(     T_u, P_u, Z_0, parameters);
    epsi_0        =     epsi*mask(1);
    epsi_B        =     epsi*mask(nstages_index);

    v_0           =     Velocity(F_u, rho_0, parameters);
    

    %Properties of the fluid in the extractor
    Z             =     Compressibility(TEMP, PRESSURE,   parameters);

    %RHO           =     rhoPB_Comp(     TEMP, PRESSURE, Z,parameters);
    RHO           = RHO_NS;

    %VELOCITY      =     Velocity(F_u, RHO, parameters);
    VELOCITY      =     VELOCITY_NS;
    v_B           =     VELOCITY(nstages_index);

    MU            = 0*ones(nstages_index,1);
    
    %% Thermal Properties
    CP            =     SpecificHeatComp(TEMP, PRESSURE, Z, RHO,                 parameters);            % [kJ/kg/K]
    CPRHOCP       =     cpRHOcp_Comp(    TEMP, PRESSURE, Z, RHO, CP, epsi.*mask, parameters);
    KRHOCP        =     kRHOcp_Comp(     TEMP, PRESSURE, Z, RHO, CP, epsi.*mask, parameters);

    %%

    xdot = [
    %{
    %--------------------------------------------------------------------
    % Concentration of extract in fluid phase | 0
    % N = 1
    - VELOCITY(1)                   ./  ( 1 - epsi .* mask(1) )                                             .*  (1 / L * nstages)       .* ( FLUID(1)                                               - 0                        ) + ...
      Dx                            ./  ( 1 - epsi .* mask(1) )                                             .* ((1 / L * nstages)^2)    .* ( FLUID(1)                 - 2*FLUID(1)                  + FLUID(2)                 ) + ...
    (epsi.*mask(1))                 ./  ( 1 - epsi .* mask(1) )                                             .* 1 / mi / lp2  * Di       .* ( SOLID(1)                                               - FLUID(1)                   * ...
    (rho_s / km ./ RHO(1) ));

    % N >= 2 =< N-1
    - VELOCITY(2:nstages_index-1)   ./  ( 1 - epsi .* mask(2:nstages_index-1) )                             .*  (1 / L * nstages   )    .* ( FLUID(2:nstages_index-1)                               - FLUID(1:nstages_index-2) ) + ...
      Dx                            ./  ( 1 - epsi .* mask(2:nstages_index-1) )                             .* ((1 / L * nstages)^2)    .* ( FLUID(1:nstages_index-2) - 2*FLUID(2:nstages_index-1)  + FLUID(3:nstages_index  ) ) + ...
    (epsi.*mask(2:nstages_index-1)) ./  ( 1 - epsi .* mask(2:nstages_index-1) )                             .* 1 ./ mi ./ lp2 .* Di     .* ( SOLID(2:nstages_index-1)                               - FLUID(2:nstages_index-1)  .* ...
    (rho_s ./ km ./ RHO(2:nstages_index-1)));                            

    % N = Nstage
    - VELOCITY(nstages_index)       ./  ( 1 - epsi .* mask(nstages_index) )                                 .*  (1 / L * nstages   )    .* ( FLUID(nstages_index)                                   - FLUID(nstages_index-1)   ) + ...
      Dx                            ./  ( 1 - epsi .* mask(nstages_index) )                                 .* ((1 / L * nstages)^2)    .* ( FLUID(nstages_index-1)   - 2*FLUID(nstages_index)      + FLUID(nstages_index)     ) + ...     
    (epsi.*mask(nstages_index))     ./  ( 1 - epsi .* mask(nstages_index) )                                 .* 1 ./ mi ./ lp2 .* Di     .* ( SOLID(nstages_index)                                   - FLUID(nstages_index  )    .* ...
    (rho_s / km ./ RHO(nstages_index)));

    %--------------------------------------------------------------------
    % Concentration of extract in solid phase | 1
     -  mask                                                                                                .* 1 ./ mi ./ lp2 .* Di     .* ( SOLID                                                  - FLUID .* (rho_s ./ km ./ RHO ) );
    %}
    zeros(2*nstages_index,1);
    %--------------------------------------------------------------------
    % Temperature | 2
    % N = 1
    - VELOCITY(1)                   ./ ( 1 - epsi.*mask(1))                 .* CPRHOCP(1)                   .*  (1 / L * nstages)       .* ( TEMP(1)                                                - T_u                      )  + ...
                                                                                KRHOCP(1)                   .* ((1 / L * nstages)^2)    .* ( T_u                     -2*TEMP(1)                     + TEMP(2)                  );

    % N >= 2 =< N-1
    - VELOCITY(2:nstages_index-1)   ./ ( 1 - epsi.*mask(2:nstages_index-1)) .* CPRHOCP(2:nstages_index-1)   .*  (1 / L * nstages)       .* ( TEMP(2:nstages_index-1)                                - TEMP(1:nstages_index-2)  )   + ...
                                                                                KRHOCP(2:nstages_index-1)   .* ((1 / L * nstages)^2)    .* ( TEMP(1:nstages_index-2) - 2*TEMP(2:nstages_index-1)    + TEMP(3:nstages_index  )  );

    % N = Nstage
   -  VELOCITY(nstages_index)       ./ ( 1 - epsi.*mask(nstages_index))     .* CPRHOCP(nstages_index)       .*  (1 / L * nstages)       .* ( TEMP(nstages_index)                                    - TEMP(nstages_index-1)    ) + ...
                                                                                KRHOCP(nstages_index)       .* ((1 / L * nstages)^2)    .* ( TEMP(nstages_index-1)   - 2*TEMP(nstages_index)        + TEMP(nstages_index)      ); %!

    %--------------------------------------------------------------------
    % Continuity - 3
    % N = 1
    - VELOCITY(1)                   ./  ( 1 - epsi.*mask(1))                                                .* (1 / L * nstages)        .* (RHO(1)                                                  - rho_0)                                    +...
    - RHO(1)                        ./  ( 1 - epsi.*mask(1))                                                .* (1 / L * nstages)        .* (VELOCITY(1)                                             - v_0)                                      +...
    - RHO(1)                        .*  VELOCITY(1)                                                         .* (1 / L * nstages)        .* ( (1 - epsi.*mask(1)).^(-1)                              - (1 - epsi_0).^(-1)                        );
                                                                                
    % N = 2:Nstage
    - VELOCITY(2:nstages_index)     ./  ( 1 - epsi.*mask(2:nstages_index))                                  .* (1 / L * nstages)        .* (RHO(2:nstages_index)                                    - RHO(1:nstages_index-1))                   +...
    - RHO(2:nstages_index)          ./  ( 1 - epsi.*mask(2:nstages_index))                                  .* (1 / L * nstages)        .* (VELOCITY(2:nstages_index)                               - VELOCITY(1:nstages_index-1))              +...
    - RHO(2:nstages_index)          .*  VELOCITY(2:nstages_index)                                           .* (1 / L * nstages)        .* ( (1 - epsi.*mask(2:nstages_index)).^(-1)                - (1 - epsi.*mask(1:nstages_index-1)).^(-1) );
                
    %--------------------------------------------------------------------
    % Momentum - 4
    %% N = 1
    ( 1 - epsi.*mask(1) ) ./ RHO(1)                                                                                                                     .*...
    ( - (1 / L * nstages) .* ( PRESSURE(1) - P_u )                                                                                                       -...
    ( VELOCITY(1)./(1 - epsi.*mask(1)) )                                                                                                                .*...
           ( ( 1./(1 - epsi.*mask(1)) ) .* (1 / L * nstages)        .* ( VELOCITY(1)                                             - v_0  )                +...
                            VELOCITY(1) .* (1 / L * nstages)        .* ( (1 - epsi.*mask(1)).^(-1)                               - (1 - epsi_0).^(-1) )) +...
    4/3*MU(1) .* ...
            ( 1./  ( 1 - epsi.*mask(1)) .*((1 / L * nstages)^2)     .* ( v_0                        -2*VELOCITY(1)               + VELOCITY(2)        )  +...
                                       2.* (1 / L * nstages)        .* (VELOCITY(1)                                              - v_0)                 .* ...
                                           (1 / L * nstages)        .* ( (1 - epsi.*mask(1)).^(-1)                               - (1 - epsi_0).^(-1) )  +...
                            VELOCITY(1) .*((1 / L * nstages)^2)     .* ( (1 - epsi_0).^(-1)         -2*(1 - epsi.*mask(1)).^(-1) + (1 - epsi.*mask(2)).^(-1) ) ));

    % N >= 2 =< N-1            
    ( 1 - epsi.*mask(2:nstages_index-1) ) ./ RHO(2:nstages_index-1)                                                                                                                                                                                            .*...
    ( - (1 / L * nstages) .* ( PRESSURE(2:nstages_index-1) - PRESSURE(1:nstages_index-2) )                                                                                                                                                                      -...
    ( VELOCITY(2:nstages_index-1)./(1 - epsi.*mask(2:nstages_index-1)) )                                                                                                                                                                                       .*...
            ( ( 1./(1 - epsi.*mask(2:nstages_index-1)) )    .* (1 / L * nstages)        .* ( VELOCITY(2:nstages_index-1)                                                                                      - VELOCITY(1:nstages_index-2)                )    +...
                        VELOCITY(2:nstages_index-1)         .* (1 / L * nstages)        .* ( (1 - epsi.*mask(2:nstages_index-1)).^(-1)                                                                        - (1 - epsi.*mask(1:nstages_index-2)).^(-1) ) )   +...
    4/3*MU(2:nstages_index-1) .* ...
                ( 1./  ( 1 - epsi.*mask(2:nstages_index-1)) .*((1 / L * nstages)^2)     .* ( VELOCITY(1:nstages_index-2)                     -2*VELOCITY(2:nstages_index-1)                                   + VELOCITY(3:nstages_index) )  +...
                                                           2.* (1 / L * nstages)        .* (VELOCITY(2:nstages_index-1)                                                                                       - VELOCITY(1:nstages_index-2)) .* ...
                                                               (1 / L * nstages)        .* ( (1 - epsi.*mask(2:nstages_index-1)).^(-1)                                                                        - (1 - epsi.*mask(1:nstages_index-2) ).^(-1) ) +...
                                VELOCITY(2:nstages_index-1) .*((1 / L * nstages)^2)     .* ( (1 - epsi.*mask(1:nstages_index-2) ).^(-1)      -2*(1 - epsi.*mask(2:nstages_index-1)).^(-1)                     + (1 - epsi.*mask(3:nstages_index)).^(-1) ) ));

    % N = Nstage
    ( 1 - epsi.*mask(nstages_index) ) ./ RHO(nstages_index)                                                                                                                                                                                       .*...
    ( - (1 / L * nstages) .* ( PRESSURE(nstages_index) - PRESSURE(nstages_index-1) )                                                                                                                                                               -...
    ( VELOCITY(nstages_index)./(1 - epsi.*mask(nstages_index)) )                                                                                                                                                                                 .* ...
                                                                                      ( ( 1./(1 - epsi.*mask(nstages_index)) )    .* (1 / L * nstages)   .* ( VELOCITY(nstages_index)               - VELOCITY(nstages_index-1)                )   +...
                                                                                                  VELOCITY(nstages_index)         .* (1 / L * nstages)   .* ( (1 - epsi.*mask(nstages_index)).^(-1) - (1 - epsi.*mask(nstages_index)).^(-1) ) )    +...
    4/3*MU(nstages_index) .* ...
            ( 1./  ( 1 - epsi.*mask(nstages_index)) .* ((1 / L * nstages)^2)    .* ( v_B                                    -2*VELOCITY(nstages_index)                  + VELOCITY(nstages_index-1)        )                                       +...
                                                   2.* (1 / L * nstages)        .* (VELOCITY(nstages_index)                                                             - VELOCITY(nstages_index-1) )                                              * ...
                                                       (1 / L * nstages)        .* ( (1 - epsi.*mask(nstages_index)).^(-1)                                              - (1 - epsi.*mask(nstages_index-1)).^(-1) )                                 +...
                           VELOCITY(nstages_index) .* ((1 / L * nstages)^2)     .* ( (1 - epsi_B).^(-1)                     -2*(1 - epsi.*mask(nstages_index)).^(-1)    + (1 - epsi.*mask(nstages_index-1)).^(-1) ) ))
    
    %--------------------------------------------------------------------
    % 5*nstage+1 = output equation
    %VELOCITY(nstages_index) * A * FLUID(nstages_index) * 1e3 ;   %kg/s - > g/s
    0;
    
    ];

 %   keyboard
end