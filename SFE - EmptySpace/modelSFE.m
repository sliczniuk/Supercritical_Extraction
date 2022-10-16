function xdot = modelSFE(x, p, N_stage, L, L_bed)
    % (t, x, u, parameters)
    % Model with (F)luid, (S)olid, (T)emperature
    % Di is a function of temperature (T), the function works with numbers,
    % vectors of numbers and vectors of symbolic variables
    % Rho (Peng-Robinson) are constant numbers

    %% Load Paramters
    T_u           =    p(1);
    P_u           =    p(2);
    F_u           =    p(3);

    parameters    =    p(4:end);

    nstages_index =    N_stage;     % Number of stages

    nstages       =    parameters{1};
    C0solid       =    parameters{2};     % Extractor initial concentration of extract
    r             =    parameters{3};     % Extractor length (m)
    epsi          =    parameters{4};     % Void bed fraction
    dp            =    parameters{5};     % Diameter of the particle (m)
    %L             =    parameters{6};     % Length of the extractor (m)
    rho_s         =    parameters{7};     %
    km            =    parameters{8};
    mi            =    parameters{9};

    Di            =    parameters{44};
    Dx            =    parameters{45};

    FLUID         =    x(0*nstages_index+1:1*nstages_index);

    SOLID         =    x(1*nstages_index+1:2*nstages_index);
    TEMP          =    x(2*nstages_index+1:3*nstages_index);
    
    %% Properties 
    A             =     pi*r^2 ;       % Cross-section of the extractor (m^2)
    rp            =     dp / 2;
    lp2           =     (rp / 3)^2;

    alpha         =     floor(nstages_index*( L_bed / L) );
    %alpha = nstages_index;
    EPSI          =     [2/3 * ones(alpha,1); ones(nstages_index-alpha,1) ];

    Z             =     Compressibility(TEMP,P_u,parameters);
    RHO           =     rhoPB_Comp(TEMP,P_u,Z,parameters);

    VELOCITY      =     Velocity(F_u,RHO,EPSI,parameters);
    %VELOCITY_EMPTY=     Velocity(F_u,RHO,parameters)*epsi;
    
    %% Extraction Kinetic
    %DIFFUSION     = axial_diffusion(x(2*nstages_index+1:3*nstages_index),P_u,F_u,RHO,parameters);
    %DIFFUSION     = Dx(x(2*nstages_index+1:3*nstages_index) ,P_u)*1e-4;
    DIFFUSION      =    Dx * ones(nstages_index,1) * 1e-6;

    %DIofT         = Di_of_T(RHO, parameters);
    DIofT         = [Di * ones(alpha,1); zeros(nstages_index-alpha,1) ] * 1e-12;
    %DIofT         =     Di*ones(nstages_index,1) * 1e-12;

    %KM           = km_of_T(RHO, parameters);
    %KM           = km(x(2*nstages_index+1:3*nstages_index) ,P_u);
    KM            =     km*ones(nstages_index,1);

    %% Thermal Properties
    CP            =     SpecificHeatComp(TEMP, P_u, Z, RHO,           parameters);
    CPRHOCP       =     cpRHOcp_Comp(    TEMP, P_u, Z, RHO, CP, EPSI, parameters);
    KRHOCP        =     kRHOcp_Comp(     TEMP, P_u, Z, RHO, CP, EPSI, parameters);

    %%
    %Ms0           =     C0solid * V * (1-epsi);
    %M             =     F_u / RHO(nstages_index) * x(4*nstages_index);
    
    %%

    xdot = [

    % Concentration of extract in fluid phase | 0
    % N = 1
    - VELOCITY(1)                  * (1 / L * nstages)                                                  * ( FLUID(1)                                               - 0                        ) + ...
      DIFFUSION(1)                 * ((1 / L * nstages)^2)                                              * ( FLUID(1)                 - 2*FLUID(1)                  + FLUID(2)                 ) + ...
    (1-EPSI(1))/EPSI(1)            * 1 / mi / lp2  * DIofT(1)                                           * ( SOLID(1)                                               - FLUID(1)                   * ...
    (rho_s / KM(1)  ./ RHO(1) ));

    % N >= 2 =< N-1
    - VELOCITY(2:nstages_index-1)  .*  (1 / L * nstages   )                                            .* ( FLUID(2:nstages_index-1)                               - FLUID(1:nstages_index-2) ) + ...
      DIFFUSION(2:nstages_index-1) .* ((1 / L * nstages)^2)                                            .* ( FLUID(1:nstages_index-2) - 2*FLUID(2:nstages_index-1)  + FLUID(3:nstages_index  ) ) + ...
    (1-EPSI(2:nstages_index-1))./EPSI(2:nstages_index-1) .* 1 ./ mi ./ lp2 .* DIofT(2:nstages_index-1)   .* ( SOLID(2:nstages_index-1)                               - FLUID(2:nstages_index-1)  .* ...
    (rho_s ./ KM(2:nstages_index-1) ./ RHO(2:nstages_index-1)));

    % N = Nstage
    - VELOCITY(nstages_index)     .*  (1 / L * nstages   )                                             .* ( FLUID(nstages_index)                                   - FLUID(nstages_index-1)   ) + ...
      DIFFUSION(nstages_index)    .* ((1 / L * nstages)^2)                                             .* ( FLUID(nstages_index-1)   - 2*FLUID(nstages_index)      + FLUID(nstages_index  )   ) + ...     
    (1-EPSI(nstages_index))./EPSI(nstages_index)          * 1 / mi / lp2 .*  DIofT(nstages_index)       .* ( SOLID(nstages_index)                                   - FLUID(nstages_index  )    .* ...
    (rho_s / KM(nstages_index)    ./ RHO(nstages_index)));

    % Concentration of extract in solid phase | 1
    - 1 / mi / lp2 .* DIofT .* ( SOLID - FLUID .* (rho_s ./ KM ./ RHO ) );

    % Temperature | 2
    % N = 1]
    - F_u / EPSI(1)                 / A *  (1 / L * nstages)    * CPRHOCP(1)                  * ( TEMP(1)                                            - T_u                     ) + ...
      KRHOCP(1)                         * ((1 / L * nstages)^2)                               * ( T_u                     -2*TEMP(1)                 + TEMP(2)                 );

    % N >= 2 =< N-1
    - F_u / EPSI(2:nstages_index-1) / A * (1 / L * nstages)    .* CPRHOCP(2:nstages_index-1) .* ( TEMP(2:nstages_index-1)                            - TEMP(1:nstages_index-2) ) + 0%...
    %  KRHOCP(2:nstages_index-1)       .* ((1 / L * nstages)^2)                               .* ( TEMP(1:nstages_index-2) -2*TEMP(2:nstages_index-1) + TEMP(3:nstages_index)   );

    % N = Nstage
   - F_u / EPSI(nstages_index)     / A *  (1 / L * nstages)    .* CPRHOCP(nstages_index)    .* ( TEMP(nstages_index)                                - TEMP(nstages_index-1)   ) + ...
      KRHOCP(nstages_index)            .* ((1 / L * nstages)^2)                              .* ( TEMP(nstages_index-1)   -2*TEMP(nstages_index)     + TEMP(nstages_index)     ); %!
    
    % 4*nstage+1 = output equation
    F_u / RHO(nstages_index) * x(nstages_index) * 1e3 ;   %kg/s - > g/s
    %F_u / RHO(nstages_index) / Ms0 * x(nstages_index) * 100;
    %M.*(1-M./Ms0) * 1e3;

    %F_u / RHO(nstages_index) * FLUID(nstages_index) ;   %kg/kg/s - Vargas
    
    ];
end