function xdot = modelSFE(x, p, mask)
    % (t, x, u, parameters)
    % Model with (F)luid, (S)olid, (T)emperature
    % Di is a function of temperature (T), the function works with numbers,
    % vectors of numbers and vectors of symbolic variables
    % Rho (Peng-Robinson) are constant numbers

    %% Load Paramters
    T_u           =    p{1};
    P_u           =    p{2};
    V_u           =    p{3};

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
    FLUID         =    x(0*nstages_index+1:1*nstages_index);
    SOLID         =    x(1*nstages_index+1:2*nstages_index);
    TEMP          =    x(2*nstages_index+1:3*nstages_index);
          
    %Properties of the fluid in the extractor
    Z             =     Compressibility(TEMP, P_u,   parameters);

    RHO           =     rhoPB_Comp(     TEMP, P_u, Z,parameters);
    
    VELOCITY      =     V_u / A;

    
    %% Thermal Properties
    CP            =     SpecificHeatComp(TEMP, P_u, Z, RHO,                 parameters);            % [kJ/kg/K]
    CPRHOCP       =     cpRHOcp_Comp(    TEMP, P_u, Z, RHO, CP, epsi.*mask, parameters);
    KRHOCP        =     kRHOcp_Comp(     TEMP, P_u, Z, RHO, CP, epsi.*mask, parameters);
    
    %% BC
    Cf_0   = 0;
    Cf_B   = FLUID(nstages_index);

    T_0    = T_u;
    T_B    = TEMP(nstages_index);

    %% Derivatives
    dCf     =         FLUID                                          - [ Cf_0;   FLUID(1:nstages_index-1)           ];
    d2Cf    = [Cf_0;  FLUID(1:nstages_index-1)   ]    - 2*FLUID      + [         FLUID(2:nstages_index)    ; Cf_B   ];

    dT      =         TEMP                                           - [ T_0;    TEMP(1:nstages_index-1)            ];
    d2T     = [T_0;   TEMP(1:nstages_index-1)    ]    - 2*TEMP       + [         TEMP(2:nstages_index)     ; T_B    ];

    dz      = L/nstages_index;
    dz2     = dz^2;
   
    %%

    xdot = [
    
    %--------------------------------------------------------------------
    % Concentration of extract in fluid phase | 0
    
    - VELOCITY      ./  ( 1 - epsi .* mask ) .* dCf  ./ dz  + ...
      Dx            ./  ( 1 - epsi .* mask ) .* d2Cf ./ dz2 +...
    (epsi.*mask)    ./  ( 1 - epsi .* mask ) .* (1 ./ mi ./ lp2 .* Di)  .* ( SOLID - FLUID .* (rho_s ./ km ./ RHO ) );
    %zeros(nstages_index,1);

    %--------------------------------------------------------------------
    % Concentration of extract in solid phase | 1
     -  mask                                 .* (1 ./ mi ./ lp2 .* Di)  .* ( SOLID - FLUID .* (rho_s ./ km ./ RHO ) );
    %zeros(nstages_index,1);
    
    %--------------------------------------------------------------------
    % Temperature | 2
    - VELOCITY      ./ ( 1 - epsi .* mask )  .* CPRHOCP   .*   dT ./ dz  +   KRHOCP  .* d2T ./ dz2;
    %zeros(nstages_index,1);

    %--------------------------------------------------------------------
    % 3*nstage+1 = output equation
    V_u * FLUID(nstages_index) * 1e3 ;   %kg/s - > g/s
    
    ];

end