function xdot = modelSFE(x, p, N_stages)
    % (t, x, u, parameters)
    % Model with (F)luid, (S)olid, (T)emperature
    % Di is a function of temperature (T), the function works with numbers,
    % vectors of numbers and vectors of symbolic variables
    % Rho (Peng-Robinson) are constant numbers

    %% TODO: T_bottom (N+1) - boundary conditions

    T_u           =    p(1);
    P_u           =    p(2);
    F_u           =    p(3);

    parameters    =    p(4:end);

    nstages_index =    N_stages;     % Number of stages

    nstages       =    parameters{1};
    C0solid       =    parameters{2};     % Extractor initial concentration of extract
    V             =    parameters{3};     % Extractor volume (m3)
    epsi          =    parameters{4};     % Void bed fraction
    dp            =    parameters{5};     % Diameter of the particle (m)
    L             =    parameters{6};     % Length of the extractor (m)
    rho_s         =    parameters{7};     %
    km            =    parameters{8};
    mi            =    parameters{9};

    Di            =    parameters{44};
    Dx            =    parameters{45};
    
    Ms0           =     C0solid * V * (1-epsi);
    
    A             =     V / L;        % Cross-section of the extractor (m^2)
    rp            =     dp / 2;
    lp2           =     (rp / 3)^2;

    Z             =     Compressibility(x(2*nstages_index+1:3*nstages_index),P_u,parameters);
    RHO           =     rhoPB_Comp(x(2*nstages_index+1:3*nstages_index),P_u,Z,parameters);
    
    %DIFFUSION     = axial_diffusion(x(2*nstages_index+1:3*nstages_index),P_u,F_u,RHO,parameters);
    %DIFFUSION     = Dx(x(2*nstages_index+1:3*nstages_index) ,P_u)*1e-4;
    DIFFUSION      =    Dx * ones(nstages_index,1) * 1e-5;

    %DIofT         = Di_of_T(RHO, parameters);
    DIofT         =     Di*ones(nstages_index,1) * 1e-6;

    %KM           = km_of_T(RHO, parameters);
    %KM           = km(x(2*nstages_index+1:3*nstages_index) ,P_u);
    KM            =     km*ones(nstages_index,1);

    CP            =     SpecificHeatComp(x(2*nstages_index+1:3*nstages_index), P_u, Z,RHO, parameters);
    CPRHOCP       =     cpRHOcp_Comp(x(2*nstages_index+1:3*nstages_index),P_u,Z,RHO,CP,parameters);
    KRHOCP        =     kRHOcp_Comp(x(2*nstages_index+1:3*nstages_index),P_u,Z,RHO,CP,parameters);
    VELOCITY      =     Velocity(F_u,RHO,parameters);
    
    %%

    xdot = [

    % Concentration of extract in fluid phase | 0
    % N = 1
    - VELOCITY(1)                  * (1 / L * nstages)          * ( x(0*nstages_index+1)                   - 0 )                                                                                + ...
      DIFFUSION(1)                 * ((1 / L * nstages)^2)      * ( 0                                      - 2*x(0*nstages_index+1)                   + x(0*nstages_index+2) )                  + ...
    (1-epsi)/epsi                  * 1 / mi / lp2  * DIofT(1)   * ( x(1*nstages_index+1)                   - x(0*nstages_index+1)                                                               * ...
    (rho_s / KM(1)  ./ RHO(1) ));

    % N >= 2 =< N-1
    - VELOCITY(2:nstages_index-1)  .*  (1 / L * nstages   )    .* ( x(0*nstages_index+2:1*nstages_index-1)                                            - x(0*nstages_index+1:1*nstages_index-2) ) + ...
      DIFFUSION(2:nstages_index-1) .* ((1 / L * nstages)^2)    .* ( x(0*nstages_index+1:1*nstages_index-2) - 2*x(0*nstages_index+2:1*nstages_index-1) + x(0*nstages_index+3:1*nstages_index) )   + ...
    (1-epsi)/epsi * 1 / mi / lp2   .* DIofT(2:nstages_index-1) .* ( x(1*nstages_index+2:2*nstages_index-1) -   x(0*nstages_index+2:1*nstages_index-1)                                           .* ...
    (rho_s / KM(2:nstages_index-1) ./ RHO(2:nstages_index-1)));

    % N = Nstage
    - VELOCITY(nstages_index)     .*  (1 / L * nstages   )     .* ( x(1*nstages_index)                     -   x(1*nstages_index-1))                                                            + ...
      DIFFUSION(nstages_index)    .* ((1 / L * nstages)^2)     .* ( x(1*nstages_index-1)                   - 2*x(1*nstages_index)                     + x(1*nstages_index) )                    + ...     
    (1-epsi)/epsi * 1 / mi / lp2  .* DIofT(nstages_index)      .* ( x(2*nstages_index)                     -   x(1*nstages_index)                                                              .* ...
    (rho_s / KM(nstages_index)    ./ RHO(nstages_index)));

    % Concentration of extract in solid phase | 1
    - 1 / mi / lp2 .* DIofT .* (x(1*nstages_index+1:2*nstages_index) - x(0*nstages_index+1:1*nstages_index) .* ...
    (rho_s / KM(1:nstages_index) ./ RHO(1:nstages_index)));

    % Temperature | 2
    % N = 1
    - F_u / epsi / A * (1 / L * nstages)  * CPRHOCP(1)  * (x(2*nstages_index+1)           - T_u) + ...
    KRHOCP(1) * ((1/L*nstages)^2) * (T_u -2*x(2*nstages_index+1) + x(2*nstages_index+2));

    % N >= 2 =< N-1
    - F_u / epsi / A * (1 / L * nstages) .* CPRHOCP(2:nstages_index-1) .* (x(2*nstages_index+2:3*nstages_index-1) - x(2*nstages_index+1:3*nstages_index-2)) + ...
    KRHOCP(2:nstages_index-1) .* ((1/L*nstages)^2) .* ( x(2*nstages_index+1:3*nstages_index-2) -2*x(2*nstages_index+2:3*nstages_index-1) + x(2*nstages_index+3:3*nstages_index) );

    % N = Nstage
    - F_u / epsi / A * (1 / L * nstages) .* CPRHOCP(nstages_index) .* (x(3*nstages_index) - x(3*nstages_index-1)) + ...
    KRHOCP(nstages_index) .* ((1/L*nstages)^2) .* ( x(3*nstages_index-1) -2*x(3*nstages_index) + x(3*nstages_index) ); %!

    % 3*nstage+1 = output equation
    F_u / RHO(nstages_index) * x(nstages_index) *1e3 ;   %kg/s - > g/s
    %F_u / RHO(nstages_index) / Ms0 * x(nstages_index) * 100;

    
    ];
end