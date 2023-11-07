function xdot = modelSFE_SA(X, p, mask, dt, Sdot)
    % (t, x, u, parameters)
    % Model with (F)luid, (S)olid, (T)emperature
    % Di is a function of temperature (T), the function works with numbers,
    % vectors of numbers and vectors of symbolic variables
    % Rho (Peng-Robinson) are constant numbers

    %% Load Paramters
    Nx = numel(X)/2;
    x = X(1:Nx);
    S = X(Nx+1:end);

    T_u           =     p{1};
    P_u           =     p{2};
    F_u           =     p{3};

    parameters    =     p(4:end);

    %nstages       =     parameters{1};
    C0solid       =     parameters{2};     % Extractor initial concentration of extract
    r             =     parameters{3};     % Extractor length (m)
    epsi          =     parameters{4};     % Void bed fraction
    dp            =     parameters{5};     % Diameter of the particle (m)
    L             =     parameters{6};     % Length of the extractor (m)
    rho_s         =     parameters{7};     %
    mode          =     parameters{8};
    km            =     1e5;
    mi            =     parameters{9};

    %Di            =     parameters{44};      Di = Di * 1e-14;
    %Dx            =     parameters{45};      Dx = Dx * 1e-6;
    %SAT           =     parameters{47};
    %shape         =     parameters{48};

    nstages_index =     numel(mask);
    
    %% Properties of the bed
    A             =     pi*r^2 ;       % Cross-section of the extractor (m^2)
    rp            =     dp / 2;
    lp2           =     (rp / 3)^2;

    %% States
    FLUID         =     x(0*nstages_index+1:1*nstages_index);
    SOLID         =     x(1*nstages_index+1:2*nstages_index);
    ENTHALPY_RHO  =     x(2*nstages_index+1:3*nstages_index);
    PRESSURE      =     x(3*nstages_index+1);

    TEMP          =     Reconstruct_T_from_enthalpy(ENTHALPY_RHO, PRESSURE, parameters);
      
    %Properties of the fluid in the extractor
    Z             =     Compressibility(TEMP, PRESSURE,    parameters);

    RHO           =     rhoPB_Comp(     TEMP, PRESSURE, Z, parameters);   
    VELOCITY      =     Velocity(F_u, RHO, parameters);
    
    %% Thermal Properties
    CP            =     SpecificHeatComp(TEMP, PRESSURE, Z, RHO,                 parameters);            % [kJ/kg/K]
    CPRHOCP       =     cpRHOcp_Comp(    TEMP, PRESSURE, Z, RHO, CP, epsi.*mask, parameters);
    KRHOCP        =     kRHOcp_Comp(     TEMP, PRESSURE, Z, RHO, CP, epsi.*mask, parameters);

    %% Extraction kientic
    Di            = Diffusion(RHO) .* 1e-14;
    shape         = Decay_Function_Coe(RHO);
    Dx            = axial_diffusion(TEMP, epsi, VELOCITY, RHO) .*1e-6;

    %% Saturation
    Csolid_percentage_left = 1 - (SOLID./C0solid);
    Csolid_percentage_left(find(~mask)) = 0;                                                % inserte zeros instead of NAN in pleces where there is no bed
    Sat_coe       =     Saturation_Concentration(Csolid_percentage_left, shape, Di);        % Inverse logistic is used to control saturation. Close to saturation point, the Sat_coe goes to zero.

    %% BC
    %Cf_0          =     if_else(F_u == 0, FLUID(1), 0);
    Cf_0          =     0;
    Cf_B          =     FLUID(nstages_index);
                                                                                            % If the sensitivity of P and F is consider, then set the input T as equal to the T inside of the extractor
                                                                                            % to avoid different small mismatch of T between the inlet and inside of the extractor
    T_0 = if_else( mode == 1, T_u, TEMP(1));
    
    T_B           =     TEMP(nstages_index);

    Z_0           =     Compressibility(T_0, PRESSURE,     parameters);
    %Z_B           =     Compressibility(T_B, PRESSURE,     parameters);

    rho_0         =     rhoPB_Comp(     T_0, P_u, Z_0, parameters);
    %rho_B         =     rhoPB_Comp(     T_B, PRESSURE, Z_B, parameters);

    u_0           =     Velocity(F_u, rho_0, parameters);
    %u_B           =     Velocity(F_u, rho_B, parameters);

    H_0           =     SpecificEnthalpy(T_0, PRESSURE, Z_0, rho_0, parameters );   

    enthalpy_rho_0 = if_else( mode == 3, ENTHALPY_RHO(1), rho_0 .* H_0 );                   % If the sensitivity of F is consider, then set the input h*rho as equal to the h*rho inside of the extractor
                                                                                            % to avoid different small mismatch betweenat the inlet and inside of the extractor
    
    %% Derivatives
    dz            = L/nstages_index;
    
    %dCfdz         = backward_diff_1_order(FLUID,Cf_0    , [],   dz);
    d2Cfdz2       = central_diff_2_order(FLUID, FLUID(1), Cf_B, dz);
    
    %dTdz          = backward_diff_1_order(TEMP,T_0, [] , dz);
    d2Tdz2        = central_diff_2_order(TEMP, T_0, T_B, dz);
        
    %dudz          = backward_diff_1_order(VELOCITY, u_0, [], dz);

    dHdz          = backward_diff_1_order(VELOCITY .* ENTHALPY_RHO, u_0 .* enthalpy_rho_0, [], dz);

    d_cons_CF_dz  = backward_diff_1_order(VELOCITY .* FLUID, u_0 .* Cf_0, [], dz);

    dPdt          = backward_diff_1_order(P_u, PRESSURE, [], dt)*1e2;
   
    re            = (Sat_coe ./ mi ./ lp2)  .* ( SOLID - FLUID .* rho_s ./ RHO ./ km );
    
    %% model
    xdot = [
    
    %--------------------------------------------------------------------
    % Concentration of extract in fluid phase | 0
   - 1            ./  ( 1 - epsi .* mask ) .* d_cons_CF_dz   + ...
    Dx            ./  ( 1 - epsi .* mask ) .* d2Cfdz2        + ...
    (epsi.*mask)  ./  ( 1 - epsi .* mask ) .* re;
    %zeros(nstages_index,1);

    %--------------------------------------------------------------------
    % Concentration of extract in solid phase | 1
    - mask                                 .* re;
    %zeros(nstages_index,1);
    
    %--------------------------------------------------------------------
    % enthalpy | 2
    - 1      ./ ( 1 - epsi .* mask )  .* dHdz +  dPdt - KRHOCP .* d2Tdz2;

    %--------------------------------------------------------------------
    % Pressure | 3
     dPdt;
    
    %--------------------------------------------------------------------
    % output equation
    %VELOCITY(nstages_index) * A * FLUID(nstages_index) * 1e3 ;   %kg/s - > g/s
    F_u ./ RHO(nstages_index) .* FLUID(nstages_index) * 1e3;

    Sdot;
    
    ];

end