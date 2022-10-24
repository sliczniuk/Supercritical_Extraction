function y = modelSFE_out2(x, nstages, rho, F_u, timeStep_in_sec)
%function y = modelSFE_out2(x, u, y_old, theta, timeStep_in_sec)

% theta = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];
%          1        2        3  4     5   6  7      8   9
    %{
    nstages       = theta{1};
    C0solid       = theta{2};                                            % Extractor initial concentration of extract
    V             = theta{3};     % Extractor volume (m3)
    epsi          = theta{4};     % Void bed fraction
    
    L             = theta{6};     % Length of the extractor (m)
    rho_s         = theta{7};     %
    
    U             = u';
    T_u           = U(1);
    P_u           = U(2);
    F_u           = U(3);   

    Z             = Compressibility(x(3*nstages),P_u,theta);
    RHO           = rhoPB_Comp(x(3*nstages),P_u,Z,theta);

    % total initial amount of solute in the solid phase
    Ms0 = C0solid * V * (1-epsi);
    
    Mf = x(nstages,:) * F_u / RHO * timeStep_in_sec;

    y = y_old + full( ( cumsum(Mf) / Ms0 )  * 100 );
    %}
    Mf = x(nstages,:) * F_u *1e3 / rho * timeStep_in_sec;

    y = cumsum(Mf) ;
    
end    