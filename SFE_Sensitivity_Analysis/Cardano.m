function [Z] = Cardaon(T, P, THETA)

    %import casadi.*
    %T = MX.sym('T',length(t));
    %P = MX.sym('P',length(p));
    %THETA = MX.sym('THETA',length(theta));

    %% Unpacking parameters
    % theta = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi, Tc, Pc, R, kappa];
    %          1        2        3  4     5   6  7      8   9   10  11  12 13

    Tc      = THETA{10};   % Critical temperature [K]
    Pc      = THETA{11};   % Critical pressure [bar]
    R       = THETA{12};   % Universal gas constant, [m3-bar/K-mol]
    kappa   = THETA{13};
    MW      = THETA{14};   % Molar mass [g/mol]
    Tr      = T ./ Tc;
    a       = 0.4572350 .* R.^2 .* Tc.^2 ./ Pc;
    b       = 0.0777961 .* R    .* Tc    ./ Pc;

    alpha = (1 + kappa .* (1 - sqrt(Tr))).^2;

    A = a .* alpha .* P ./ R.^2 ./ T.^2;
    B = b .* P ./ R ./ T;

    %%
    
    UC = - (1 - B)              ;
    SC =   (A - 2.*B - 3.*B.^2) ;
    TC = - ( A .* B - B.^2 - B.^3);

    PC = (3.*SC - UC.^2) ./ 3;
    QC = (2.*UC.^3 ./ 27) - UC.*SC./3 + TC;

    DC = (PC./3).^3 + (QC./2).^2;

    %% if D>0
    Z = ( sqrt(DC) - QC./2 ).^(1/3) - PC ./ (3 .* (sqrt(DC)-QC./2).^(1/3) ) - UC./3;

    %% if D<0
    %theta = sqrt(-PC^3 / 27);
    %phi   = acos(-QC/2/theta);

    %Z1    = 2*theta^(1/3) * cos(phi/3) - UC/3;
    %Z2    = 2*theta^(1/3) * cos(phi/3 + 2*pi/3) - UC/3;
    %Z3    = 2*theta^(1/3) * cos(phi/3 + 4*pi/3) - UC/3;

    %Z = [Z1, Z2, Z3];

end