startup

import casadi.*

Parameters_table        = readtable('Parameters.csv') ;        % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});

DATA_set                = {'LUKE_T40_P200', 'LUKE_T50_P200', 'LUKE_T40_P300_org', 'LUKE_T50_P300_org'};

V_Flow                  = 0.4;                                              % Volumetric flow rate l/mi
r                       = Parameters{3};                                    % Radius of the extractor  [m]
epsi                    = Parameters{4};                                    % Fullness [-]
L                       = Parameters{6};                                    % Total length of the extractor [m]

RE = []; RHO = []; MU = []; V = []; D_L = [];

for ii=1:numel(DATA_set)
    % Dataset
    DATA                = DATA_set{ii};
    LabResults          = xlsread([DATA,'.xlsx']);

    T                   = LabResults(1,1)+273.15;
    P                   = LabResults(1,2);

    Z                   = Compressibility( T, P,         Parameters );
    rho                 = rhoPB_Comp(      T, P, Z,      Parameters );
    RHO = [RHO,rho];

    F                   = V_Flow * rho * 1e-3 / 60 ;

    v                   = Velocity(F,rho,Parameters);
    V  = [V,v];

    mu                  = Viscosity(T,rho);
    MU = [MU,mu];

    dp = 0.15;

    Re                  = dp * v * rho * epsi / mu;
    RE = [RE,Re];

    d_L                 = 2*0.001*v/(0.2+0.011*Re^0.48);
    D_L = [D_L, d_L];

end

DATA_estimated = readtable('estimation.csv');
D_L_estimated  = DATA_estimated{3,2:end};

[trendLine, SLine]   = polyfit(RE, D_L_estimated, 1 );
p_est                = polyval(trendLine, RE);
p_y                  = polyval(trendLine, linspace(min(RE), max(RE) ) );

subplot(3,3,3)
scatter(RE, D_L_estimated, 'filled')
hold on
plot(linspace(min(RE), max(RE) ),p_y,'LineWidth',2);
hold off

pp=sprintf('$D^M_e=(%6.4g) R_e + (%6.4g)$',trendLine');

oo = sprintf('$R^2 = %g$\n', round(power(corr2(D_L_estimated, p_est),2),2) ) ;

title( [oo,pp], 'FontSize', 8)

xlabel('$Re~[-]$')
ylabel('$D^M_e [m^2/s] \times 1e-6$')

set(gcf,'PaperOrientation','landscape'); print(figure(1),['Trend_Lines_D_e.pdf'],'-dpdf','-bestfit');
close all;