clc, close all
clear all
addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

DATA = {'LUKE_T40_P200.xlsx', 'LUKE_T50_P200.xlsx', 'LUKE_T40_P300.xlsx', 'LUKE_T50_P300.xlsx'};

%% load parameters and set number of stages
Parameters_table     = readtable('Parameters.csv') ;
Parameters           = Parameters_table{:,3};
Parameters_sym       = MX(Parameters_table{:,3});

nstages              = 5;                                           

which_k              = [8,44];

%% Set time of the simulation
simulationTime       = 150;                                                 % Minutes
timeStep             = 1;                                                   % Minutes
SamplingTime         = 5;                                                   % Minutes
delayTime            = 0;                                                   % Minutes

timeStep_in_sec      = timeStep * 60;                                       % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                 = [0 Time_in_sec/60];                                  % Minutes

N_Sample             = find(Time == SamplingTime) - 1;
N_Delay              = delayTime / timeStep;
N_Time               = length(Time_in_sec);

%% symbolic variables
Nx = 3*nstages+1;
Nu = 3;
Nk = numel(which_k);
Np = Nu + numel(Parameters_sym);

% Create symbolic variables
x  = MX.sym('x', Nx);
u  = MX.sym('u', Nu);

% Store the results of the optimization
KOUT = nan(Nk,numel(DATA));

% Create the solver
OPT_solver              = casadi.Opti();
nlp_opts                = struct;
nlp_opts.ipopt.max_iter = 100;
ocp_opts                = {'nlp_opts', nlp_opts};
OPT_solver.solver(         'ipopt'   , nlp_opts)

% Descision variables
K                       = OPT_solver.variable(Nk);
k_lu                    = [ zeros(Nk,1) , 100*ones(Nk,1) ];

%% Set parameters
m_ref                    = 78;                                           % g of product obtained from a kg of biomass
                 
C0fluid                  = 0;                                            % Extractor initial concentration of extract
                                                                         % Fluid phase kg / m^3
                 
V                        = 0.00165;                                      % Volume of the extractor  [m3] 
L                        = 0.095;                                        % Length of the extractor [m]
epsi                     = 2/3;                                          % Porosity [-] 
                 
C0solid                  = m_ref *1e-3 / (V*(1-epsi));                   % Solid phase kg / m^3
%{             
dp                   = 0.00010;                                      % Diameter of the particle [m] - Vargas
rho_s                = 1300.0;                                       % Densisty of the solid phase [kg / m^3] - FC
km                   = 0.75;                                         % Partition coefficient (?)
             
mi                   = 1/2;                                          % Geometric shape coefficient (?) 

% Assign new values of parameters to the Parameters vector 
%                       nstages, C0solid, V, epsi, dp, L, rho_s, km, mi                
Parameters(1:9)      = [nstages, C0solid, V, epsi, dp, L, rho_s, K(1), mi];
%}
Parameters_sym(which_k)  = K;
p                        = [u; Parameters_sym];

k0 = [0.01; 1];

%{
nstages = 2;                                            %

m_ref   = 78;                                           % g of product obtained from a kg of biomass

C0fluid = 0;                                            % Extractor initial concentration of extract
% Fluid phase kg / m^3

V       = 0.00165;                                      % Volume of the extractor  [m3] - Vargas
L       = 0.095;                                        % Length of the extractor [m] - Vargas
epsi    = 2/3;                                          % Porosity [-] - Vargas

C0solid = m_ref *1e-3 / (V*(1-epsi));                   % Solid phase kg / m^3

dp      = 0.00010;                                      % Diameter of the particle [m] - Vargas
rho_s   = 1300.0;                                       % Densisty of the solid phase [kg / m^3] - FC
km      = 0.29;                                         % Partition coefficient (?)

mi      = 1/2;                                          % Geometric shape coefficient (?) - Vargas

Ti      = [313.15, 323.15, 333.15, 343.15];
Di      = [3.03E-11, 3.66E-11, 1.97E-10, 2.73E-10];
T_kp    = [313.15, 323.15, 333.15, 343.15];
dkp     = [8.13E-01 6.94E-01 1.04E-01 6.67E-02];

Tc      = 304.1;                                        % Critical temperature [K]
Pc      = 73.8;                                         % Critical pressure [bar]
omega   = 0.228;                                        % Acentric factor [-]
R       = 8.314e-5;                                     % Universal gas constant, [m3-bar/K-mol]

kappa   = 0.37464 + 1.54226 * omega - 0.26992 * omega^2;

MW      = 44e-3;                                        % Molar weight of CO2 [Kg/mol]

[EA_Di, betah_Di] = Arrhenius_Di(Ti,Di,R);
[EA_km, betah_km] = Arrhenius_km(T_kp,dkp,R,rho_s);

CP_0 =   4.18686;
CP_A =  19.80;        %4.5980;
CP_B =   7.344E-2;    %0.0125;
CP_C = - 5.602E-5;    %2.86E-06;
CP_D =   1.715E-8;    %-2.70E-09;

cpSolid = 1.5E3;      % J / K / Kg

% Axials diffusions correlations
a_axial = 0.1;
b_axial = 0.011;
c_axial = 0.48;

% Heats Conductivitys
A1_cond = -105.161;
A2_cond =  0.9007;
A3_cond =  0.0007;
A4_cond =  3.50E-15;
A5_cond =  3.76E-10;
A6_cond =  0.7500;
A7_cond =  0.0017;

% Viscositys
A1_visc = -1.146067E-1;
A2_visc =  6.978380E-7;
A3_visc =  3.976765E-10;
A4_visc =  6.336120E-2;
A5_visc = -1.166119E-2;
A6_visc =  7.142596E-4;
A7_visc =  6.519333E-6;
A8_visc = -3.567559E-1;
A9_visc =  3.180473E-2;

%                 1        2    3    4    5  6     7   8   9   10  11  12   13   14   15      16       17    18    19    20    21    22      23        24      25        26       27      28       29        30        31       32       33       34      35       36        37      38       39       40      41        42       43   
parameters = {nstages, C0solid, V, epsi, dp, L, rho_s, km, mi, Tc, Pc, R, kappa, MW, EA_Di, betah_Di, CP_0, CP_A, CP_B, CP_C, CP_D, EA_km, betah_km, cpSolid, a_axial, b_axial, c_axial, A1_cond, A2_cond, A3_cond, A4_cond, A5_cond, A6_cond, A7_cond, A1_visc, A2_visc, A3_visc, A4_visc, A5_visc, A6_visc, A7_visc, A8_visc, A9_visc};
name       = {'nstages', 'C0solid', 'V', 'epsi', 'dp', 'L', 'rho_s', 'km', 'mi', 'Tc', 'Pc', 'R', 'kappa', 'MW', 'EA_Di', 'betah_Di', 'CP_0', 'CP_A', 'CP_B', 'CP_C', 'CP_D', 'EA_km', 'betah_km', 'cpSolid', 'a_axial', 'b_axial', 'c_axial', 'A1_cond', 'A2_cond', 'A3_cond', 'A4_cond', 'A5_cond', 'A6_cond', 'A7_cond', 'A1_visc', 'A2_visc', 'A3_visc', 'A4_visc', 'A5_visc', 'A6_visc', 'A7_visc', 'A8_visc', 'A9_visc'};
%}

%% Set Integrator
tic
%f_r = @(x, u, k) modelSFE_Regression(x, u, k, parameters);
f = @(x, p) modelSFE(x, p, nstages);
% Integrator
F = buildIntegrator(f, [Nx,Np] , timeStep_in_sec);
toc

%% Loop over datasets

%for i=1%:numel(DATA)
i = 1;

    LabResults = xlsread(DATA{i});

    T0homog   = LabResults(1,1)+273.15;
    feedPress = LabResults(1,2);
    rho       = LabResults(1,3);
    data      = LabResults(:,5)';

    V_Flow = 0.39;

    %%
    feedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars

    feedFlow  = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec

    %% Inital values
    x0 = [C0fluid*ones(nstages,1);
          C0solid*ones(nstages,1);
          T0homog*ones(nstages,1);
          0];

    uu = [feedTemp', feedPress', feedFlow'];

    X = MX(Nx,N_Time+1);
    X(:,1) = x0;

    %%
    tic
    for j=1:N_Time
        X(:,j+1)=F(X(:,j), [uu(j,:)'; Parameters_sym] );
    end
    toc

    %%
    Yield_estimate = X(3*nstages+1,:);
    Yield_estimate = [zeros(1,N_Delay) Yield_estimate(1:end-N_Delay)];
    Yield_estimate = Yield_estimate(1:N_Sample:end);

    J = (data-Yield_estimate ) * diag(1:1:1) * (data-Yield_estimate )'; 

    %%
    for nk=1:Nk
          OPT_solver.subject_to( k_lu(nk,1) <= K(nk,:) <= k_lu(nk,2) )
    end

    %%
    OPT_solver.minimize(J);

    OPT_solver.set_initial(K, k0);

    %%
    % tic
    try
        sol = OPT_solver.solve();
        kout = sol.value(K)
    catch
        kout = OPT_solver.debug.value(K);
    end
    toc

    Parameters(which_k) = kout;
    Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

    [xx_out] = simulateSystem(F, [], x0, Parameters_opt  );

    %%
    hold on
    plot(Time,xx_out(end,:))
    plot(Time(1:N_Sample:end),data,'o')
    hold off
    
    %%
%{
    [yy_out, tt_out, xx_out] = simulateSystem_ParameterEstimation(F, g, x0, uu, kout );

    %%
    subplot(2,2,i)
    hold on
    plot(Time,xx_out(end,:))
    plot(Time(1:N_Sample:end),data,'o');
    hold off
    ylabel('Mass of product [g]')
    xlabel('Time [min]')
    caption = sprintf('T = %g, P = %g, \\rho = %g', T0homog, feedPress(1), rho );
    title(caption);

    KOUT(:,i) = kout;

end
%save fit4.mat
%}
%%
%{
for i=1:numel(DATA)
    LabResults = xlsread(DATA{i});

    T0homog   = LabResults(1,1)+273.15;
    feedPress = LabResults(1,2);
    rho       = LabResults(1,3);
    data      = LabResults(:,5)';

    feedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars

    feedFlow  = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec


    % Inital values
    x0 = [C0fluid*ones(nstages,1);
        C0solid*ones(nstages,1);
        T0homog*ones(nstages,1);
        0];

    uu = [feedTemp', feedPress', feedFlow'];

    [yy_out, tt_out, xx_out] = simulateSystem_ParameterEstimation(F_r, g, x0, uu, KOUT(:,i) );

    subplot(2,2,i)
    hold on
    plot(xx_out(end,:))
    plot(Time(1:N_Sample:end),data,'o');
    hold off
    ylabel('Mass of product [g]')
    xlabel('Time [min]')
    caption = sprintf('T = %g, P = %g, \\rho = %g', T0homog, feedPress(1), rho );
    title(caption);
end
%}