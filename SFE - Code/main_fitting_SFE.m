clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

DATA = {'Data_40_200.mat', 'Data_50_200.mat', 'Data_40_300.mat', 'Data_50_300.mat'};

%% Parameters
nstages = 5;                                                   %

C0fluid = 0;                                            % Extractor initial concentration of extract
% Fluid phase kg / m^3
C0solid = 32;                                           % Solid phase kg / m^3
V       = 0.01/4;                                       % Volume of the extractor  [m3] - Vargas
epsi    = 0.70;                                         % Porosity [-] - Vargas
dp      = 0.001;                                        % Diameter of the particle [m] - Vargas
L       = 0.70;                                         % Length of the extractor [m] - Vargas
rho_s   = 1087.2;                                       % Densisty of the solid phase [kg / m^3] - Vargas
km      = 0.29;                                         % Partition coefficient (?)
mi      = 1/5;                                          % Geometric shape coefficient (?) - Vargas
Ti      = [313.15, 323.15, 333.15, 343.15];
Di      = [ 3.03E-11, 3.66E-11, 1.97E-10, 2.73E-10];
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

% axial diffusion correlation
a_axial = 0.1;
b_axial = 0.011;
c_axial = 0.48;

% Heat Conductivity
A1_cond = -105.161;
A2_cond =  0.9007;
A3_cond =  0.0007;
A4_cond =  3.50E-15;
A5_cond =  3.76E-10;
A6_cond =  0.7500;
A7_cond =  0.0017;

% Viscosity
A1_visc = -1.146067E-1;
A2_visc =  6.978380E-7;
A3_visc =  3.976765E-10;
A4_visc =  6.336120E-2;
A5_visc = -1.166119E-2;
A6_visc =  7.142596E-4;
A7_visc =  6.519333E-6;
A8_visc = -3.567559E-1;
A9_visc =  3.180473E-2;

sigma    = 1;
sigma_km = 0.3;
mu_km    = 0.1;
Di_slack = 1;

%                 1        2    3    4    5  6     7   8   9   10  11  12   13   14   15      16       17    18    19    20    21    22      23        24      25        26       27      28       29        30        31       32       33       34      35       36        37      38       39       40      41        42       43      44        45   
parameters = {nstages, C0solid, V, epsi, dp, L, rho_s, km, mi, Tc, Pc, R, kappa, MW, EA_Di, betah_Di, CP_0, CP_A, CP_B, CP_C, CP_D, EA_km, betah_km, cpSolid, a_axial, b_axial, c_axial, A1_cond, A2_cond, A3_cond, A4_cond, A5_cond, A6_cond, A7_cond, A1_visc, A2_visc, A3_visc, A4_visc, A5_visc, A6_visc, A7_visc, A8_visc, A9_visc, Di_slack, sigma };
which_parameter= {8,44,45};
theta = parameters;
close all

%% Time
simulationTime       = 150;                                                 % Minutes
delayTime            = 2;                                                   % Minutes
timeStep             = 1;                                                   % Minutes
timeStep_in_sec      = timeStep *60;                                        % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                 = [0 Time_in_sec/60];                                  % Minutes
SamplingTime         = 5;                                                   % Minutes
N_Sample             = find(abs(Time-SamplingTime) < timeStep/10)-1;
N_Delay              = delayTime / timeStep;

%% Casadi variables and Model

% Create symbolic variables
u  = MX.sym('u', 3);
x  = MX.sym('x', 3*nstages+1);
k  = MX.sym('k', length(which_parameter));

%Variables
Nx = 3*nstages+1;
Nu = 3;
Nk = length(which_parameter);
Ny = 1;

%% Model
% Model Equations
f = @(x, u, k) modelSFE_SS(x, u, k, parameters, which_parameter);
%g = @(x) modelSFE_out(x, parameters);
g = @(x, u, y_old) modelSFE_out2(x, u, y_old, parameters, timeStep_in_sec);

% Integrator
F = buildIntegrator_ParameterEstimation(f, [Nx,Nu,Nk] , timeStep_in_sec);

%% dummy parameters
RHO = [];
k0 = [0.1;1;1];
KOUT = k0;

Operating_Conditions_experiments = nan(2,numel(DATA));

 %%
for i = 1:numel(DATA)
    load(DATA{i});
    Operating_Conditions_experiments(:,i) = [T0homog;feedPress(1)];

    %T0homog   = 50 + 273.15;                             % Extractor initial temperature (pseudo-homogeneous)
    feedTemp  = T0homog      * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress(1) * ones(1,length(Time_in_sec));  % Bars

    Z = Compressibility(T0homog, feedPress(1), parameters);
    rho = full(rhoPB_Comp(T0homog, feedPress(1), Z, parameters));
    RHO = [RHO, rho];

    feedFlow  = 0.4 * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec

    %% Inital values
    x0 = [C0fluid*ones(nstages,1);
        C0solid*ones(nstages,1);
        T0homog*ones(nstages,1);
        0];

    uu = [feedTemp', feedPress', feedFlow'];

    data2 = data;

    %%
    OCP = struct('Nx', Nx, 'Nu', Nu, 'Nk', Nk, 'Ny', Ny, 'N', length(Time_in_sec), ...
        'x_lu', [],  ...                % 0*ones(Nx,1) inf*ones(Nx,1)
        'k_lu', [0*ones(Nk,1) inf*ones(Nk,1) ],  ...                % 0*ones(Nk,1) [1.5;inf;inf]; 0*ones(Nk,1) inf*ones(Nk,1)
        'x_eq', [], ...
        'k_eq', [], ...
        'N_Sample', N_Sample,...
        'N_Delay', N_Delay,...
        'F'   , F, ...
        'L'   , [], ...       %
        'LS'  , @(x) sum( (x-data2) .^2 ),...
        'MSE', @(x,sigma_MSE) -(-length(data2)/2*log(2*pi*sigma_MSE^2) - 1/(2*sigma_MSE^2) * sum( (x-data2) .^2 )), ...                 %'Lf'  , @(x,yd) sum( (g(x)-yd).^2 )   );   % 0.5*( sum(g(x)-yd) )'*( sum(g(x)-yd) ) %'Lf'  , @(x,yd) sum( (g(x)-yd).^2 )   );   % 0.5*( sum(g(x)-yd) )'*( sum(g(x)-yd) )
        'MAP', @(x,km,sigma_MSE) -( -length(data2)/2*log(2*pi*sigma_MSE^2) - 1/(2*sigma_MSE^2) * sum( (x-data2) .^2 ) - 1/2*log(2*pi*sigma_km^2) - 1/(2*sigma_km^2) * sum( (km-mu_km) .^2 )  ) ...
        );


    %% Solve optimization problem
    [~, kout] = singleShooting_ParameterEstimation(OCP, x0, uu',  KOUT(:,1),  parameters);
    KOUT = [KOUT, kout];
    
end
KOUT = KOUT(:,2:end);

%% plot the parameter fitting
plot_fitting(DATA, Time_in_sec, RHO, F, g, x0, KOUT, which_parameter, N_Delay, N_Sample, theta, 'print_of');

%% fit the parameters from different operating conditions
[Regression_SFE] = linear_regression_SFE(KOUT,RHO,Operating_Conditions_experiments, 'print_of');

%% Model with regression
%regression_fitting_results(DATA, Time_in_sec, nstages, RHO, x, u, k, F, g, x0, KOUT, which_parameter, N_Delay, N_Sample, theta, Regression_SFE, 'print_off')