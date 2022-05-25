clc, clear all, close all, clearvars -global
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

%% Parameters
nstages = 3;                                                   %

dataset = 'Vargas';

switch dataset
    case 'Reverchon'
        C0solid = 9;                                            % Extractor initial concentration of extract - Reverchon
        V       = 4e-4;                                         % Volume of the extractor (400 ml)
        epsi    = 0.8;                                          % Porosity [-]
        dp      = 0.75e-3;                                      % Diameter of the particle [m]
        L       = 0.165;                                        % Length of the extractor [m]
        rho_s   = 413.25;                                       % Densisty of the solid phase [kg / m^3]
        km      = 0.29;                                         % Partition coefficient (?)
        mi      = 3/5;                                          % Geometric shape coefficient (?)
        Ti      = (273.15 + [25, 40, 70]);                      % Data temperatures of some biomass
        Di      = [ 3.04, 2.88, 4.64] * 10^-12;                 % Data diffusivities of some biomass - King
        %Di      = [ 6, 6, 6] * 10^-13;                         % Data diffusivities of some biomass - Reverchon
        T_kp    = [313.15, 323.15, 333.15, 343.15];
        dkp     = (rho_s .*km ./[485.15 285 235.4 212.4]);
        feedFlow  = 0.53/3600;                                  % Kg / sec
        disp('Reverchon dataset')
    case 'Vargas'
        C0solid = 27;                                           % Vargas
        V       = 4.173e-6;                                     % Volume of the extractor  [m3] - Vargas
        epsi    = 0.669;                                        % Porosity [-] - Vargas
        dp      = 5e-4;                                         % Diameter of the particle [m] - Vargas
        L       = 0.045;                                        % Length of the extractor [m] - Vargas
        rho_s   = 1087.2;                                       % Densisty of the solid phase [kg / m^3] - Vargas
        km      = 0.29;                                         % Partition coefficient (?)
        mi      = 1/5;                                          % Geometric shape coefficient (?) - Vargas
        Ti      = [313.15, 323.15, 333.15, 343.15];
        Di      = [ 3.03E-11, 3.66E-11, 1.97E-10, 2.73E-10];
        T_kp    = [313.15, 323.15, 333.15, 343.15];
        dkp     = [8.13E-01 6.94E-01 1.04E-01 6.67E-02];
        feedFlow  = 0.035/3600;                                 % Kg / sec - vargas
        disp('Vargas dataset')
end

T0homog = 50 + 273.15;                                          % Extractor initial temperature (pseudo-homogeneous)
C0fluid = 0;                                                    % Extractor initial concentration of extract
                                                                % Fluid phase kg / m^3
Tc      = 304.1;                                                % Critical temperature [K]
Pc      = 73.8;                                                 % Critical pressure [bar]
omega   = 0.228;                                                % Acentric factor [-]
R       = 8.314e-5;                                             % Universal gas constant, [m3-bar/K-mol]

kappa   = 0.37464 + 1.54226 * omega - 0.26992 * omega^2;

MW      = 44e-3;                                                % Molar weight of CO2 [Kg/mol]

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

parameters = {nstages, C0solid, V, epsi, dp, L, rho_s, km, mi, Tc, Pc, R, kappa, MW, EA_Di, betah_Di, CP_0, CP_A, CP_B, CP_C, CP_D, EA_km, betah_km, cpSolid, a_axial, b_axial, c_axial, A1_cond, A2_cond, A3_cond, A4_cond, A5_cond, A6_cond, A7_cond, A1_visc, A2_visc, A3_visc, A4_visc, A5_visc, A6_visc, A7_visc, A8_visc, A9_visc };

%% Time
simulationTime       =  30;                                      % Minutes 
timeStep             =  10;                                      % Minutes  
timeStep_in_sec      = timeStep*60;                               % Seconds
Time_in_sec          = (0:timeStep:simulationTime)*60;            % Seconds
Time                 = Time_in_sec/60;                            % Minutes

%% 
feedTemp  = T0homog  * ones(1,length(Time_in_sec));  % Kelvin
feedPress = 90       * ones(1,length(Time_in_sec));  % Bars
feedFlow  = feedFlow * ones(1,length(Time_in_sec));  % 

%% Casadi variables and Model 

% Create symbolic variables
u     = MX.sym('u', 3);
x     = MX.sym('x', 3*nstages);
k     = MX.sym('k', 1);

%Variables
Nx = 3*nstages;
Nu = 3;
Nk = 1;
Ny = 1;

%% Model
% Model Equations
f = @(x, u, k) modelSFE_SS(x, u, k, parameters);
g = @(x) modelSFE_out(x, parameters);

% Integrator
F = buildIntegrator_ParameterEstimation(f, [Nx,Nu,Nk] , timeStep_in_sec);   

%% Inital values
x0 = [C0fluid*ones(nstages,1);
      C0solid*ones(nstages,1);
      T0homog*ones(nstages,1)];

uu = [feedTemp', feedPress', feedFlow'];

k0 = [EA_Di];

%% Simulate system with initial parameter
[yy, tt, xx] = simulateSystem_ParameterEstimation(F, g, x0, uu, k0);
data = yy*0.95;
fprintf('Initial simulation is finished\n')

%%
OCP = struct('Nx', Nx, 'Nu', Nu, 'Nk', Nk, 'Ny', Ny, 'N', length(uu), ...
    'x_lu', [0*ones(Nx,1) inf*ones(Nx,1)],  ...
    'k_lu', [ ],  ...
    'x_eq', [], ...
    'k_eq', [], ...
    'F'   , F, ...
    'L'   , [  ], ...       % 
    'Lf'  , @(x,yd) sum( (g(x)-yd) .* (g(x)-yd)  )   ); %'Lf'  , @(x,yd) sum( (g(x)-yd).^2 )   );   % 0.5*( sum(g(x)-yd) )'*( sum(g(x)-yd) ) 

%% Solve LSE
[~, kout] = singleShooting_ParameterEstimation(OCP, x0, uu', data, k0)

%% Simulate system with estiamted parameter
[yy_out, tt_out, xx_out] = simulateSystem_ParameterEstimation(F, g, x0, uu, kout);

%% Plot
hold on 
plot( data, 'o')
plot( yy )
plot( yy_out,'.-')  
hold off