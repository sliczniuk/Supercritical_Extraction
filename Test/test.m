clc, clear all, close all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

%% States

% Create symbolic variables
x     = MX.sym('x', 3);                                   % states
u     = MX.sym('u', 1);                                   % controls
k     = MX.sym('k', 2);                                   % parameters to estimate

% Number of Variables
Nx = 3;
Nu = 1;
Nk = 2;
Ny = 1;

%% Time
simulationTime       = 20/60;                                     % Minutes 
timeStep             = 0.01 ;                                     % Minutes  
timeStep_in_sec      = timeStep*60;                               % Seconds
Time_in_sec          = (0:timeStep:simulationTime)*60;            % Seconds
Time                 = Time_in_sec/60;                            % Minutes

%% Parameters
k0 = [0.1, 0.5];
theta = 1;

%% Controls
uu = zeros(length(Time_in_sec),1);

%% Model
% Model Equations
f = @(x, u, k) test_model(x, u, k, theta);
g = @(x) model_out(x);

% Integrator
F = buildIntegrator_ParameterEstimation(f, [Nx,Nu,Nk] , timeStep_in_sec);       

% Initial conditions
x0 = [1; 0; 0.1];
k_initial = [0.08 0.55];

% Simulate system with initial parameter
[yy, tt, xx] = simulateSystem_ParameterEstimation(F, g, x0, uu, k0);
data = yy;

 %% Optimization

OCP = struct('Nx', Nx, 'Nu', Nu, 'Nk', Nk, 'Ny', Ny, 'N', length(uu), ...
    'x_lu', [0*ones(Nx,1) inf*ones(Nx,1)],  ...
    'k_lu', [ ],  ...
    'x_eq', [], ...
    'k_eq', [], ...
    'F'   , F, ...
    'L'   , [  ], ...       % 
    'Lf'  , @(x,yd) 0.5*( sum(g(x)-yd) )'*( sum(g(x)-yd) )    );

%% Solve LSE
[~, kout] = singleShooting_ParameterEstimation(OCP, x0, uu', data, k_initial)

%% Simulate system with estiamted parameter
[yy_out, tt_out, xx_out] = simulateSystem_ParameterEstimation(F, g, x0, uu, kout);

%% Plot
hold on 
plot( yy, 'o')
plot( xx')
plot( xx_out','.-')  
hold off



