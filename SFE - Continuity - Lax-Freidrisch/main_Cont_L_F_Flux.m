clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

Parameters_table        = readtable('Parameters.csv') ;        % Fulle table with prameters

%% Set time of the simulation
simulationTime          = 150;                                        % Minutes

timeStep                = simulationTime/500;                               % Minutes

timeStep_in_sec         = timeStep * 60;                                       % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                    = [0 Time_in_sec/60];                                  % Minutes

N_Time                  = length(Time);

%% Specify parameters to estimate
nstages                 = 300;

before  = 0.135;        nstagesbefore   = 1:floor(before*nstages);
bed     = 0.165;        nstagesbed      = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
                        nstagesafter    = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

which_k                 = [8, 44];

%% Set parameters
mSOL_s                   = 100;                                         % g of product in biomass
mSOL_f                   = 0;                                           % g of biomass in fluid

%C0fluid                 = 1;                                           % Extractor initial concentration of extract - Fluid phase kg / m^3

V                       = 0.01;                                         %
r                       = 0.075;                                        % Radius of the extractor  [m3]
L                       = V / pi / r^2;                                 % Total length of the extractor [m]
L_stages                = linspace(0, L, nstages);
A                       = pi*r^2;                                       % Extractor cross-section
epsi                    = 0.1 ;                                         % Fullness [-]

%--------------------------------------------------------------------
V_slice                 = (L/nstages) * pi * r^2;
V_before                = V_slice * numel(nstagesbefore);
V_after                 = V_slice * numel(nstagesafter);
V_bed                   = V_slice * numel(nstagesbed);                  % Volume of the fixed bed [m3]

V_before_solid = repmat(V_before * 0 / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_solid    = repmat(V_bed * epsi / numel(nstagesbed),    numel(nstagesbed),1);
V_after_solid  = repmat(V_after * 0  / numel(nstagesbed),    numel(nstagesafter),1);

V_solid = [V_before_solid; V_bed_solid; V_after_solid];

V_before_fluid  = repmat(V_before * 1       / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_fluid     = repmat(V_bed * (1 - epsi) / numel(nstagesbed),    numel(nstagesbed),1);
V_after_fluid   = repmat(V_after * 1        / numel(nstagesafter),  numel(nstagesafter),1);

V_fluid = [V_before_fluid; V_bed_fluid; V_after_fluid];

%--------------------------------------------------------------------
dp                      = 0.00010;                                      % Diameter of the particle [m] - Vargas
rho_s                   = 1250.0;                                       % Densisty of the solid phase [kg / m^3] - FC
km                      = 0.29;                                         % Partition coefficient (?)

mi                      = 1/2;                                          % Geometric shape coefficient (?)

C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;            % Solid phase kg / m^3

C0fluid                 = mSOL_f * 1e-3 / (V_before + V_bed * (1-epsi) + V_after);

m0fluid(nstagesbefore) = C0fluid * V_before / numel(nstagesbefore);
m0fluid(nstagesbed)    = C0fluid * (V_bed * (1 - epsi)) / numel(nstagesbed);
m0fluid(nstagesafter)  = C0fluid * V_after / numel(nstagesafter);

%%
Nx                      = 5*nstages+1;
Nu                      = 3 + numel( Parameters_table{:,3} );
Nk                      = numel(which_k);

%% symbolic variables
% Create symbolic variables
x                       = MX.sym('x', Nx);
u                       = MX.sym('u', Nu);

%% Velocity profile: linear change
u_in  = 1.00;
u_out = 0.8;

%% Set Integrator
f_cons                  = @(x, u) modelSFE_Cont_L_F_Flux(x, u, bed_mask, u_in, u_out, timeStep_in_sec, 'conservative');

% Integrator
F_cons                  = buildIntegrator(f_cons, [Nx,Nu] , timeStep_in_sec);

%
f_flux                  = @(x, u) modelSFE_Cont_L_F_Flux(x, u, bed_mask, u_in, u_out, timeStep_in_sec, 'flux');

% Integrator
F_flux                  = buildIntegrator(f_flux, [Nx,Nu] , timeStep_in_sec);

%%
V_Flow     = 0.39;
T0homog    = 40+273.15;
feedPress  = 300 ;
k0         = [1, 0.1];

Parameters          = Parameters_table{:,3};
Parameters(1:9)     = [nstages, C0solid, r, epsi, dp, L, rho_s, km, mi];

Parameters(which_k) = k0;

%%
rho        = rhoPB_Comp(T0homog, feedPress.* linspace( u_in, u_out, nstages)', Compressibility(T0homog,feedPress.* linspace( u_in, u_out, nstages)', num2cell(Parameters) ), num2cell(Parameters) );

% Set operating conditions
feedTemp   = T0homog   * ones(1,N_Time) + 0;  % Kelvin
%feedTemp( round(numel(Time)/4) : round(numel(Time)/2) )   = T0homog + 0;
%feedTemp( round(numel(Time)/2)  : end )   = T0homog +20 ;

feedPress  = feedPress * ones(1,N_Time) + 0 ;  % Bars
%feedPress(round(numel(Time)/2):round(3*numel(Time)/4))   = feedPress(1) - 10;

%feedFlow   = V_Flow * rho * 1e-3 / 60 * ones(1,N_Time);  % l/min -> kg/min -> Kg / sec
feedFlow    = V_Flow * 1e-3 / 60 * ones(1,N_Time);  % l/min -> kg/min -> m3 / sec
%feedFlow(round(numel(Time)/20):round(numel(Time)/10))   = 2*feedFlow(1) ;

uu         = [feedTemp', feedPress', feedFlow'];

Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

%% Initial conditions
x0_cons     = [
            C0fluid * ones(nstages,1) .* ( 1 - epsi .* bed_mask );
            C0solid * bed_mask;
            T0homog * ones(nstages,1);
            rho     .* ( 1 - epsi .* bed_mask );        
            feedFlow(1) / A * ones(nstages,1)  ;
            0;
            ];

%% Simulate system
[xx_0_cons] = simulateSystem(F_cons, [], x0_cons, Parameters_opt  );
[xx_0_flux] = simulateSystem(F_flux, [], x0_cons, Parameters_opt  );

%%
figure(1)
PlotResults(xx_0_cons, Time, nstagesbed, Parameters, uu, bed_mask, epsi, 'conservative')
sgtitle ('conservative')

figure(2)
PlotResults(xx_0_flux, Time, nstagesbed, Parameters, uu, bed_mask, epsi, 'conservative')
sgtitle ('flux')

figure(3)
subplot(1,2,1)
plotyy(Time, xx_0_cons(Nx,:), Time, 1e3 * (sum(xx_0_cons(0*nstages+1:1*nstages,:) .* V_fluid ./ ( 1 - epsi .* bed_mask ) ) + sum(xx_0_cons(1*nstages+1:2*nstages,:) .* V_solid)) + xx_0_cons(Nx,:))
title ('conservative')

subplot(1,2,2)
plotyy(Time, xx_0_flux(Nx,:), Time, 1e3 * (sum(xx_0_flux(0*nstages+1:1*nstages,:) .* V_fluid ./ ( 1 - epsi .* bed_mask ) ) + sum(xx_0_flux(1*nstages+1:2*nstages,:) .* V_solid)) + xx_0_flux(Nx,:))
title ('flux')


%%
%xx_0(end,end)
xx_0_cons(end,end)