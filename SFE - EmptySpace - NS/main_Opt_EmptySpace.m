clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

%DATA = {'LUKE_T40_P200.xlsx', 'LUKE_T50_P200.xlsx', 'LUKE_T40_P300.xlsx', 'LUKE_T50_P300.xlsx'};
%DATA = {'LUKE_T40_P300.xlsx'};

Parameters_table        = readtable('Parameters.csv') ;        % Fulle table with prameters

%% Set time of the simulation
simulationTime          = 250;                                                 % Minutes
SamplingTime            = 5;                                                   % Minutes

timeStep                = 1/4;                                                 % Minutes

timeStep_in_sec         = timeStep * 60;                                       % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                    = [0 Time_in_sec/60];                                  % Minutes
%--------------------------------------------------------------------
SAMPLE                  = [0:SamplingTime:150];
SAMPLE(1)               = 0;
%{
N_Sample = [];
for i = 1:numel(SAMPLE)
    N_Sample = [N_Sample ; find(Time == SAMPLE(i))];
end
if numel(N_Sample) ~= numel(SAMPLE)
    keyboard
end
%}

N_Time                  = length(Time_in_sec);

%% Specify parameters to estimate
nstages                 = 100;

before  = 0.14;         nstagesbefore   = 1:floor(before*nstages);
bed     = 0.16;         nstagesbed      = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
nstagesafter    = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

which_k                 = [8, 44];

%% Set parameters
mSOL_s                   = 100;                                           % g of product in biomass
mSOL_f                   = 0;                                           % g of biomass in fluid

%C0fluid                 = 1;                                           % Extractor initial concentration of extract - Fluid phase kg / m^3

V                       = 0.01;                                         %
r                       = 0.075;                                        % Radius of the extractor  [m3]
L                       = V / pi / r^2;                                 % Total length of the extractor [m]
epsi                    = 1/3;                                          % Fullness [-]

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

%% Set Integrator
f                       = @(x, u) modelSFE(x, u, bed_mask);

% Integrator
F                       = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

%%
V_Flow     = 0.39;
T0homog    = 40+273.15;
feedPress  = 300 ;
k0         = [0.1, 0.5];

%%
rho        = rhoPB_Comp(T0homog, feedPress, Compressibility(T0homog,feedPress,table2cell(Parameters_table(:,3))), table2cell(Parameters_table(:,3)));

% Set operating conditions
feedTemp   = T0homog   * ones(1,length(Time_in_sec))  ;  % Kelvin
feedTemp(round(numel(Time)/4):round(numel(Time)/2))   = feedTemp(1) + 10;

feedPress  = feedPress * ones(1,length(Time_in_sec));  % Bars
%feedPress(round(numel(Time)/3):round(2*numel(Time)/3))   = feedPress(1) + 50;

feedFlow   = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec
%feedFlow(round(numel(Time)/3):round(2*numel(Time)/3))   = 2*feedFlow(1) ;

uu         = [feedTemp', feedPress', feedFlow'];

% Initial conditions
x0         = [C0fluid * ones(nstages,1);
            C0solid * bed_mask;
            T0homog*ones(nstages,1);
            rho*ones(nstages,1);                                                                                           % rho*ones(nstages,1);
            Velocity(feedFlow(1), rho(1), table2cell(Parameters_table(:,3)) ) * ones(nstages,1);                           % Velocity(feedFlow(1), rho(1), table2cell(Parameters_table(:,3)) ) * ones(nstages,1)
            0;
            ];

Parameters          = Parameters_table{:,3};
Parameters(1:9)     = [nstages, C0solid, r, epsi, dp, L, rho_s, km, mi];

Parameters(which_k) = k0;
Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

[xx_0] = simulateSystem(F, [], x0, Parameters_opt  );

%%
ind = 1:numel(Time);

NAME = {'C_f','C_s','T','Continuity','Momentum'};

for i=0:numel(NAME)-1

    subplot(2,3,i+1)
    imagesc(Time,1:nstages,xx_0(i*nstages+1:(i+1)*nstages,:)); colorbar
    hold on
    yline(nstagesbed(1),'w')
    yline(nstagesbed(end),'w')
    hold off
    title(NAME{i+1})
    
    end

subplot(2,3,6)
plotyy(Time, xx_0(end,:), Time, 1e3 * (sum(xx_0(0*nstages+1:1*nstages,ind) .* V_fluid) + sum(xx_0(1*nstages+1:2*nstages,ind) .* V_solid)) + xx_0(Nx,ind))

%%
xx_0(end,end)
