%clc, close all
%clear all

startup;
%p = Pushbullet(pushbullet_api);

%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

%% Load paramters
DATA_set                = {'LUKE_T40_P200', 'LUKE_T50_P200', 'LUKE_T40_P300', 'LUKE_T50_P300'};
Parameters_table        = readtable('Parameters.csv') ;        % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});
DATA                    = DATA_set{1};

%% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 150;
simulationTime          = PreparationTime + ExtractionTime;

timeStep                = 1/2;                                                 % Minutes

timeStep_in_sec         = timeStep * 60;                                       % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                    = [0 Time_in_sec/60];                                  % Minutes

N_Time                  = length(Time_in_sec);

%% Specify parameters to estimate
nstages                 = Parameters{1};

% Bed geometry
before                  = 0.1;          
bed                     = 0.165;        

nstagesbefore           = 1:floor(before*nstages);
nstagesbed              = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
nstagesafter            = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

%% Number of variables
Nx                      = 3 * nstages+1;
Nu                      = 3 + numel( Parameters );

%% symbolic variables
x                       = MX.sym('x', Nx);
u                       = MX.sym('u', Nu);

%% Set Integrator
f                       = @(x, u) modelSFE_uniform_U(x, u, bed_mask);

% Integrator
F                       = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

%V                       = 0.010;                                       % deafault = 0.01
r                       = Parameters{3};                                % Radius of the extractor  [m]
epsi                    = Parameters{4};                                % Fullness [-]
L                       = Parameters{6};                                % Total length of the extractor [m]
V_Flow                  = 0.39;

L_nstages               = linspace(0,L,nstages);
V                       = L  * pi * r^2;                                 % Total volume of the extractor [m3]
A                       = pi *      r^2;                                 % Extractor cross-section

%%
LabResults              = xlsread([DATA,'.xlsx']);

T0homog                 = LabResults(1,1)+273.15;
feedPress               = LabResults(1,2);

data_org                = LabResults(:,5)';
data                    = diff(data_org);

%% Set parameters
mSOL_s                   = 78;                                          % g of product in biomass
mSOL_f                   = 0;                                           % g of biomass in fluid

%C0fluid                 = 1;                                           % Extractor initial concentration of extract - Fluid phase kg / m^3

%--------------------------------------------------------------------
V_slice                 = (L/nstages) * pi * r^2;

V_before                = V_slice * numel(nstagesbefore);
V_after                 = V_slice * numel(nstagesafter);
V_bed                   = V_slice * numel(nstagesbed);                  % Volume of the fixed bed [m3]

V_before_solid          = repmat(V_before * 0          / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_solid             = repmat(V_bed    * epsi       / numel(nstagesbed)   , numel(nstagesbed)   ,1);
V_after_solid           = repmat(V_after  * 0          / numel(nstagesbed)   , numel(nstagesafter) ,1);

V_solid                 = [V_before_solid; V_bed_solid; V_after_solid];

V_before_fluid          = repmat(V_before * 1          / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_fluid             = repmat(V_bed    * (1 - epsi) / numel(nstagesbed)   , numel(nstagesbed)   ,1);
V_after_fluid           = repmat(V_after  * 1          / numel(nstagesafter) , numel(nstagesafter) ,1);

V_fluid                 = [V_before_fluid; V_bed_fluid; V_after_fluid];

%--------------------------------------------------------------------
                
C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;            % Solid phase kg / m^3

C0fluid                 = mSOL_f * 1e-3 / (V_before + V_bed * (1-epsi) + V_after);

%%
Z                       = Compressibility( T0homog, feedPress,         Parameters );
rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );

enthalpy                = SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters);                

% Set operating conditions
feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;  % Kelvin
%feedTemp(round(numel(Time)/10):round(numel(Time)/4))   = feedTemp(1) + 10;

feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;  % Bars
%feedPress(round(numel(Time)/10):round(numel(Time)/4))   = feedPress(1) - 100;

feedFlow                = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> km3/s
%feedFlow(1:N_Sample(1)) = 0;% linspace(feedFlow(1)/10,feedFlow(1),numel(feedFlow(1:N_Sample(1))));

uu                      = [feedTemp', feedPress', feedFlow'];

% Initial conditions
x0                      = [ C0fluid * ones(nstages,1);
                            C0solid * bed_mask;
                            T0homog * ones(nstages,1);
                            0;
                            ];

%%
%Parameters(which_k) = DATA_K_OUT(1:numel(which_k),ii);
Parameters_opt = [uu repmat(cell2mat(Parameters),1,N_Time)'];
[xx_out] = simulateSystem(F, [], x0, Parameters_opt  );

%%
subplot(2,2,1);imagesc(xx_out(1:250,:))
subplot(2,2,2);imagesc(xx_out(251:500,:))
subplot(2,2,3);imagesc(xx_out(501:750,:))
subplot(2,2,4);plot(xx_out(end,:))
