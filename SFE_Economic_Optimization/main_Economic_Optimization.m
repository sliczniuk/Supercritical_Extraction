startup;
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

Parameters_table        = readtable('Parameters.csv') ;                                 % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                              % Parameters within the model 

%% Load paramters
m_total                 = 200*80;

%V_Flow                  = 0.4;                                                          % Volumetric flow rate l/min

% Bed geometry
before                  = 0.05;                                                         % Precentage of length before which is empty
bed                     = 0.90;                                                         % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 150;                                                          % Minutes
timeStep                = 1;                                                            % Minutes
SamplingTime            = 5;                                                            % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                                % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;                        % Seconds
Time                    = [0 Time_in_sec/60];                                           % Minutes

N_Time                  = length(Time_in_sec);

%% Specify parameters to estimate
nstages                 = Parameters{1};

nstagesbefore           = 1:floor(before*nstages);
nstagesbed              = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
nstagesafter            = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

%% Number of variables
Nx                      = 3 * nstages+2;                                                % 3*Nstages(C_f, C_s, H) + P(t) + yield
Nu                      = 3 + numel( Parameters );                                      % T_in, P, F + numel(Parameters)

%% Extractor geometry
r                       = Parameters{3};                                                % Radius of the extractor  [m]
epsi                    = Parameters{4};                                                % Fullness [-]
L                       = Parameters{6};                                                % Total length of the extractor [m]

L_nstages               = linspace(0,L,nstages);
V                       = L  * pi * r^2;                                                % Total volume of the extractor [m3]
A                       = pi *      r^2;                                                % Extractor cross-section

%--------------------------------------------------------------------
V_slice                 = (L/nstages) * pi * r^2;

V_before                = V_slice * numel(nstagesbefore);
V_after                 = V_slice * numel(nstagesafter);
V_bed                   = V_slice * numel(nstagesbed);                                  % Volume of the fixed bed [m3]

V_before_solid          = repmat(V_before * 0          / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_solid             = repmat(V_bed    * epsi       / numel(nstagesbed)   , numel(nstagesbed)   ,1);
V_after_solid           = repmat(V_after  * 0          / numel(nstagesbed)   , numel(nstagesafter) ,1);

V_solid                 = [V_before_solid; V_bed_solid; V_after_solid];

V_before_fluid          = repmat(V_before * 1          / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_fluid             = repmat(V_bed    * (1 - epsi) / numel(nstagesbed)   , numel(nstagesbed)   ,1);
V_after_fluid           = repmat(V_after  * 1          / numel(nstagesafter) , numel(nstagesafter) ,1);

V_fluid                 = [V_before_fluid; V_bed_fluid; V_after_fluid];

L_bed_after_nstages     = L_nstages(nstagesbed(1):end);
L_bed_after_nstages     = L_bed_after_nstages - L_bed_after_nstages(1);
L_end                   = L_bed_after_nstages(end);

%% symbolic variables
x                       = MX.sym('x', Nx);
u                       = MX.sym('u', Nu);

%% Set Integrator
f                       = @(x, u) modelSFE(x, u, bed_mask, timeStep_in_sec);
xdot                    = modelSFE(x, u, bed_mask, timeStep_in_sec);
F                       = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

%% Set operating conditions
feedPress               = 200;                                                          % pressure in the extractor [bar]
massFlow                = 0.02;                                                         % mass flow rate in kg/s

%% compressor
[T_out_compressor, W_com, Comp_cost] = Pump_estimation(10+273, 70, feedPress, massFlow, Parameters);

%% heat exchanger
[T_w_out, T_c_out, HX_cost] = Heat_Exchanger_estimation(T_out_compressor, feedPress, Parameters);
T0homog                     = full(T_w_out);            

%% Extractor
[Extractor_cost]        = Extractor_estimation(feedPress,Parameters);

%% Controls
Z                       = Compressibility( T0homog, feedPress,         Parameters );
rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );

enthalpy_rho            = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );

feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;                 % Kelvin

feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;                 % Bars

feedFlow                = massFlow * ones(1,length(Time_in_sec));                       % kg/s

uu                      = [feedTemp', feedPress', feedFlow'];

%% Re range: 584(40C, 300bar) < 631(50C, 300bar) < 718(40C 200bar) < 795(50C, 200bar) for empty pipe acc to Aspen
VELOCITY                = Velocity(massFlow, rho, Parameters);
mu                      = Viscosity(T0homog, rho);
Re                      = (2*r) .* VELOCITY ./ (mu / rho);
keyboard

%% Set inital state and inital conditions
msol_max                = m_total;                                                      % g of product in solid and fluid phase
mSol_ratio              = 1;

mSOL_s                  = msol_max*mSol_ratio;                                          % g of product in biomass
mSOL_f                  = msol_max*(1-mSol_ratio);                                      % g of biomass in fluid

C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;                            % Solid phase kg / m^3
Parameters{2}           = C0solid;

G                       =@(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;

m_fluid                 = G(L_bed_after_nstages)*( L_bed_after_nstages(2) );            % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
m_fluid                 = [zeros(1,numel(nstagesbefore)) m_fluid];
C0fluid                 = m_fluid * 1e-3 ./ V_fluid';

% Initial conditions
x0                      = [ C0fluid'                        ;
                            C0solid  * bed_mask             ;
                            enthalpy_rho * ones(nstages,1)  ;
                            feedPress(1)                    ;
                            0                               ;
                            ];

%% Run plain simulation
Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];
[xx_0]                 = simulateSystem(F, [], x0, Parameters_init_time );
