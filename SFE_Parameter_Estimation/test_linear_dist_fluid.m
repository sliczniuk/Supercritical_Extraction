%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');

clc; clear all

addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

Parameters_table        = readtable('Parameters.csv') ;        % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});

%% Load paramters
DATA_set                = {'LUKE_T40_P200', 'LUKE_T50_P200', 'LUKE_T40_P300_org', 'LUKE_T50_P300_org'};
%DATA_set                = {'LUKE_T50_P300'};

which_k                 = [8,  44, 45,          ];                          % Select which parameters are used for fitting
k0                      = [1,  3 , 1 ,  77, 0.65, 0.3];                          % Give inital value for each parameter
Nk                      = numel(which_k)+3;                                 % Parameters within the model + m_max, m_ratio, sigma
k_lu                    = [ [0;0;0;77;0;0] , [1;inf;inf;100;1;inf] ];

Iteration_max           = 100;                                              % Maximum number of iterations for optimzer
Time_max                = 24;                                               % Maximum time of optimization in [h]

V_Flow                  = 0.4;                                              % Volumetric flow rate l/min

% Bed geometry
before                  = 0.1;                                              % Precentage of length before which is empty
bed                     = 0.165;                                            % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 150;
timeStep                = 0.25;                                             % Minutes
SamplingTime            = 5;                                                % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);

SAMPLE                  = [PreparationTime:SamplingTime:simulationTime];

% Check if the number of data points is the same for both the dataset and the simulation
N_Sample                = [];
for i = 1:numel(SAMPLE)
    N_Sample            = [N_Sample ; find(round(Time,3) == round(SAMPLE(i))) ];
end
if numel(N_Sample) ~= numel(SAMPLE)
    keyboard
end

DATA_K_OUT              = nan(Nk,numel(DATA_set));                          % Store Parameters obatined from all fits (par num x num exper)

%% Specify parameters to estimate
nstages                 = Parameters{1}/2;

nstagesbefore           = 1:floor(before*nstages);
nstagesbed              = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
nstagesafter            = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

%% Number of variables
Nx                      = 3 * nstages+2;                    % 3*Nstages(C_f, C_s, H) + P(t) + yield
Nu                      = 3 + numel( Parameters );

%% Set Integrator
f                       = @(x, u) modelSFE(x, u, bed_mask, timeStep_in_sec);

% Integrator
F                       = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

r                       = Parameters{3};                                % Radius of the extractor  [m]
epsi                    = Parameters{4};                                % Fullness [-]
L                       = Parameters{6};                                % Total length of the extractor [m]

L_nstages               = linspace(0,L,nstages);
V                       = L  * pi * r^2;                                % Total volume of the extractor [m3]
A                       = pi *      r^2;                                % Extractor cross-section

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

L_bed_after_nstages         = L_nstages(nstagesbed(1):end);
L_bed_after_nstages         = L_bed_after_nstages - L_bed_after_nstages(1);
L_end                       = L_bed_after_nstages(end);

%% Set parameters
msol_max                    = 78;                                             % g of product in solid and fluid phase
mSol_ratio                  = 0.75;

mSOL_s                      = msol_max*mSol_ratio;                                        % g of product in biomass
mSOL_f                      = msol_max*(1-mSol_ratio);                                    % g of biomass in fluid

%--------------------------------------------------------------------

C0solid                     = mSOL_s * 1e-3 / ( V_bed * epsi)  ;            % Solid phase kg / m^3

G                           =@(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;

m_fluid                     = G(L_bed_after_nstages)*( L_bed_after_nstages(2) );                 % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
m_fluid                     = [zeros(1,numel(nstagesbefore)) m_fluid];         
C0fluid                     = m_fluid * 1e-3 ./ V_fluid';