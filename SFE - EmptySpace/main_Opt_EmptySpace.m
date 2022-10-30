clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

DATA = {'LUKE_T40_P300.xlsx'};
Parameters_table        = readtable('Parameters.csv') ;        % Fulle table with prameters

%% Set time of the simulation
PreparationTime         = 5;                                                  % Minutes
ExtractionTime          = 150;                                                 % Minutes
simulationTime          = PreparationTime + ExtractionTime;

timeStep                = 1/4;                                                 % Minutes

timeStep_in_sec         = timeStep * 60;                                       % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                    = [0 Time_in_sec/60];                                  % Minutes

N_Time                  = length(Time_in_sec);

%%
SamplingTime            = 5;                                                   % Minutes

SAMPLE   = [PreparationTime:SamplingTime:simulationTime];
SAMPLE(1) = PreparationTime;

N_Sample = [];
for i = 1:numel(SAMPLE)
    N_Sample = [N_Sample ; find(Time == SAMPLE(i))];
end
if numel(N_Sample) ~= numel(SAMPLE)
    keyboard
end

%% Specify parameters to estimate
nstages                 = 100;

before  = 0.14;         nstagesbefore   = 1:floor(before*nstages);
bed     = 0.165;        nstagesbed      = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
nstagesafter    = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

which_k                 = [8, 44, 45];

%% Set parameters
mSOL_s                   = 78;                                           % g of product in biomass
mSOL_f                   = 0;                                           % g of biomass in fluid

%C0fluid                 = 1;                                           % Extractor initial concentration of extract - Fluid phase kg / m^3

V                       = 0.01;                                         %
r                       = 0.075;                                        % Radius of the extractor  [m3]
L                       = V / pi / r^2;                                 % Total length of the extractor [m]
A                       = pi*r^2;                                       % Extractor cross-section
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
Nx                      = 3*nstages+1;
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

k0         = [1, 0.1, 10];


LabResults = xlsread(DATA{1});
        
T0homog   = LabResults(1,1)+273.15;
feedPress = LabResults(1,2);

data_org  = LabResults(:,5)';
data      = diff(data_org);

%%

% Set operating conditions
feedTemp   = T0homog   * ones(1,length(Time_in_sec)) + 0 ;  % Kelvin
%feedTemp(round(numel(Time)/10):round(numel(Time)/4))   = feedTemp(1) - 20;

feedPress  = feedPress * ones(1,length(Time_in_sec)) + 0 ;  % Bars
%feedPress(round(numel(Time)/3):round(2*numel(Time)/3))   = feedPress(1) - 10;

feedFlow   = V_Flow * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> km3/s
feedFlow(1:N_Sample(1))   = 0 ;

uu         = [feedTemp', feedPress', feedFlow'];

% Initial conditions
x0         = [
            C0fluid * ones(nstages,1);
            C0solid * bed_mask;
            T0homog*ones(nstages,1);
            0;
            ];

Parameters          = Parameters_table{:,3};
Parameters(1:9)     = [nstages, C0solid, r, epsi, dp, L, rho_s, km, mi];

%% load parameters and set number of stages

Parameters_sym       = MX(Parameters_table{:,3});           % Vector of paraneters in the form casadi vector

% Create the solver
OPT_solver                  = casadi.Opti();

nlp_opts                    = struct;
%nlp_opts.ipopt.max_iter     = 10;
nlp_opts.ipopt.max_cpu_time = 4000;
ocp_opts                    = {'nlp_opts', nlp_opts};
OPT_solver.solver(             'ipopt'   , nlp_opts)

% Descision variables
k                       = OPT_solver.variable(Nk);
k_lu                    = [ zeros(Nk,1) , inf*ones(Nk,1) ];
% Constraints
for nk=1:Nk
    OPT_solver.subject_to( k_lu(nk,1) <= k(nk,:) <= k_lu(nk,2) );
end

%% Assign new values of parameters to the Parameters vector
%                       nstages, C0solid, V, epsi, dp, L, rho_s, km, mi
Parameters_sym(1:9)     = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];

% Decide which parameters are decision variabales

Parameters_sym(which_k) = k;

% Store symbolic results of the simulation
X = MX(Nx,N_Time+1);
X(:,1) = x0;

% Symbolic integration
for j=1:N_Time
    X(:,j+1)=F(X(:,j), [uu(j,:)'; Parameters_sym] );
end

%% Find the measurment from the simulation
Yield_estimate = diff(X(3*nstages+1,N_Sample));

%% Create the cost function
J = (data-Yield_estimate ) * diag(1:1:1) * (data-Yield_estimate )';
fJ = Function('fJ', {k}, {J} );

%% Set opt and inital guess
OPT_solver.minimize(J);

OPT_solver.set_initial(k, k0);

%% Solve the opt
try
    sol = OPT_solver.solve();
    kout = full(sol.value(k));
catch
    kout = OPT_solver.debug.value(k);

end

%% Simulate system
Parameters(which_k) = k0;
Parameters_opt = [uu repmat(Parameters,1,N_Time)'];
[xx_0] = simulateSystem(F, [], x0, Parameters_opt  );

Parameters(which_k) = kout;
Parameters_opt = [uu repmat(Parameters,1,N_Time)'];
[xx_out] = simulateSystem(F, [], x0, Parameters_opt  );
%%

hold on 
plot(Time, xx_0(end,:));
plot(Time, xx_out(end,:));
plot(SAMPLE, data_org,'o');
xline(5)
hold off