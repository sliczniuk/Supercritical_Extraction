clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

DATA = {'LUKE_T40_P200.xlsx', 'LUKE_T50_P200.xlsx', 'LUKE_T40_P300.xlsx', 'LUKE_T50_P300.xlsx'};

%% load parameters and set number of stages
Parameters_table     = readtable('Parameters.csv') ;        % Fulle table with prameters
Parameters           = Parameters_table{:,3};               % 
Parameters_sym       = MX(Parameters_table{:,3});           % Vector of paraneters in the form casadi vector

%% Specify parameters to estimate
nstages              = 5;       
which_k              = [8  ,44];
k0                   = [0.3, 1];

%% 
Nx                   = 3*nstages+1;
Nu                   = 3 + numel(Parameters);
Nk                   = numel(which_k);

%% Set time of the simulation
simulationTime       = 150;                                                 % Minutes
timeStep             = 1/2;                                                 % Minutes
SamplingTime         = 5;                                                   % Minutes
delayTime            = 0;                                                   % Minutes

timeStep_in_sec      = timeStep * 60;                                       % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                 = [0 Time_in_sec/60];                                  % Minutes

N_Sample             = find(Time == SamplingTime) - 1;
N_Delay              = delayTime / timeStep;
N_Time               = length(Time_in_sec);

%% symbolic variables

% Create symbolic variables
x                       = MX.sym('x', Nx);
u                       = MX.sym('u', Nu);

% Store the results of the optimization
KOUT                    = nan(Nk,numel(DATA));

% Create the solver
OPT_solver              = casadi.Opti();
nlp_opts                = struct;
nlp_opts.ipopt.max_iter = 20;
ocp_opts                = {'nlp_opts', nlp_opts};
OPT_solver.solver(         'ipopt'   , nlp_opts)

% Descision variables
k                       = OPT_solver.variable(Nk);
k_lu                    = [ zeros(Nk,1) , inf*ones(Nk,1) ];

%% Set parameters
m_ref                   = 78;                                           % g of product obtained from a kg of biomass
                     
C0fluid                 = 0;                                            % Extractor initial concentration of extract
                                                                         % Fluid phase kg / m^3
                     
V                       = 0.00165;                                      % Volume of the extractor  [m3] 
L                       = 0.095;                                        % Length of the extractor [m]
epsi                    = 2/3;                                          % Porosity [-] 
                 
C0solid                 = m_ref *1e-3 / (V*(1-epsi));                   % Solid phase kg / m^3

dp                      = 0.00010;                                      % Diameter of the particle [m] - Vargas
rho_s                   = 1300.0;                                       % Densisty of the solid phase [kg / m^3] - FC
km                      = 0.29;                                         % Partition coefficient (?)
             
mi                      = 1/2;                                          % Geometric shape coefficient (?) 

% Assign new values of parameters to the Parameters vector 
%                       nstages, C0solid, V, epsi, dp, L, rho_s, km, mi                
Parameters(1:9)         = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];
Parameters_sym(1:9)     = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];

% Decide which parameters are decision variabales
Parameters_sym(which_k) = k;


%% Set Integrator
tic
%f_r = @(x, u, k) modelSFE_Regression(x, u, k, parameters);
f = @(x, u) modelSFE(x, u, nstages);

% Integrator
F = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);
toc

%% Loop over datasets

i = 4;

% load dataset
LabResults = xlsread(DATA{i});

T0homog   = LabResults(1,1)+273.15;
feedPress = LabResults(1,2);
rho       = LabResults(1,3);
data      = LabResults(:,5)';

V_Flow = 0.39;

% Set operating conditions
feedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars

feedFlow  = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec

uu = [feedTemp', feedPress', feedFlow'];

% Inital conditions
x0 = [C0fluid*ones(nstages,1);
    C0solid*ones(nstages,1);
    T0homog*ones(nstages,1);
    0];

% Store symbolic results of the simulation
X = MX(Nx,N_Time+1);
X(:,1) = x0;


%% Symbolic integration
tic
for j=1:N_Time
    X(:,j+1)=F(X(:,j), [uu(j,:)'; Parameters_sym] );
end
toc

%% Find the measurment from the simulation
Yield_estimate = X(3*nstages+1,:);
Yield_estimate = [zeros(1,N_Delay) Yield_estimate(1:end-N_Delay)];
Yield_estimate = Yield_estimate(1:N_Sample:end);

%% Create the cost function
J = (data-Yield_estimate ) * diag(1:1:1) * (data-Yield_estimate )';
fJ = Function('fJ', {k}, {J} );

%% Constraints
for nk=1:Nk
    OPT_solver.subject_to( k_lu(nk,1) <= k(nk,:) <= k_lu(nk,2) )
end

%% Set opt and inital guess
OPT_solver.minimize(J);

OPT_solver.set_initial(k, k0);

%% Solve the opt
tic
try
    sol = OPT_solver.solve();
    kout = sol.value(k)
catch
    kout = OPT_solver.debug.value(k);
end
toc

%% Simulate system with obtained parameters
Parameters(which_k) = kout;
Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

[xx_out] = simulateSystem(F, [], x0, Parameters_opt  );

%%
hold on
plot(Time,xx_out(end,:))
plot(Time(1:N_Sample:end),data,'o')
hold off
