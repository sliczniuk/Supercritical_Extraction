clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

%% load parameters and set number of stages
Parameters_table     = readtable('Parameters.csv') ;        % Fulle table with prameters
Parameters           = Parameters_table{:,3};               %
Parameters_sym       = MX(Parameters_table{:,3});           % Vector of paraneters in the form casadi vector

%% Specify parameters to estimate
nstages              = 50;
which_k              = [8, 44, 45];
k0                   = [1, 100, 0];

%%
Nx                   = 3*nstages+1;
Nu                   = 3 + numel(Parameters);
Nk                   = numel(which_k);

%% Set time of the simulation
ExtractionTime       = 150;                                                 % Minutes
PreparationTime      = 0;                                                   % Minutes
simulationTime       = PreparationTime + ExtractionTime;                    % Minutes

timeStep             = 1/2;                                                 % Minutes
SamplingTime         = 5;                                                   % Minutes

timeStep_in_sec      = timeStep * 60;                                       % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                 = [0 Time_in_sec/60];                                  % Minutes

N_Time               = length(Time_in_sec);


%% symbolic variables

% Create symbolic variables
x                       = MX.sym('x', Nx);
u                       = MX.sym('u', Nu);

% Create the solver
OPT_solver              = casadi.Opti();
nlp_opts                = struct;
nlp_opts.ipopt.max_iter = 50;
ocp_opts                = {'nlp_opts', nlp_opts};
OPT_solver.solver(         'ipopt'   , nlp_opts)

% Descision variables
k                       = OPT_solver.variable(Nk);
k_lu                    = [ zeros(Nk,1) , inf*ones(Nk,1) ];
% Constraints
for nk=1:Nk
    OPT_solver.subject_to( k_lu(nk,1) <= k(nk,:) <= k_lu(nk,2) );
end

%% Set parameters
m_ref                   = 2.1;                                          % g of product obtained from a biomass

C0fluid                 = 0;                                            % Extractor initial concentration of extract
% Fluid phase kg / m^3

V                       = 12e-5;                                        % Volume of the extractor  [m3]
L                       = 0.09;                                         % Length of the extractor [m]
epsi                    = 0.4;                                          % Porosity [-]

C0solid                 = m_ref *1e-3 / (V*(1-epsi));                   % Solid phase kg / m^3

dp                      = 0.85 * 1e-3;                                  % Diameter of the particle [m] - Vargas
rho_s                   = 1800.0;                                       % Densisty of the solid phase [kg / m^3] - FC
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
T0homog   = 40+273.15;
feedPress = 300;

% Set operating conditions
feedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars

feedFlow  = 1 * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec

uu = [feedTemp', feedPress', feedFlow'];

% Inital conditions
x0 = [C0fluid*ones(nstages,1);
    C0solid*ones(nstages,1);
    T0homog*ones(nstages,1);
    0;
    ];

%% DATA
data1          = [0.06052631578947365, 0.23245614035087725, 0.3587719298245615, 0.47807017543859665, 0.6043859649122808, 0.7473684210526317, 0.8464912280701755, 0.9070175438596493, 0.9438596491228072, 0.955263157894737];
SAMPLE1        = [11, 28, 40, 50, 62, 75, 89, 102, 122, 142];
    
data2          = [0.08859649122807023, 0.1885964912280702, 0.30175438596491233, 0.5166666666666668, 0.7228070175438598, 0.8640350877192984, 0.9228070175438599, 0.9385964912280703, 0.9508771929824563, 0.9622807017543862, 0.9728070175438598, 0.9815789473684212, 0.9894736842105265];
SAMPLE2        = [16, 21, 25, 33, 38, 43, 48, 53, 60, 65, 72, 77, 83];
    
DATA           = {data1, data2};
SAMPLE         = {SAMPLE1, SAMPLE2};
    
%% Solve and simulate
KOUT = [];
for i=1:numel(DATA)

    data   = DATA{i};
    sample = SAMPLE{i};

    N_Sample   = [];
    for i = 1:numel(sample)
        N_Sample = [N_Sample ; find(Time == sample(i))];
    end

    if numel(N_Sample) ~= numel(data)
        keyboard
    end

    % Store symbolic results of the simulation
    X = MX(Nx,N_Time+1);
    X(:,1) = x0;

    % Symbolic integration
    tic
    for j=1:N_Time
        X(:,j+1)=F(X(:,j), [uu(j,:)'; Parameters_sym] );
    end
    toc

    % Find the measurment from the simulation
    Yield_estimate = X(end,N_Sample) / m_ref;

    % Create the cost function
    J = (data-Yield_estimate ) * diag(1:1:1) * (data-Yield_estimate )';

    %
    % Set opt and inital guess
    OPT_solver.minimize(J);

    OPT_solver.set_initial(k, k0);

    % Solve the opt
    tic
    try
        sol = OPT_solver.solve();
        kout = full(sol.value(k))
    catch
        kout = OPT_solver.debug.value(k);
    end
    toc

    KOUT = [KOUT, kout];

end
%%
id = 1;
for i=1:numel(DATA)

    data   = DATA{i};
    sample = SAMPLE{i};

    kout = KOUT(:,i);

    N_Sample   = [];
    for i = 1:numel(sample)
        N_Sample = [N_Sample ; find(Time == sample(i))];
    end

    % Simulate system with obtained parameters

    Parameters(which_k) = k0;
    Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

    [xx_0] = simulateSystem(F, [], x0, Parameters_opt  );

    Parameters(which_k) =KOUT;
    Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

    [xx_out] = simulateSystem(F, [], x0, Parameters_opt  );

    %
    subplot(numel(DATA),3,id)
    imagesc(Time,1:nstages,xx_out(1:nstages,:)); colorbar
    xlabel('Time [min]')
    ylabel('Stages')
    title('Fluid Concentration')
    ylabel(['T=',mat2str(T0homog),'P=',mat2str(feedPress(1))])
    axis tight

    subplot(numel(DATA),3,id+1)
    imagesc(Time,1:nstages,xx_out(1*nstages+1:2*nstages,:)); colorbar
    xlabel('Time [min]')
    ylabel('Stages')
    title('Solid Concentration')
    axis tight

    subplot(numel(DATA),3,id+2)
    hold on
    plot(Time,xx_out(end,:)/m_ref)
    plot(Time,xx_0(end,:)/m_ref,'--')
    plot(sample,data,'o','MarkerSize',5)
    hold off

    legend('Estiamted','Inital','Data')
    legend off
    xlabel('Time [min]')
    ylabel('Mass of the extract [g]')
    title(['k=',sprintf('%0.3g, ',kout)])
    axis tight

    id=id+3;

end

%%
%set(gcf,'PaperOrientation','landscape')
%print(gcf, 'FitPreparationTime_orgV.pdf', '-dpdf','-r700','-fillpage')
%%
%}