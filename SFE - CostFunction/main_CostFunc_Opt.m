clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

DATA = {'LUKE_T40_P200.xlsx', 'LUKE_T50_P200.xlsx', 'LUKE_T40_P300.xlsx', 'LUKE_T50_P300.xlsx'};
%DATA = {'LUKE_T50_P200.xlsx'};

Parameters_table     = readtable('Parameters.csv') ;        % Fulle table with prameters

%% Set time of the simulation
simulationTime       = 150;                                                 % Minutes
SamplingTime         = 5;                                                   % Minutes

timeStep             = 1/2;                                                 % Minutes

timeStep_in_sec      = timeStep * 60;                                       % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                 = [0 Time_in_sec/60];                                  % Minutes

SAMPLE   = [0:SamplingTime:simulationTime];
SAMPLE(1) = 0;

N_Sample = [];
for i = 1:numel(SAMPLE)
    N_Sample = [N_Sample ; find(Time == SAMPLE(i))];
end
if numel(N_Sample) ~= numel(SAMPLE)
    keyboard
end

N_Time               = length(Time_in_sec);

%% Specify parameters to estimate
nstages              = 10;
which_k              = [8, 44];

%% Set parameters
m_ref                = 78;                                           % g of product obtained from a kg of biomass

V                    = 0.010;                                    % Volume of the extractor  [m3]
L                    = 0.6;                                       % Length of the extractor [m]
epsi                 = 0.90;                                      % Porosity [-]

a=1;
C0solid              =a*m_ref *1e-3 / (V*(1-epsi));            % Solid phase kg / m^3
C0fluid              =(1-a)*m_ref *1e-3 / (V*(epsi));            % Extractor initial concentration of extract - Fluid phase kg / m^3

dp                   = 0.00010;                                      % Diameter of the particle [m]
rho_s                = 1250.0;                                       % Densisty of the solid phase [kg / m^3] - FC
km                   = 0.29;                                         % Partition coefficient (?)

mi                   = 1/2;                                          % Geometric shape coefficient (?)

V_Flow               = 0.39;

Parameters = Parameters_table{:,3};
Parameters(1:9)     = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];

%% Loop over datasets

kp_Set = [0.03]; %[0.01, 0.1, 0.5, 1, 2];
Di_Set = [0.2]; %[0.01, 0.1, 1, 2];

KOU      = nan( numel(which_k), numel(kp_Set), numel(Di_Set), numel(DATA) );
JJ       = nan(                 numel(kp_Set), numel(Di_Set), numel(DATA) );
J_STATUS = cell(                numel(kp_Set), numel(Di_Set), numel(DATA) );

%%
Nx                   = 3*nstages+1;
Nu                   = 3 + numel( Parameters_table{:,3} );
Nk                   = numel(which_k);

%% symbolic variables

% Create symbolic variables
x                       = MX.sym('x', Nx);
u                       = MX.sym('u', Nu);

%% Set Integrator
f = @(x, u) modelSFE(x, u, nstages);

% Integrator
F = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

%%

    for d=1:numel(Di_Set)

        %clc; [i, d]

        for kp=1:numel(kp_Set)

            parfor id=1:numel(DATA)
        
             % load dataset
             LabResults = xlsread(DATA{id});
        
             T0homog   = LabResults(1,1)+273.15;
             feedPress = LabResults(1,2);
             rho       = LabResults(1,3);
             data      = LabResults(:,5)';
             data      = diff(data);
        
             % Set operating conditions
             feedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
             feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars
        
             feedFlow  = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec
        
             uu = [feedTemp', feedPress', feedFlow'];
        
             % Inital conditions
             x0 = [C0fluid*ones(nstages,1);
                 C0solid*ones(nstages,1);
                 T0homog*ones(nstages,1);
                 0;
                 ];
        
            k0 = [kp_Set(kp), Di_Set(d)]
            
            %% load parameters and set number of stages
            Parameters_sym       = MX(Parameters_table{:,3});           % Vector of paraneters in the form casadi vector

            % Create the solver
            OPT_solver                  = casadi.Opti();
            
            nlp_opts                    = struct;
            nlp_opts.ipopt.max_iter     = 10;
            %nlp_opts.ipopt.max_cpu_time = 300;
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

            % store the results
            KOUT(:,kp,d,id) = kout;
            JJ(kp,d,id) = full(fJ(kout));
            J_STATUS{kp,d,id} = OPT_solver.return_status;

            %clc

            end

        end

        %save Fit_Di_km_Data2.mat

    end


%% Set time of the simulation
%{
ExtractionTime       = 40;                                                 % Minutes
simulationTime       = PreparationTime + ExtractionTime;                    % Minutes

timeStep_in_sec      = timeStep * 60;                                       % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                 = [0 Time_in_sec/60];                                  % Minutes

N_Sample             = find(Time == SamplingTime) - 1;
N_Preparation        = find(Time == PreparationTime) - 1;

N_Time               = length(Time_in_sec);

%% Set Integrator
tic
%f_r = @(x, u, k) modelSFE_Regression(x, u, k, parameters);
f = @(x, u) modelSFE(x, u, nstages);

% Integrator
F = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);
toc
%}

%% Simulate system with obtained parameters

    id = 1;
    for i = 1:numel(DATA)
    
        % load dataset
        LabResults = xlsread(DATA{i});
    
        T0homog   = LabResults(1,1)+273.15;
        feedPress = LabResults(1,2);
        rho       = LabResults(1,3);
        data      = LabResults(:,5)';
    
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
    
        kout  = KOUT(:,1,1,i);

        k0 = [0.03, 0.2]; 
    
        %
        Parameters(which_k) = k0;
        Parameters_opt = [uu repmat(Parameters,1,N_Time)'];
    
        [xx_0] = simulateSystem(F, [], x0, Parameters_opt  );
    
        Parameters(which_k) =kout;
        Parameters_opt = [uu repmat(Parameters,1,N_Time)'];
    
        [xx_out] = simulateSystem(F, [], x0, Parameters_opt  );
    
        %
        %figure(i)
        set(gcf,'PaperOrientation','landscape')

        subplot(numel(DATA),3,id)
        imagesc(Time,1:nstages,xx_out(1:nstages,:)); colorbar
        xlabel('Time [min]')
        ylabel('Stages')
        title('Fluid Concentration')
        ylabel(['T=',mat2str(T0homog),', P=',mat2str(feedPress(1))])
        axis tight
    
        subplot(numel(DATA),3,id+1)
        imagesc(Time,1:nstages,xx_out(1*nstages+1:2*nstages,:)); colorbar
        xlabel('Time [min]')
        ylabel('Stages')
        title('Solid Concentration')
        axis tight
    
        subplot(numel(DATA),3,id+2)
    
        hold on
        plot(Time,xx_out(end,:))
        plot(Time,xx_0(end,:),'--')
        plot(SAMPLE,data,'o','MarkerSize',5)
        hold off
    
        legend('Estiamted','Inital','Data')
        legend off
        xlabel('Time [min]')
        ylabel('Mass of the extract [g]')
        %title(['k=',sprintf('%0.3g, %0.3g, %0.3g, Cost Function=%0.3g, \n STATUS=%s',kout, JJ(v,i), J_STATUS{v,i} )  ])
        axis tight
        ylim([0, m_ref+2])

        %print(gcf, ['Bounded_Fit_V_',mat2str(V),'_Kp_k0_',mat2str(k0),'.pdf'], '-dpdf','-r700','-fillpage')
    
        id = id + 3;
    
    end

%%

%%
%}