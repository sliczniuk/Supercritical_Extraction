clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

DATA = {'VARGAS_T40_P90.xlsx','VARGAS_T50_P90.xlsx','VARGAS_T60_P90.xlsx','VARGAS_T70_P90.xlsx'};
DATA_to_check = [4];

%% load parameters and set number of stages
Parameters_table     = readtable('Parameters.csv') ;        % Fulle table with prameters
Parameters           = Parameters_table{:,3};               % 
Parameters_sym       = MX(Parameters_table{:,3});           % Vector of paraneters in the form casadi vector

%% Specify parameters to estimate
nstages              = 50;       
which_k              = [8, 44, 45];
K0                   = [0.8133, 3.03e+1, 0; 
                        0.6944, 3.66e+1, 0; 
                        0.1041, 1.97e+2, 0; 
                        0.0667, 2.73e+2, 0];
%Dx_0                 = 0;

%% 
Nx                   = 3*nstages+1;
Nu                   = 3 + numel(Parameters);
Nk                   = numel(which_k);

%% Set time of the simulation
ExtractionTime       = 60;                                                 % Minutes
PreparationTime      = 0;                                                   % Minutes
simulationTime       = PreparationTime + ExtractionTime;                    % Minutes

timeStep             = 1/2;                                                 % Minutes
SamplingTime         = 5;                                                   % Minutes

timeStep_in_sec      = timeStep * 60;                                       % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                 = [0 Time_in_sec/60];                                  % Minutes

SAMPLE               = [0, 3, 6, 9, 12, 15, 20, 30, 40, 60];
N_Sample             = [];
for i = 1:numel(SAMPLE)
    N_Sample = [N_Sample ; find(Time == SAMPLE(i))];
end

N_Preparation        = find(Time == PreparationTime) - 1;

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
C0fluid = 0;                                            % Extractor initial concentration of extract
                                                        % Fluid phase kg / m^3

V       = 4.173e-6;                                     % Volume of the extractor  [m3] - Vargas
C0solid = 25;                                           % Solid phase kg / m^3
epsi    = 0.669;                                        % Porosity [-] - Vargas
dp      = 5e-4;                                         % Diameter of the particle [m] - Vargas
L       = 0.045;                                        % Length of the extractor [m] - Vargas
rho_s   = 1087.2;                                       % Densisty of the solid phase [kg / m^3] - Vargas
km      = 0.29;                                         % Partition coefficient (?)
mi      = 1/5;                                          % Geometric shape coefficient (?) - Vargas

% Assign new values of parameters to the Parameters vector 
%                       nstages, C0solid, V, epsi, dp, L, rho_s, km, mi                
Parameters(1:9)         = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];
Parameters_table{1:9,3} = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi]';
Parameters_sym(1:9)     = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];

% Decide which parameters are decision variabales

Parameters_sym(which_k) = k;

V_Flow = 3.34e-8;                                                       % m3/s

%% Set Integrator
tic
%f_r = @(x, u, k) modelSFE_Regression(x, u, k, parameters);
f = @(x, u) modelSFE(x, u, nstages);

% Integrator
F = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);
toc

%% Loop over datasets

%km_range = linspace(0.01, 2, 200);
%Di_range = linspace(0.01, 1, 204);
%cJ       = nan(numel(km_range),numel(Di_range),4);
KOUT     = nan(Nk,4);

%%

for i=DATA_to_check
    clc
    i

    k0 = [K0(i,:)];
    % load dataset
    LabResults = xlsread(DATA{i});
    
    T0homog   = LabResults(1,1)+273.15;
    feedPress = LabResults(1,2);
    rho       = LabResults(1,3);
    data      = LabResults(:,5)' *1e-2;
    %data      = data(1:(N_Time-N_Preparation)/N_Sample+1);
        
    % Set operating conditions
    feedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars
    
    feedFlow  = V_Flow * rho * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec
    feedFlow(1:N_Preparation) = 0;
    
    uu = [feedTemp', feedPress', feedFlow'];
    
    % Inital conditions
    x0 = [C0fluid*ones(nstages,1);
        C0solid*ones(nstages,1);
        T0homog*ones(nstages,1);    
        0;
        ];
%{
    Parameters(which_k) = k0;
    Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

    [xx_0] = simulateSystem(F, [], x0, Parameters_opt  );

    subplot(3,1,1)
    imagesc(Time,1:nstages,xx_0(1:nstages,:)); colorbar

    subplot(3,1,2)
    imagesc(Time,1:nstages,xx_0(nstages+1:2*nstages,:)); colorbar

    subplot(3,1,3)
    hold on
    plot(Time,xx_0(end,:))
    plot(SAMPLE,data ,'o' )
    plot(SAMPLE,xx_0(end,N_Sample),'*')
    hold off

end
%}

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
    Yield_estimate = X(end,N_Sample) / 0.0015;
    %Yield_estimate = modelSFE_out2(X, nstages, rho, feedFlow(1), timeStep_in_sec);

    %Yield_estimate = Yield_estimate(N_Sample)   ;
    
    % Create the cost function
    J = (data-Yield_estimate ) * diag(1:1:1) * (data-Yield_estimate )';

    %{
    fJ = Function('fJ', {k}, {J} );

    % Brute force    
    for kk = 1:numel(km_range)
        tic
        parfor dd = 1:numel(Di_range)
            cJ(kk,dd,i) = full(fJ([km_range(kk), Di_range(dd),1]));
        end
        toc
    end

    cJi = cJ(:,:,i);
    [r,c] = find(cJi == min( cJi(:) ));
    k0    = [km_range(r), Di_range(c),Dx_0];
    %}
    
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

    % store the results
    KOUT(:,i) = kout;
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
for i = DATA_to_check

    clc
    i

    k0 = [K0(i,:)];
    % load dataset
    LabResults = xlsread(DATA{i});
    
    T0homog   = LabResults(1,1)+273.15;
    feedPress = LabResults(1,2);
    rho       = LabResults(1,3);
    data      = LabResults(:,5)' *1e-2;
    %data      = data(1:(N_Time-N_Preparation)/N_Sample+1);
        
    % Set operating conditions
    feedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars
    
    feedFlow  = V_Flow * rho * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec
    feedFlow(1:N_Preparation) = 0;
    
    uu = [feedTemp', feedPress', feedFlow'];
    
    % Inital conditions
    x0 = [C0fluid*ones(nstages,1);
        C0solid*ones(nstages,1);
        T0homog*ones(nstages,1);    
        0;
        ];


    % find min value in the matrix
    %cJi = cJ(:,:,i);
    %[r,c] = find(cJi == min( cJi(:) ));
    %k0    = [km_range(r), Di_range(c)];
    kout  = KOUT(:,i);

    %
    Parameters(which_k) = k0;
    Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

    [xx_0] = simulateSystem(F, [], x0, Parameters_opt  );

    Parameters(which_k) =kout;
    Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

    [xx_out] = simulateSystem(F, [], x0, Parameters_opt  );

    %
    subplot(numel(DATA_to_check),3,id)
    imagesc(Time,1:nstages,xx_out(1:nstages,:)); colorbar
    xlabel('Time [min]')
    ylabel('Stages')
    title('Fluid Concentration')
    ylabel(['T=',mat2str(T0homog),'P=',mat2str(feedPress(1))])
    axis tight

    subplot(numel(DATA_to_check),3,id+1)
    imagesc(Time,1:nstages,xx_out(1*nstages+1:2*nstages,:)); colorbar
    xlabel('Time [min]')
    ylabel('Stages')
    title('Solid Concentration')
    axis tight

    %{
    subplot(numel(DATA_to_check),4,4*i-1)
    surf(Di_range,km_range,cJi,'EdgeColor','none','FaceAlpha',0.5); colorbar
    hold on
    scatter(Di_range(c),km_range(r),cJi(r,c),'filled','k','SizeData',10)
    scatter(kout(2),kout(1),full(fJ(kout)),'k','SizeData',10)
    hold off
    title(sprintf('k= %0.3g, %0.3g, %0.3g',kout))
    axis square tight
    xlabel('Di')
    ylabel('km')
    view(2);
    
    subplot(numel(DATA_to_check),4,id+2)
    imagesc(Time,1:nstages,xx_out(3*nstages+1:4*nstages,:)); colorbar
    xlabel('Time [min]')
    ylabel('Stages')
    title('Solid Concentration')
    axis tight
%}

    subplot(numel(DATA_to_check),3,id+2)
    
    hold on
    plot(Time,xx_out(end,:))
    plot(Time,xx_0(end,:),'--')
    plot(SAMPLE,data ,'o' )
    hold off

    legend('Estiamted','Inital','Data')
    legend off
    xlabel('Time [min]')
    ylabel('Mass of the extract [g]')
    title(['k=',sprintf('%0.3g, ',kout)])
    axis tight
    %ylim([0, m_ref+2])

    id = id + 3;

end
%%
%set(gcf,'PaperOrientation','landscape')
%print(gcf, 'FitPreparationTime_orgV.pdf', '-dpdf','-r700','-fillpage')
%%
%}