clc, close all
clear all
addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

%DATA = {'LUKE_T40_P200.xlsx', 'LUKE_T50_P200.xlsx', 'LUKE_T40_P300.xlsx', 'LUKE_T50_P300.xlsx'};
DATA = {'LUKE_T40_P300.xlsx'};

Parameters_table        = readtable('Parameters.csv') ;        % Fulle table with prameters

%% Set time of the simulation
simulationTime          = 150;                                                 % Minutes
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

%% Loop over datasets

kp_Set                  = [0.1];
Di_Set                  = [0.5];

KOU                     =  nan( numel(which_k), numel(kp_Set), numel(Di_Set), numel(DATA));
JJ                      =  nan(                 numel(kp_Set), numel(Di_Set), numel(DATA));
J_STATUS                = cell(                 numel(kp_Set), numel(Di_Set), numel(DATA));

%%
Nx                      = 4*nstages+1;
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
feedTemp   = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
feedTemp(round(numel(Time)/4):round(numel(Time)/2))   = feedTemp(1) + 50;

feedPress  = feedPress * ones(1,length(Time_in_sec));  % Bars
%feedPress(round(numel(Time)/3):round(2*numel(Time)/3))   = feedPress(1) + 50;

feedFlow   = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec
%feedFlow(round(numel(Time)/3):round(2*numel(Time)/3))   = 2*feedFlow(1) ;

uu         = [feedTemp', feedPress', feedFlow'];

% Initial conditions
x0         = [C0fluid * ones(nstages,1);
    C0solid * bed_mask;
    T0homog*ones(nstages,1);
    zeros(nstages,1);
    0;
    ];

Parameters          = Parameters_table{:,3};
Parameters(1:9)     = [nstages, C0solid, r, epsi, dp, L, rho_s, km, mi];

Parameters(which_k) = k0;
Parameters_opt = [uu repmat(Parameters,1,N_Time)'];

[xx_0] = simulateSystem(F, [], x0, Parameters_opt  );

%%
ind = 1:numel(Time);

figure(1)
subplot(2,2,1)
imagesc(Time,1:nstages,xx_0(1*nstages+1:2*nstages,:)); colorbar
hold on
yline(nstagesbed(1),'w')
yline(nstagesbed(end),'w')
hold off
title('C_s')

subplot(2,2,2)
imagesc(Time,1:nstages,xx_0(0*nstages+1:1*nstages,:)); colorbar
hold on
yline(nstagesbed(1),'w')
yline(nstagesbed(end),'w')
hold off
title('C_f')

subplot(2,2,3)
imagesc(Time,1:nstages,xx_0(2*nstages+1:3*nstages,:)); colorbar;
hold on
yline(nstagesbed(1),'w')
yline(nstagesbed(end),'w')
hold off
title('T')

subplot(2,2,4)
imagesc(Time,1:nstages,xx_0(3*nstages+1:4*nstages,:)); colorbar;
hold on
yline(nstagesbed(1),'w')
yline(nstagesbed(end),'w')
hold off
title('Continuity ')

figure(2)
plotyy(Time, xx_0(end,:), Time, 1e3 * (sum(xx_0(0*nstages+1:1*nstages,ind) .* V_fluid) + sum(xx_0(1*nstages+1:2*nstages,ind) .* V_solid)) + xx_0(4*nstages+1,ind))

%%
xx_0(end,end)

%{
for i=1:numel(DATA)
    for d=1:numel(Di_Set)

        clc; [i,d]

        for kp=1:numel(kp_Set)

% load dataset
     LabResults = xlsread(DATA{i});

     V_Flow     = 0.39;

     T0homog    = LabResults(1,1)+273.15;
     feedPress  = LabResults(1,2) ;
     %rho        = LabResults(1,3) ;
     rho        = rhoPB_Comp(T0homog, feedPress, Compressibility(T0homog,feedPress,table2cell(Parameters_table(:,3))), table2cell(Parameters_table(:,3)));

     data_org   = LabResults(:,5)';
     data       = diff(data_org);

     % Set operating conditions
     feedTemp   = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
     feedTemp(round(numel(Time)/4):round(numel(Time)/2))   = feedTemp(1) + 50;

     feedPress  = feedPress * ones(1,length(Time_in_sec));  % Bars
     %feedPress(round(numel(Time)/3):round(2*numel(Time)/3))   = feedPress(1) + 50;
     
     feedFlow   = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec
     %feedFlow(round(numel(Time)/3):round(2*numel(Time)/3))   = 2*feedFlow(1) ;

     uu         = [feedTemp', feedPress', feedFlow'];

     % Initial conditions
     x0         = [C0fluid * ones(nstages,1);
                  C0solid * bed_mask;
                  T0homog*ones(nstages,1);
                  zeros(nstages,1);
                  0;
                  ];



            % ind = numel(Time); 1e3 * (sum(xx_0(0*nstages+1:1*nstages,ind) .* V_fluid) + sum(xx_0(1*nstages+1:2*nstages,ind) .* V_solid)) + xx_0(3*nstages+1,ind)

            
            %% Create the solver
            OPT_solver              = casadi.Opti();
            
            nlp_opts                = struct;
            nlp_opts.ipopt.max_iter = 5;
            ocp_opts                = {'nlp_opts', nlp_opts};
            OPT_solver.solver(         'ipopt'   , nlp_opts)

            % Descision variables
            k                       = OPT_solver.variable(Nk);
            k_lu                    = [ zeros(Nk,1) , inf*ones(Nk,1) ];
            % Constraints
            for nk=1:Nk
                OPT_solver.subject_to( k_lu(nk,1) <= k(nk,:) <= k_lu(nk,2) );
            end

            %% load parameters and set number of stages
            
            Parameters           = Parameters_table{:,3};               %
            Parameters_sym       = MX(Parameters_table{:,3});           % Vector of paraneters in the form casadi vector

            %% Assign new values of parameters to the Parameters vector
            
            %                       nstages, C0solid, V, epsi, dp, L, rho_s, km, mi
            Parameters(1:9)         = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];
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
            KOUT(:,kp,d,i) = kout;
            JJ(kp,d,i) = full(fJ(kout));
            J_STATUS{kp,d,i} = OPT_solver.return_status;

end
%save Fit_Di_km.mat
%clc
%disp('Saved')
end
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
%{
for v=1:numel(kp_Set)

    k0 = [kp_Set(v), 1];

    id = 1;
    for i = DATA_to_check
    
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
    
        kout  = KOUT(:,v,i);
    
        %
        Parameters(which_k) = k0;
        Parameters_opt = [uu repmat(Parameters,1,N_Time)'];
    
        [xx_0] = simulateSystem(F, [], x0, Parameters_opt  );
    
        Parameters(which_k) =kout;
        Parameters_opt = [uu repmat(Parameters,1,N_Time)'];
    
        [xx_out] = simulateSystem(F, [], x0, Parameters_opt  );
    
        %
        figure(v)
        set(gcf,'PaperOrientation','landscape')

        subplot(numel(DATA_to_check),3,id)
        imagesc(Time,1:nstages,xx_out(1:nstages,:)); colorbar
        xlabel('Time [min]')
        ylabel('Stages')
        title('Fluid Concentration')
        ylabel(['T=',mat2str(T0homog),', P=',mat2str(feedPress(1))])
        axis tight
    
        subplot(numel(DATA_to_check),3,id+1)
        imagesc(Time,1:nstages,xx_out(1*nstages+1:2*nstages,:)); colorbar
        xlabel('Time [min]')
        ylabel('Stages')
        title('Solid Concentration')
        axis tight
    
        subplot(numel(DATA_to_check),3,id+2)
    
        hold on
        plot(Time,xx_out(end,:))
        plot(Time,xx_0(end,:),'--')
        plot(SAMPLE,data,'o','MarkerSize',5)
        hold off
    
        legend('Estiamted','Inital','Data')
        legend off
        xlabel('Time [min]')
        ylabel('Mass of the extract [g]')
        title(['k=',sprintf('%0.3g, %0.3g, %0.3g, Cost Function=%0.3g, \n STATUS=%s',kout, JJ(v,i), J_STATUS{v,i} )  ])
        axis tight
        ylim([0, m_ref+2])

        print(gcf, ['Bounded_Fit_V_',mat2str(V),'_Kp_k0_',mat2str(k0),'.pdf'], '-dpdf','-r700','-fillpage')
    
        id = id + 3;
    
    end

end
%%
save Bounded_Fit_Kp.mat
%%
%}
%}