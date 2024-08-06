startup;
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

%% Create the solver
Iteration_max               = 1000;                                         % Maximum number of iterations for optimzer
Time_max                    = 96;                                           % Maximum time of optimization in [h]

nlp_opts                    = struct;
nlp_opts.ipopt.max_iter     = Iteration_max;
nlp_opts.ipopt.max_cpu_time = Time_max*3600;
%nlp_opts.ipopt.hessian_approximation ='limited-memory';
%nlp_opts.ipopt.acceptable_tol  = 1e-4;
%nlp_opts.ipopt.acceptable_iter = 5;

%%
CASE = 'T';

Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

LabResults              = xlsread('wpd_datasets.xlsx');
which_dataset           = 4;

%SAMPLE                  = LabResults(21:34,1);
data_org                = LabResults(21:34,which_dataset+1)';

%% Load paramters
m_total                 = 3.0;

% Bed geometry
before                  = 0.04;                                             % Precentage of length before which is empty
bed                     = 0.92;                                             % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 300;
timeStep                = 2.5;                                               % Minutes
OP_change_Time          = 20;                                              % Minutes
SAMPLING_TIME           = 5;

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);
OP_change               = OP_change_Time:OP_change_Time:ExtractionTime;

% Check if the number of data points is the same for both the dataset and the simulation
%%{
SAMPLE = SAMPLING_TIME:SAMPLING_TIME:ExtractionTime;
N_Sample                = [];
for i = 1:numel(SAMPLE)
    N_Sample            = [N_Sample ; find(round(Time,3) == round(SAMPLE(i))) ];
end
if numel(N_Sample) ~= numel(SAMPLE)
    keyboard
end
%}

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
Nx                      = 3 * nstages+2;                                    % 3*Nstages(C_f, C_s, H) + P(t) + yield
Nu                      = 3 + numel( Parameters );                          % T_in, P, F + numel(Parameters)

%% Extractor geometry
r                       = Parameters{3};                                    % Radius of the extractor  [m]
epsi                    = Parameters{4};                                    % Fullness [-]
dp                      = Parameters{5};                                    % Paritcle diameter
L                       = Parameters{6};                                    % Total length of the extractor [m]

L_nstages               = linspace(0,L,nstages);
V                       = L  * pi * r^2;                                    % Total volume of the extractor [m3]
A                       = pi *      r^2;                                    % Extractor cross-section

%--------------------------------------------------------------------
V_slice                 = (L/nstages) * pi * r^2;

V_before                = V_slice * numel(nstagesbefore);
V_after                 = V_slice * numel(nstagesafter);
V_bed                   = V_slice * numel(nstagesbed);                      % Volume of the fixed bed [m3]

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

% Integrator
F                       = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

%% Set inital state and inital conditions
msol_max                = m_total;                                          % g of product in solid and fluid phase
mSol_ratio              = 1;

mSOL_s                  = msol_max*mSol_ratio;                              % g of product in biomass
mSOL_f                  = msol_max*(1-mSol_ratio);                          % g of biomass in fluid

C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;                % Solid phase kg / m^3
Parameters{2}           = C0solid;

G                       =@(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;

m_fluid                 = G(L_bed_after_nstages)*( L_bed_after_nstages(2) ); % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
m_fluid                 = [zeros(1,numel(nstagesbefore)) m_fluid];
C0fluid                 = m_fluid * 1e-3 ./ V_fluid';

num_of_exp = 6;
ppm = ParforProgressbar(num_of_exp);

KOUT_MATRIX = nan(num_of_exp, 2*numel(OP_change));   % decision var + inital obj + final obj
OBJ_MATRIX  = nan(num_of_exp, 2);

parfor ii = 1:num_of_exp
    
    %% Set operating conditions
    T0homog                 = 35+273;                                           % K
    feedPress               = 150;                                              % MPa -> bar
    Flow                    = 5;                                                % kg/s
    
    Z                       = Compressibility( T0homog, feedPress,         Parameters );
    rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );
    
    enthalpy_rho            = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );
    
    feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;     % Kelvin
    
    feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;     % Bars
    
    feedFlow                = Flow * ones(1,length(Time_in_sec)) * 1e-5;        % kg/s
    
    uu                      = [feedTemp', feedPress', feedFlow'];
    
    MU                      =  Viscosity(T0homog,rho);
    VELOCITY                =  Velocity(Flow, rho, Parameters);
    RE                      =  dp .* rho .* VELOCITY ./ MU;
    
    % Initial conditions
    x0                      = [ C0fluid'                         ;
                                C0solid         * bed_mask       ;
                                enthalpy_rho    * ones(nstages,1);
                                feedPress(1)                     ;
                                0                                ;
                                ];
    
    %% Set the inital simulation and plot it against the corresponding dataset
    %Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];
    %[xx_0]                 = simulateSystem(F, [], x0, Parameters_init_time );
    
    %% Solve the optimization problem with inital guess
    %hold on
    %plot(Time,xx_0(end,:))
    %plot(SAMPLE,data_org,'o')
    %hold off

    OPT_solver                  = casadi.Opti();
    ocp_opts                    = {'nlp_opts', nlp_opts};
    OPT_solver.solver(             'ipopt'   , nlp_opts)

    %% Define decision variables

    switch CASE 
        case 'T'
            T0homog         = OPT_solver.variable(numel(OP_change))';
                              OPT_solver.subject_to( 30+273 <= T0homog <= 40+273 );

            feedTemp        = repmat(T0homog,OP_change_Time/timeStep,1);
            feedTemp        = feedTemp(:)';
            feedTemp        = [ feedTemp, T0homog(end)*ones(1,N_Time - numel(feedTemp)) ];    
        
        case 'P'
            Press           = OPT_solver.variable(numel(OP_change))';
                              OPT_solver.subject_to( 100 <= Press <= 200 );

            feedPress       = repmat(Press,OP_change_Time/timeStep,1);
            feedPress       = feedPress(:)';
            feedPress       = [ feedPress, Press(end)*ones(1,N_Time - numel(feedPress)) ];    
    end
    
    Flow                    = OPT_solver.variable(numel(OP_change))';
                              OPT_solver.subject_to( 3.33 <= Flow <= 6.67 )
    
    T_0                     = feedTemp(1);   
    
    Z                       = Compressibility( T_0, feedPress(1),         Parameters );
    rho                     = rhoPB_Comp(      T_0, feedPress(1), Z,      Parameters );
    enthalpy_rho            = rho.*SpecificEnthalpy(T_0, feedPress(1), Z, rho, Parameters ) ;
    
    %feedFlow                = Flow * ones(1,length(Time_in_sec));               % kg/s
    feedFlow                = repmat(Flow,OP_change_Time/timeStep,1) * 1e-5;
    feedFlow                = feedFlow(:)';
    feedFlow                = [ feedFlow, Flow(end)*ones(1,N_Time - numel(feedFlow)) ];    
    
    uu                      = [feedTemp', feedPress', feedFlow'];
    
    %% Initial conditions
    x0                      = [ C0fluid'                         ;
                                C0solid         * bed_mask       ;
                                enthalpy_rho    * ones(nstages,1);
                                feedPress(1)                     ;
                                0                                ;
                                ];
    
    sigma                   = 0.12;
    which_theta             = [44:46];
    
    % Store symbolic results of the simulation
    Parameters_sym          = MX(cell2mat(Parameters));
    Parameters_sym_t        = OPT_solver.parameter(numel(which_theta));
    Parameters_sym(which_theta)   = Parameters_sym_t;
    
    X                       = MX(Nx,N_Time+1);
    X(:,1)                  = x0;
    
    %% Symbolic integration
    for jj=1:N_Time
        X(:,jj+1)=F(X(:,jj), [uu(jj,:)'; Parameters_sym] );
    end
    
    %% Obtain symbolic yield
    Yield_estimate          = X(Nx,[1; N_Sample]);
    data_obs                = diff(Yield_estimate);
    
    %% Finad jacobian abd FI matrix
    S                       = jacobian(data_obs, Parameters_sym_t);
    
    S_norm                  = S' .* abs(cell2mat(Parameters(which_theta))) ./ repmat(data_obs,numel(which_theta) );
    FI                      = (S_norm * S_norm');
    FI                      = FI ./ (sigma^2);  
    
    if numel(FI) ~= numel(which_theta)^2
        keyboard
    end
    
    D_opt = -log(myDet(FI));
    
    %% Defin intial guesses
    switch CASE
        case 'T'
            %T0 = linspace(30,40,numel(T0homog))+273;
            T0 = ( (40-30).*rand(1,numel(T0homog)) + 30 ) + 273;
        case 'P'
            P0 = ( (200 -100) .*rand(1,numel(Press)) + 100  ) ;
    end
    %F0 = linspace(3.33,6.67,numel(Flow));
    F0 = ( (6.67-3.33).*rand(1,numel(Flow)) + 3.33 ) ;
    
    %% Solve the optimization problem
    OPT_solver.set_value(Parameters_sym_t, cell2mat(Parameters(which_theta)));
    OPT_solver.minimize(D_opt);
    switch CASE
        case 'T'
            OPT_solver.set_initial([T0homog, Flow], [T0, F0] );
        case 'P'
            OPT_solver.set_initial([Press, Flow], [P0, F0] );
    end
    
    %%
    try
        sol  = OPT_solver.solve();
        switch CASE
            case 'T'
                KOUT = full(sol.value([T0homog, Flow])) 
            case 'P'
                KOUT = full(sol.value([Press, Flow])) 
        end
    catch
        switch CASE
            case 'T'
                KOUT = OPT_solver.debug.value([T0homog, Flow])
            case 'P'
                KOUT = OPT_solver.debug.value([Press, Flow])
        end
    end

    OBJ = OPT_solver.stats.iterations.obj;
    
    OBJ_MATRIX(ii,:)  = OBJ([1,end]);
    KOUT_MATRIX(ii,:) = KOUT;

    pause(100/num_of_exp);
    % increment counter to track progress
    ppm.increment();

end

delete(ppm);

%% Save data
%{
switch CASE_o
    case 'T'
        save DOE_T.mat
    case 'P'
        save DOE_P.mat
end
%}
%% scatter plot to sumerize the optimization results
AA = sort(OBJ_MATRIX(:,2));

figure()
scatterhist(OBJ_MATRIX(:,2),OBJ_MATRIX(:,1), 'Kernel', 'on', 'MarkerSize',4, 'Color','k')
xlabel('Final value of the objective function')
ylabel('Inital value of the objective function')
hold on
xline(AA(1))
xline(AA(3))
patch([AA(1) AA(3) AA(3) AA(1)], [max(OBJ_MATRIX(:,1)) max(OBJ_MATRIX(:,1)) min(OBJ_MATRIX(:,1)) min(OBJ_MATRIX(:,1))],'red','FaceAlpha',0.25, 'EdgeColor','none')
hold off
axis tight

%exportgraphics(figure(1), ['scatter_P.png'], "Resolution",300);
%close all

%% plot 3 optimal profile
for i = 1:3
    figure()
    II = find(OBJ_MATRIX(:,2) == AA(i));

    yyaxis right
    stairs([0 OP_change], [KOUT_MATRIX(II,1:numel(OP_change)) KOUT_MATRIX(II,numel(OP_change))], 'LineWidth', 2)
    ylabel('Pressure bar')

    yyaxis left
    stairs([0 OP_change], [KOUT_MATRIX(II,numel(OP_change)+1:end) KOUT_MATRIX(II,end)], 'LineWidth', 2)
    ylabel('Flow rate kg/s $\cdot 10^{-5}$')

    xlabel('Time min')

    %exportgraphics(figure(1), ['Profile_P_',num2str(i),'.png'], "Resolution",300);
    %close all
end

%%
%{
%% Set operating conditions
T0homog                 = T0;
Flow                    = F0;

feedTemp                = repmat(T0homog,OP_change_Time/timeStep,1);
feedTemp                = feedTemp(:)';
%feedTemp                = [ feedTemp, T0homog(end)*ones(1,N_Time - numel(feedTemp)) ];    

T_0                     = feedTemp(1);   

Z                       = Compressibility( T_0, feedPress(1),         Parameters );
rho                     = rhoPB_Comp(      T_0, feedPress(1), Z,      Parameters );
enthalpy_rho            = rho.*SpecificEnthalpy(T_0, feedPress(1), Z, rho, Parameters ) ;

%feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;     % Bars

feedFlow                = repmat(Flow,OP_change_Time/timeStep,1);
feedFlow                = feedFlow(:)';
%feedFlow                = [ feedFlow, Flow(end)*ones(1,N_Time - numel(feedFlow)) ];  
feedFlow                = feedFlow * 1e-5;        % kg/s

uu                      = [feedTemp', feedPress', feedFlow'];

% Initial conditions
x0                      = [ C0fluid'                         ;
                            C0solid         * bed_mask       ;
                            enthalpy_rho    * ones(nstages,1);
                            feedPress(1)                     ;
                            0                                ;
                            ];

%% Set the inital simulation and plot it against the corresponding dataset
Parameters_opt_0       = [uu repmat(cell2mat(Parameters),1,N_Time)'];
[xoxo_0]               = simulateSystem(F, [], x0, Parameters_opt_0 );

%% Set operating conditions
T0homog                 = KOUT(1:numel(T0));
Flow                    = KOUT(numel(F0)+1:end);

feedTemp                = repmat(T0homog,OP_change_Time/timeStep,1);
feedTemp                = feedTemp(:)';
%feedTemp                = [ feedTemp, T0homog(end)*ones(1,N_Time - numel(feedTemp)) ];    

T_0                     = feedTemp(1);   

Z                       = Compressibility( T_0, feedPress(1),         Parameters );
rho                     = rhoPB_Comp(      T_0, feedPress(1), Z,      Parameters );
enthalpy_rho            = rho.*SpecificEnthalpy(T_0, feedPress(1), Z, rho, Parameters ) ;

%feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;     % Bars

feedFlow                = repmat(Flow,OP_change_Time/timeStep,1);
feedFlow                = feedFlow(:)';
%feedFlow                = [ feedFlow, Flow(end)*ones(1,N_Time - numel(feedFlow)) ];  
feedFlow                = feedFlow * 1e-5;        % kg/s

uu                      = [feedTemp', feedPress', feedFlow'];

% Initial conditions
x0                      = [ C0fluid'                         ;
                            C0solid         * bed_mask       ;
                            enthalpy_rho    * ones(nstages,1);
                            feedPress(1)                     ;
                            0                                ;
                            ];

%% Set the inital simulation and plot it against the corresponding dataset
Parameters_opt         = [uu repmat(cell2mat(Parameters),1,N_Time)'];
[xoxo]                 = simulateSystem(F, [], x0, Parameters_opt );

%% Plot Yield curves for inital and opt system and controls

hold on
plot(Time, full(xoxo(end,:)), 'LineWidth', 2)
plot(Time, full(xoxo_0(end,:)), 'LineWidth', 2)
hold off
legend('Optimized solution','Inital guess', 'Location','Best', 'box', 'off')
ylabel('Yield g')
xlabel('Time min')
%set(gcf,'PaperOrientation','landscape'); print(figure(1),['1.pdf'],'-dpdf','-bestfit')
exportgraphics(figure(1), ['1.png'], "Resolution",300);
close all

hold on
stairs([0 SAMPLE],[KOUT(1:30) KOUT(30)], 'LineWidth', 2)
stairs([0 SAMPLE],[T0 T0(end)], 'LineWidth', 2)
hold off
ylabel('Temp K')
xlabel('Time min')
%legend('Optimized solution','Inital guess', 'Location','Best', 'box', 'off')
%set(gcf,'PaperOrientation','landscape'); print(figure(1),['2.pdf'],'-dpdf','-bestfit')
exportgraphics(figure(1), ['2.png'], "Resolution",300);
close all

hold on
stairs([0 SAMPLE],[KOUT(31:end) KOUT(end)]* 1e-5, 'LineWidth', 2)
stairs([0 SAMPLE],[F0 F0(end)]* 1e-5, 'LineWidth', 2)
hold off
ylabel('Mass flow rate kg/s')
xlabel('Time min')
%legend('Optimized solution','Inital guess', 'Location','Best', 'box', 'off')
%set(gcf,'PaperOrientation','landscape'); print(figure(1),['3.pdf'],'-dpdf','-bestfit')
exportgraphics(figure(1), ['3.png'], "Resolution",300);
close all

%}

%}