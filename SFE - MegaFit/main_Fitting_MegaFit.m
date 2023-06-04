startup; datetime("now")
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

Parameters_table        = readtable('Parameters.csv') ;        % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});
Parameters_sym          = MX(cell2mat(Parameters));

%% Create the solver
OPT_solver                  = casadi.Opti();

Iteration_max               = 100;                                          % Maximum number of iterations for optimzer
Time_max                    = 12.0;                                         % Maximum time of optimization in [h]

nlp_opts                    = struct;
nlp_opts.ipopt.max_iter     = Iteration_max;
nlp_opts.ipopt.max_cpu_time = Time_max*3600;     

ocp_opts                    = {'nlp_opts', nlp_opts};
OPT_solver.solver(             'ipopt'   , nlp_opts)

which_k                 = [     8, 45,          ];                      % Select which parameters are used for fitting
k0                      = [     1,  1, 0.8,  0.4];                      % Give inital value for each parameter
k_lu                    = [ [0.05;  0;   0;    0], ...
                          [  100; inf;   1; inf]];
Nk                      = numel(which_k)+2;                                 % Parameters within the model + (m_max), m_ratio, sigma

m_total                 = 90; 

%ratio                   = OPT_solver.variable(2);
                          %OPT_solver.subject_to( 0 <= sat <= 1 );
%Parameters_sym(44)      = ratio(1);
%Parameters_sym(47)      = ratio(2);
                          %OPT_solver.subject_to( 0 <= DI_correlation(800,ratio) );
                          %OPT_solver.subject_to( 0 <= DI_correlation(950,ratio) );

%% Load paramters
DATA_set                = {'LUKE_T40_P200', 'LUKE_T50_P200', 'LUKE_T40_P300_org', 'LUKE_T50_P300_org'};

V_Flow                  = 0.4;                                              % Volumetric flow rate l/min

% Bed geometry
before                  = 0.1;                                              % Precentage of length before which is empty
bed                     = 0.165;                                            % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 200;
ExtractionTime_SA       = 300;
timeStep                = 1;                                                % Minutes
SamplingTime            = 5;                                                % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);

SAMPLE                  = [PreparationTime:SamplingTime:150];

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
nstages                 = Parameters{1};

nstagesbefore           = 1:floor(before*nstages);
nstagesbed              = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
nstagesafter            = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

%% Number of variables
Nx                      = 3 * nstages+2;                                % 3*Nstages(C_f, C_s, H) + P(t) + yield
Nu                      = 3 + numel( Parameters );

%% Extractor geometry
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

L_bed_after_nstages     = L_nstages(nstagesbed(1):end);
L_bed_after_nstages     = L_bed_after_nstages - L_bed_after_nstages(1);
L_end                   = L_bed_after_nstages(end);

%%
Courant_Number = nan(1,numel(DATA_set));

%parpool(4)

cost = 0;
K    = [];
RHO  = [];

for ii=1:numel(DATA_set)

    
    %% symbolic variables
    x                           = MX.sym('x', Nx);
    u                           = MX.sym('u', Nu);
    
    %% Set Integrator
    f                           = @(x, u) modelSFE(x, u, bed_mask, timeStep_in_sec);
    xdot                        = modelSFE(x, u, bed_mask, timeStep_in_sec);
    
    % Integrator
    F                           = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

    % Descision variables
    k                           = OPT_solver.variable(Nk);

    % Dataset
    DATA                        = DATA_set{ii};
    LabResults                  = xlsread([DATA,'.xlsx']);

    %% Set parameters
    msol_max                    = m_total;%k(numel(which_k)+1);                                              % g of product in solid and fluid phase
    mSol_ratio                  = k(numel(which_k)+1);
    sigma                       = k(numel(which_k)+2);

    mSOL_s                      = msol_max*mSol_ratio;                                              % g of product in biomass
    mSOL_f                      = msol_max*(1-mSol_ratio);                                          % g of biomass in fluid

    %%
    Parameters_sym(which_k)     = k(1:numel(which_k));    

    %--------------------------------------------------------------------

    C0solid                     = mSOL_s * 1e-3 / ( V_bed * epsi)  ;                                 % Solid phase kg / m^3
    Parameters_sym(2)           = C0solid;
    
    G                           =@(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;
    
    m_fluid                     = G(L_bed_after_nstages)*( L_bed_after_nstages(2) );                 % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
    m_fluid                     = [zeros(1,numel(nstagesbefore)) m_fluid];         
    C0fluid                     = m_fluid * 1e-3 ./ V_fluid';

    %%
    
    T0homog                     = LabResults(1,1)+273.15;
    feedPress                   = LabResults(1,2);

    Z                           = Compressibility( T0homog, feedPress,         Parameters );
    rho                         = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );
    RHO                         = [RHO; rho];

    enthalpy_rho                = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );

    % Set operating conditions
    feedTemp                    = T0homog   * ones(1,length(Time_in_sec)) + 0 ;  % Kelvin

    feedPress                   = feedPress * ones(1,length(Time_in_sec)) + 0 ;  % Bars

    feedFlow                    = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/s

    Courant_Number(ii)          = Velocity(feedFlow(1), rho, Parameters) * timeStep_in_sec / (L_nstages(2));

    uu                          = [feedTemp', feedPress', feedFlow'];

    % Initial conditions
    x0                          = [ C0fluid'                        ;
                                    C0solid      * bed_mask         ;
                                    enthalpy_rho * ones(nstages,1)  ;
                                    feedPress(1)                    ;
                                    0                               ;
                                    ];

    %%
    % Store symbolic results of the simulation
    X                           = MX(Nx,N_Time+1);
    X(:,1)                      = x0;

    % Symbolic integration
    for j=1:N_Time
        X(:,j+1)=F(X(:,j), [uu(j,:)'; Parameters_sym] );
    end

    %% Find the measurment from the simulation
    Yield_estimate              = X(Nx,N_Sample);
    Yield_estimate_diff         = diff(Yield_estimate);

    %% load dataset
    data_org                    = LabResults(:,5)';
    data                        = diff(data_org);

    %% Create the cost function
    %J                         = (data-Yield_estimate_diff ) * diag([1:numel(data)/2, numel(data)/2:-1:1]) * (data-Yield_estimate_diff )';
    J                           = (data-Yield_estimate_diff ) * diag(1) * (data-Yield_estimate_diff )';
    J_L                         = -numel(data)./2 .* ( log(2.*pi) + log(sigma^2) ) - J./(2.*sigma^2);
    J_L                         = - J_L;
    %J_L     = 1./sqrt( ( 2*3.14.*sigma.^2).^(numel(data)) ) .* exp( -J./(2.*sigma.^2) )

    %fJ  = Function('fJ', {k}, {J} );

    cost = cost + J_L;
    K = [K; k];
    
    %% Constraints
    for nk=1:Nk
        OPT_solver.subject_to( k_lu(nk,1) <= k(nk,:) <= k_lu(nk,2) );
    end

    OPT_solver.subject_to( Yield_estimate(end)+1 <= X(Nx,end) );

end

%% Set opt and inital guess
OPT_solver.minimize(cost);
% 
K0 = repmat(k0,1,numel(DATA_set));
%K0 = [0.02141,-14.5, K0];
%K  = [ratio;K];

OPT_solver.set_initial(K, K0);
% 
% %% Solve the opt
tic
% 
try
        sol = OPT_solver.solve();
        KOUT = full(sol.value(K));
catch
        KOUT = full(OPT_solver.debug.value(K));
end
datetime("now")

%m_total    = KOUT(1:2);
%KOUT       = KOUT(3:end);
DATA_K_OUT = reshape(KOUT,Nk,[]);

%Parameters{44} = ratio(1);
%Parameters{47} = ratio(2);

%% save parameter estimation results in csv format

%format shortE

NAME = {'$k_m[-]$', '$D_i[m^2/s] \times 1e-14$', '$D_e^M[m^2/s] \times 1e-6$', 'sat', '$\tau[-]$', '$\sigma[-]$'};
name = {'km', 'Di', 'De', 'sat', 'tau', 'sigma'};
Title = {'T40P200', 'T50P200', 'T40P300', 'T50P300'};

TT   = cell2table(NAME');

DATA_K_OUT_order_of_mag = DATA_K_OUT;
DATA_K_OUT_order_of_mag(2,:) = DATA_K_OUT_order_of_mag(2,:);
DATA_K_OUT_order_of_mag(3,:) = DATA_K_OUT_order_of_mag(3,:);
%DATA_K_OUT_order_of_mag(4,:) = DATA_K_OUT_order_of_mag(4,:) * 1e-6;

for i=1:numel(DATA_set)
    TT.(i+1)=round(DATA_K_OUT_order_of_mag(:,i),5);
end

TT.Properties.VariableNames = ["Parameter","$40[C] 200[bar]$","$50[C] 200[bar]$","$40[C] 300[bar]$","$50[C] 300[bar]$"];
writetable(TT,'estimation.csv')
TT

%format short

%% save dataset in csv format

TTT = table([0:5:150]');

for ii=1:numel(DATA_set)
    DATA                    = DATA_set{ii};
    LabResults              = xlsread([DATA,'.xlsx']);
    data_org                = LabResults(:,5);
    TTT.(ii+1)=data_org;
end

TTT.Properties.VariableNames = ["Time [min]","$40[C] 200[bar]$","$50[C] 200[bar]$","$40[C] 300[bar]$","$50[C] 300[bar]$"];
writetable(TTT,'Dataset.csv')
TTT

%%  Inrease the simulation for sensitivity analysis
% simulationTime          = PreparationTime + 600;
% 
% timeStep_in_sec         = timeStep * 60;                                    % Seconds
% Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
% Time                    = [0 Time_in_sec/60];                               % Minutes
% 
% N_Time                  = length(Time_in_sec);

name_v                  = 'F';

RHO   = [];
U     = [];
for ii=1:numel(DATA_set)


    DATA                    = DATA_set{ii};
    LabResults              = xlsread([DATA,'.xlsx']);

    data_org                = LabResults(:,5)';
    data                    = diff(data_org);

    KOUT                    = DATA_K_OUT(:,ii);
    
    % Set operating conditions
    T0homog                 = LabResults(1,1)+273.15;
    feedPress               = LabResults(1,2);

    Z                       = Compressibility( T0homog, feedPress,         Parameters );
    rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );
    RHO                     = [RHO, rho];

    enthalpy_rho            = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );

    
    feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;  % Kelvin

    feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;  % Bars

    feedFlow                = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/s

    uu                      = [feedTemp', feedPress', feedFlow'];

    % Solve for the inital guess
    for i=1:numel(which_k)
        Parameters{which_k(i)}  = k0(i);
    end

    msol_max                = m_total;%k0(numel(which_k)+1);                                     % g of product in solid and fluid phase
    mSol_ratio              = k0(numel(which_k)+1);

    mSOL_s                  = msol_max*mSol_ratio;                                               % g of product in biomass
    mSOL_f                  = msol_max*(1-mSol_ratio);                                           % g of biomass in fluid

    C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;                                 % Solid phase kg / m^3
    Parameters{2}           = C0solid;

    G                       =@(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;

    m_fluid                 = G(L_bed_after_nstages)*( L_bed_after_nstages(2) );                 % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
    m_fluid                 = [zeros(1,numel(nstagesbefore)) m_fluid];         
    C0fluid                 = m_fluid * 1e-3 ./ V_fluid';

    % Initial conditions
    x0                      = [ C0fluid'                        ;
                                C0solid  * bed_mask             ;
                                enthalpy_rho * ones(nstages,1)  ;
                                feedPress(1)                    ;
                                0                               ;
                                ];

    Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];
    [xx_0]                 = simulateSystem(F, [], x0, Parameters_init_time );
 

    % Solve for the optimzed parameters
    Parameters_opt = Parameters;
    for i=1:numel(which_k)
        Parameters_opt{which_k(i)}  = KOUT(i);
    end

    % 
    msol_max_opt            = m_total;%KOUT(numel(which_k)+1);                                    % g of product in solid and fluid phase
    mSol_ratio_opt          = KOUT(numel(which_k)+1);
    mSOL_s_opt              = msol_max_opt*mSol_ratio_opt;                                        % g of product in biomass
    mSOL_f_opt              = msol_max_opt*(1-mSol_ratio_opt);                                    % g of biomass in fluid

    C0solid_opt             = mSOL_s_opt * 1e-3 / ( V_bed * epsi)  ;                              % Solid phase kg / m^3
    Parameters_opt{2}       = C0solid_opt;

    G                       =@(x) -(2*mSOL_f_opt / L_end^2) * (x-L_end) ;
    
    m_fluid_opt             = G(L_bed_after_nstages)*( L_bed_after_nstages(2) );                  % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
    m_fluid_opt             = [zeros(1,numel(nstagesbefore)) m_fluid_opt];         
    C0fluid_opt             = m_fluid_opt * 1e-3 ./ V_fluid';

    % Initial conditions
    x0_opt                  = [ C0fluid_opt'                   ;
                                C0solid_opt  * bed_mask        ;
                                enthalpy_rho * ones(nstages,1) ;
                                feedPress(1);
                                0;
                                ];

    Parameters_opt_time     = [uu repmat(cell2mat(Parameters_opt),1,N_Time)'];
    [xx_out]                = simulateSystem(F, [], x0_opt, Parameters_opt_time);

    Plot_Fit

    %

    [S,p,Sdot]              = Sensitivity(x, xdot, u, [3] );

    % Initial conditions
        x0_SA                   = [ C0fluid';
                                    C0solid  * bed_mask        ;
                                    enthalpy_rho * ones(nstages,1) ;
                                    feedPress(1);
                                    0;
                                    zeros(length(S)-length(xdot),1)
                                    ];

        f_SA = @(S, p) Sdot(S, p, bed_mask);
        Results = Integrator_SS(Time*60, x0_SA, S, p, Sdot, Parameters_opt_time );
        Res = Results(Nx+1:end,:) ;

        subplot(3,2,1)
        imagesc(Time, L_nstages, Res(1*nstages+1:2*nstages,:)); cb = colorbar;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'w--',{'end of bed','beginning of bed'}, 'Interpreter', 'latex')
        hold off
        cb.Label.String = ['$\frac{d c_s}{d',name_v,'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 14;
        xlabel('Time [min]'); ylabel('Length [m]'); title(sprintf('$\\rho_f$=%4.2f',rho))

        subplot(3,2,2)
        imagesc(Time, L_nstages, Res(2*nstages+1:3*nstages,:) ); cb = colorbar;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'w--',{'end of bed','beginning of bed'}, 'Interpreter', 'latex')
        hold off
        cb.Label.String = ['$\frac{d h}{d',name_v,'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 14;
        xlabel('Time [min]'); ylabel('Length [m]')

        subplot(3,2,3)
        imagesc(Time, L_nstages, Res(0*nstages+1:1*nstages,:)); cb = colorbar;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'w--',{'end of bed','beginning of bed'}, 'Interpreter', 'latex')
        hold off
        cb.Label.String = ['$\frac{d c_f}{d',name_v,'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 14;
        xlabel('Time [min]'); ylabel('Length [m]')

        subplot(3,2,4)
        plot(Time, Res(end-1,:))
        xlabel('Time [min]'); ylabel(['$\frac{d P_{t-1}}{d',name{ii},'}$'])

        subplot(3,2,5)
        hold on
        plot(Time, Res(end,:)); colorbar;
        xlabel('Time [min]'); ylabel(['$\frac{d y}{d',name_v,'}$'])
        hold off
        colormap jet

        set(gcf,'PaperOrientation','landscape'); print(figure(1),['DI_',Title{ii},'.pdf'],'-dpdf','-bestfit')
        close all; 

        U = [U; Res(end,:)];


end

plot(Time, U);
xlabel('Time [min]'); ylabel(['$\frac{d y}{d',name_v,'}$'])
legend(string(RHO))
set(gcf,'PaperOrientation','landscape'); print(figure(1),[name_v,'_all.pdf'],'-dpdf','-bestfit'); close all

%%

Plot_Trend_Line

%%

T = [40 50 40 50];
P = [200 200 300 300];
for ii=1:Nk
    km = DATA_K_OUT(ii,:);
    createSurfFit(T, P, km, NAME{ii}, name{ii});
end
