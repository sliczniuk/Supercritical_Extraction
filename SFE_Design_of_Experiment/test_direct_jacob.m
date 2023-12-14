startup;
delete(gcp('nocreate'));
%p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

excel_file = 'output_P250.xls';
numeEval = 69;

parfor ii=61:numeEval

Parameters_table                = readtable('Parameters.csv') ;             % Table with prameters
Parameters                      = num2cell(Parameters_table{:,3});          % Parameters within the model + (m_max), m_ratio, sigma

%% Create the solver
Iteration_max                   = 100;                                       % Maximum number of iterations for optimzer
Time_max                        = 12.0;                                      % Maximum time of optimization in [h]

nlp_opts                        = struct;
nlp_opts.ipopt.max_iter         = Iteration_max;
nlp_opts.ipopt.max_cpu_time     = Time_max*3600;
nlp_opts.ipopt.hessian_approximation ='limited-memory';
nlp_opts.ipopt.tol              = 1e-6;
nlp_opts.ipopt.acceptable_tol   = 1e-4;
%nlp_opts.ipopt.acceptable_iter = 5;

OPT_solver                      = casadi.Opti();
ocp_opts                        = {'nlp_opts', nlp_opts};
OPT_solver.solver(                  'ipopt'   , nlp_opts)

%% Load paramters
m_total                 = 80;

V_Flow                  = 0.4;                                              % Volumetric flow rate l/min

% Bed geometry
before                  = 0.1;                                              % Precentage of length before which is empty
bed                     = 0.165;                                            % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 150;
DeadTime                = 25;
timeStep                = 1;                                                % Minutes
timeShited              = 0;                                                % Minutes
SamplingTime            = 5;                                                % Minutes
OP_change_Time          = 5;                                                % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);
SAMPLE                  = SamplingTime:SamplingTime:ExtractionTime;
OP_change               = timeShited:OP_change_Time:(ExtractionTime-DeadTime);

% Check if the number of data points is the same for both the dataset and the simulation
N_Sample                = [];
for i = 1:numel(SAMPLE)
    N_Sample            = [N_Sample ; find(round(Time,3) == round(SAMPLE(i))) ];
end
if numel(N_Sample) ~= numel(SAMPLE)
    keyboard
end

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

%% Set operating conditions
% bounds for T and F
T_max                   = 50+273;
T_min                   = 40+273;

% set P
feedPress               = 250;
feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;     % Bars

% set T
T0homog                 = OPT_solver.variable(numel(OP_change))';
                          OPT_solver.subject_to( T_min <= T0homog <= T_max );

T_0                     = T0homog(1);   
feedTemp                = repmat(T0homog,OP_change_Time/timeStep,1);
feedTemp                = feedTemp(:)';
feedTemp                = [ feedTemp, T0homog(end)*ones(1,N_Time - numel(feedTemp)) ];    

% set h
Z                       = Compressibility( T_0, feedPress(1),         Parameters );
rho                     = rhoPB_Comp(      T_0, feedPress(1), Z,      Parameters );

enthalpy_rho            = rho.*SpecificEnthalpy(T_0, feedPress(1), Z, rho, Parameters );

feedFlow                = V_Flow .* rho .* 1e-3 ./ 60 * ones(1,length(Time_in_sec)) ;  % l/min -> kg/s

% set vector of controls
uu                      = [feedTemp', feedPress', feedFlow'];

%% Set inital state and inital conditions
msol_max                = m_total;                                          % g of product in solid and fluid phase
mSol_ratio              = 0.7;

mSOL_s                  = msol_max*mSol_ratio;                              % g of product in biomass
mSOL_f                  = msol_max*(1-mSol_ratio);                          % g of biomass in fluid

C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;                % Solid phase kg / m^3
Parameters{2}           = C0solid;

G                       = @(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;

m_fluid                 = G(L_bed_after_nstages)*( L_bed_after_nstages(2) ); % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
m_fluid                 = [zeros(1,numel(nstagesbefore)) m_fluid];
C0fluid                 = m_fluid * 1e-3 ./ V_fluid';

% Initial conditions
x0                      = [ C0fluid'                        ;
                            C0solid  * bed_mask             ;
                            enthalpy_rho * ones(nstages,1)  ;
                            feedPress(1)                    ;
                            0                               ;
                            ];

%%
sigma                   = 0.12;
which_theta             = [44:47];

% Store symbolic results of the simulation
Parameters_sym          = MX(cell2mat(Parameters));
Parameters_sym_t        = OPT_solver.parameter(numel(which_theta));
Parameters_sym(which_theta)   = Parameters_sym_t;

X                       = MX(Nx,N_Time+1);
X(:,1)                  = x0;
% Symbolic integration
for jj=1:N_Time
    X(:,jj+1)=F(X(:,jj), [uu(jj,:)'; Parameters_sym] );
end

%% Find the measurment from the simulation
Yield_estimate          = X(Nx,[1; N_Sample]);
data_obs                = diff(Yield_estimate);

S                       = jacobian(data_obs, Parameters_sym_t);

S_norm                  = S' .* abs(cell2mat(Parameters(which_theta))) ./ repmat(data_obs,numel(which_theta) );
FI                      = (S_norm * S_norm');
FI                      = FI ./ (sigma^2);  

if numel(FI) ~= numel(which_theta)^2
    keyboard
end

D_opt = -log(myDet(FI));

OPT_solver.set_value(Parameters_sym_t, cell2mat(Parameters(which_theta)));

OPT_solver.minimize(D_opt);

T0 = (T_max-T_min).*rand(1,numel(T0homog)) + T_min;

OPT_solver.set_initial(T0homog, T0 );

try
    sol  = OPT_solver.solve();
    KOUT = full(sol.value([T0homog]))
catch
    KOUT = OPT_solver.debug.value([T0homog])
end

%% Reconstruct T profile
FF     = Function('FF', {[T0homog, Parameters_sym_t']}, {X});
XX     = full(FF( [KOUT, cell2mat(Parameters(which_theta))'] ));

TT_rec = full(Reconstruct_T_from_enthalpy(XX(3*nstages,2:end)',250, Parameters))';

FF     = Function('FF', {T0homog}, {feedTemp});
TT     = full(FF( [KOUT] ));
TT_0   = full(FF( [T0] ));

%% save data

OBJ = OPT_solver.stats.iterations.obj;
    
writematrix('T0',excel_file,'sheet',ii,'Range','A1');
writematrix('T_opt',excel_file,'sheet',ii,'Range','B1');
writematrix('T_out',excel_file,'sheet',ii,'Range','C1');
writematrix('OBJ',excel_file,'sheet',ii,'Range','D1');
writematrix(OPT_solver.return_status,excel_file,'sheet',ii,'Range','E1');

writematrix([TT_0; TT; TT_rec]',excel_file,'sheet',ii,'Range','A2:C200');
writematrix(OBJ',excel_file,'sheet',ii,'Range','D2');
%}
end

%% read data
OBJ_data = [];
OBJ_init = [];
Group    = {};

TT_opt   = [];
TT_out   = [];

figure(1)
%subplot(2,1,1)
%hold on
for ii=1:numeEval
    data = readcell(excel_file,'sheet',ii);
    OBJ_init = [OBJ_init, data{2,4}];
    OBJ_data = [OBJ_data, data{end,4}];

    TT_opt = [TT_opt, cell2mat(data(2:151,2))];
    TT_out = [TT_out, cell2mat(data(2:151,3))];

    Group = [Group; data{1,5}];

    %plot(cell2mat(data(2:end,3)));
end
%hold off

%subplot(2,1,2)
plot(TT_opt)

figure(2)
scatterhist(OBJ_init, OBJ_data, 'Kernel','on','Location','SouthEast', 'Direction','out');

figure(3)
h = scatterhist(OBJ_init, OBJ_data, 'Group', Group, 'Kernel','off','Location','SouthEast', 'Direction','out');
hold on
for ii = 1:numeEval
    text(OBJ_init(ii),OBJ_data(ii),num2str(ii))
end
yy = yline(min(OBJ_data),'-');
hold off
L = legend;
L.String(end) = {'Mini(obj)'};
legend boxoff 
%delete(h(2));
fontsize(16,"points")
set(gcf,'PaperOrientation','landscape'); print(figure(3),['Multiple_shot_DOE_P250.pdf'],'-dpdf','-bestfit')

T_in  = TT_opt(:,[16,33,45,57]);
T_out = TT_out(:,[16,33,45,57]);

figure(4)

plot([0:149],T_in, 'Linewidth', 2);

xlabel('Time[min]')
ylabel('Temperature [K]')

ylim([35, 55]+273)

fontsize(24,"points")
set(gcf,'PaperOrientation','landscape'); print(figure(4),['Best_T_in_P250.pdf'],'-dpdf','-bestfit')

figure(5)
plot([0:149],T_out, 'Linewidth', 2);

xlabel('Time[min]')
ylabel('Temperature [K]')

ylim([35, 55]+273)

fontsize(24,"points")
set(gcf,'PaperOrientation','landscape'); print(figure(5),['Best_T_out_P250.pdf'],'-dpdf','-bestfit')

close all

%}





















