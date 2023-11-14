startup;
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

%% Create the solver
Iteration_max               = 50;                                            % Maximum number of iterations for optimzer
Time_max                    = 20.0;                                           % Maximum time of optimization in [h]

nlp_opts                    = struct;
nlp_opts.ipopt.max_iter     = Iteration_max;
nlp_opts.ipopt.max_cpu_time = Time_max*3600;
%nlp_opts.ipopt.acceptable_tol  = 1e-4;
%nlp_opts.ipopt.acceptable_iter = 5;

OPT_solver                  = casadi.Opti();
ocp_opts                    = {'nlp_opts', nlp_opts};
OPT_solver.solver(             'ipopt'   , nlp_opts)

%% Load paramters
m_total                 = 80;

V_Flow                  = 0.4;                                              % Volumetric flow rate l/min

% Bed geometry
before                  = 0.1;                                              % Precentage of length before which is empty
bed                     = 0.165;                                            % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 60;
timeStep                = 1;                                               % Minutes
SamplingTime            = 10;                                              % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);
SAMPLE                  = SamplingTime:SamplingTime:ExtractionTime;
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

% Set operating conditions
T0homog                 = OPT_solver.variable(numel(N_Sample))';
OPT_solver.subject_to( 40+273 <= T0homog <= 50+273 );

feedPress               = 250;

Z                       = Compressibility( T0homog(1), feedPress,         Parameters );
rho                     = rhoPB_Comp(      T0homog(1), feedPress, Z,      Parameters );

enthalpy_rho            = rho.*SpecificEnthalpy(T0homog(1), feedPress, Z, rho, Parameters );

feedTemp                = repmat(T0homog,N_Time/numel(N_Sample),1);
feedTemp                = feedTemp(:)';

feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;  % Bars

feedFlow                = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/s

uu                      = [feedTemp', feedPress', feedFlow'];

%% Set inital state and inital conditions
msol_max                = m_total;                                          % g of product in solid and fluid phase
mSol_ratio              = 0.7;

mSOL_s                  = msol_max*mSol_ratio;                              % g of product in biomass
mSOL_f                  = msol_max*(1-mSol_ratio);                          % g of biomass in fluid

C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;                % Solid phase kg / m^3
Parameters{2}           = C0solid;

G                       =@(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;

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
choose_param            = [44:47];
% Store symbolic results of the simulation
Parameters_sym          = MX(cell2mat(Parameters));
Parameters_sym_t        = OPT_solver.parameter(numel(choose_param));
Parameters_sym(choose_param)   = Parameters_sym_t;

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
%FF                          = Function('FF', {[T0homog, Parameters_sym_t']}, {S});
%S                           = FF([cell2mat(Parameters(choose_param))']);

norm_factor             = cell2mat(Parameters(choose_param)) ./ repmat(data_obs,numel(choose_param),1);
S_norm                  = S .* abs(norm_factor)' ;

FI                      = 0;

for ii=1:numel(data_obs)
    FI                  = FI + (S_norm(ii,:) * S_norm(ii,:)');
end

D_opt = -myDet(FI);

OPT_solver.set_value(Parameters_sym_t, cell2mat(Parameters(choose_param)));

OPT_solver.minimize(D_opt);

OPT_solver.set_initial(T0homog, [linspace(50,40,numel(N_Sample))+273]);

try
    sol  = OPT_solver.solve();
    KOUT = full(sol.value(T0homog));
catch
    KOUT = OPT_solver.debug.value(T0homog);
end





























