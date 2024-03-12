startup;
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

excel_file = 'Chamomile_Di_Gamma_2.xls';
rng(69)

%%
Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

LabResults              = xlsread('dataset_2.xlsx');

N_trial                 = 1;
Iteration_max           = 15;                                              % Maximum number of iterations for optimzer
Time_max                = 12;                                               % Maximum time of optimization in [h]

nlp_opts                    = struct;
nlp_opts.ipopt.max_iter     = Iteration_max;
nlp_opts.ipopt.max_cpu_time = Time_max*3600;
nlp_opts.ipopt.acceptable_tol     = 1e-4;
nlp_opts.ipopt.tol      = 1e-6;
%nlp_opts.ipopt.hessian_approximation ='limited-memory';

which_k                 = [      44,   46];              % Select which parameters are used for fitting
k_lu                    = [ [     0;    0], ...
                            [   inf;  inf] ];
Nk                      = numel(which_k)+0;                                   % Parameters within the model + sigma

DATA_K_OUT              = nan(Nk,N_trial);                          % Store Parameters obatined from all fits (par num x num exper)
OBJ_OUT                 = nan(1,N_trial);

AA = xlsread('Regression.xlsx');
DI = [0.701472275, 1.331443055, 2.239307889, 2.711813187, 1.32629228, 1.485504345, 1.73827467, 2.59502961, 0.48656241, 1.363499511, 0.72227, 0.756214019];
GG = [4.274825704, 2.189390368, 2.552240039, 1.365163176, 2.830760407, 2.573487374, 1.642279591, 1.906200052, 4.287215235, 2.723682117, 3.82240, 3.35589348];

%% Load paramters
m_total                 = 3.0;

% Bed geometry
before                  = 0.04;                                             % Precentage of length before which is empty
bed                     = 0.92;                                              % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 600;
timeStep                = 10;                                                % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);

SAMPLE                  = LabResults(6:19,1);

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

for jj = 1
    which_dataset           = jj;

    data_org                = LabResults(6:19,which_dataset+1)';
    data_diff               = diff(data_org);

    %noise                  = -0.1 + (0.1 + 0.1) .* rand(1,numel(data_diff));
    %data_diff              = data_diff + noise;

    %% Set operating conditions
    T0homog                 = LabResults(2,which_dataset+1);                    % K
    feedPress               = LabResults(3,which_dataset+1) * 10;               % MPa -> bar
    Flow                    = LabResults(4,which_dataset+1) * 1e-5 ;            % kg/s

    Z                       = Compressibility( T0homog, feedPress,         Parameters );
    rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );

    enthalpy_rho            = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );

    feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;     % Kelvin

    feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;     % Bars

    feedFlow                = Flow * ones(1,length(Time_in_sec));               % kg/s

    uu                      = [feedTemp', feedPress', feedFlow'];

    % Initial conditions
    x0                      = [ C0fluid'                         ;
                                C0solid         * bed_mask       ;
                                enthalpy_rho    * ones(nstages,1);
                                feedPress(1)                     ;
                                0                                ;
                                ];

    %Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];
    %[xx_0]                 = simulateSystem(F, [], x0, Parameters_init_time );
    %{\
    %%
    %for ii=1:N_trial
        %k0                      = 0.01 + (100-0.01) .* rand(Nk,1);

        %k0                      = ones(Nk,1);
        k0                      = [0.72, 2.6];

        OPT_solver              = casadi.Opti();
        ocp_opts                = {'nlp_opts', nlp_opts};
        OPT_solver.solver(         'ipopt'   , nlp_opts)

        % Descision variables
        k                       = OPT_solver.variable(Nk);

        % Store symbolic results of the simulation
        Parameters_sym          = MX(cell2mat(Parameters));
        Parameters_sym(which_k) = k(1:numel(which_k));

        X                       = MX(Nx,N_Time+1);
        X(:,1)                  = x0;

        % Symbolic integration
        for j=1:N_Time
            X(:,j+1)=F(X(:,j), [uu(j,:)'; Parameters_sym] );
        end

        %% Find the measurment from the simulation
        Yield_estimate         = X(Nx,N_Sample);
        Yield_estimate_diff    = diff(Yield_estimate);

        %% load dataset
        Yield_estimate_diff_norm = Yield_estimate_diff ./ max(data_diff) ;
        data_diff_norm           = data_diff           ./ max(data_diff) ;
        residuals                = (data_diff_norm - Yield_estimate_diff_norm );

        %% Create the cost function
        %sigma                  = k0(end);
        J                      = residuals * diag(1) * residuals';
        sigma                  = J/(numel(residuals) - numel(which_k));

        J_L                    = - (numel(residuals)./2 .* log(2.*pi)) - (numel(residuals)./2 .* log(sigma)) - J./(2.*sigma);  % normal
        J_L                    = - J_L;

        %% Constraints
        for nk=1:Nk
            OPT_solver.subject_to( k_lu(nk,1) <= k(nk,:) <= k_lu(nk,2) );
        end

        OPT_solver.minimize(J_L);

        OPT_solver.set_initial(k, k0)

        try
            sol = OPT_solver.solve();
            KOUT = full(sol.value(k));
        catch
            KOUT = OPT_solver.debug.value(k)
        end

        OBJ = OPT_solver.stats.iterations.obj;

        %DATA_K_OUT(:,ii) = KOUT;
        %OBJ_OUT(ii)      = OBJ(end);

%    end

    %% Save data
    %writematrix([OBJ(end); KOUT],excel_file,'sheet',which_dataset,'Range','A1:Z200');
    %}
%end
%% Plots
%{
for ii=1:N_trial
    %KOUT = DATA_K_OUT(:,ii);
    KOUT = KOUT;
    %KOUT = [ 0.72, 2.5 ];
    %KOUT = [ DI(jj) , GG(jj) ];
    Parameters_opt = Parameters;
    for i=1:numel(which_k)
        Parameters_opt{which_k(i)}  = KOUT(i);
    end

    Parameters_init_time   = [uu repmat(cell2mat(Parameters_opt),1,N_Time)'];
    [xx_0]                 = simulateSystem(F, [], x0, Parameters_init_time );

    data_org = cumsum([0 data_diff]);
    figure(jj)
    subplot(2,1,1)
    hold on; plot(Time,xx_0(end,:), 'LineWidth',2, 'DisplayName', [num2str(AA(2,jj)),'[K],',num2str(AA(3,jj)),'[MPa]']); plot(SAMPLE, data_org,'ko', 'LineWidth',2, 'HandleVisibility','off' ); hold off
    subplot(2,1,2)
    hold on; plot(SAMPLE(2:end),diff(xx_0(end,N_Sample))); plot(SAMPLE(2:end), data_diff,'o'); hold off
    xlabel('t [min]')
    ylabel('$\frac{dy}{dt}~\left[ \frac{g}{min} \right]$')
    %set(gcf,'PaperOrientation','landscape'); print(figure(1),['Fit_Di_Gamma_dataset_',num2str(jj),'_org.pdf'],'-dpdf','-bestfit')
    %close all
    
end
end
%xlabel('t [min]')
%ylabel('y [g]')

%legend boxoff 
%lgd = legend('Location','northoutside', 'Orientation','horizontal');
%lgd.FontSize = 10;
%set(gca,'FontSize',12)
%}

%% plot the parameter space
%{\
import casadi.*
FF = Function('FF',{k},{J_L});

XX = linspace(0.1,2.5,50);
YY = linspace(0.1,15,60);
[X,Y] = meshgrid(XX,YY);

Z = nan(numel(XX),numel(YY)); 

for ii=1:numel(XX)
    ii
    parfor jj=1:numel(YY)
        Z(ii,jj) = full(FF([XX(ii),YY(jj)]));
    end
end

colormap turbo
[M,c] = contourf(X,Y,Z',100);
c.LineColor = 'none';

xlabel('$D_i^R \times 10^{-13}$')
ylabel('$\Upsilon$')

hold on
plot(KOUT(1),KOUT(2),'ro')
hold off

hbar = colorbar;
ylabel(hbar,'$-Ln(L)$', 'interpreter', 'latex');
set(gca,'FontSize',12)

%yscale log

exportgraphics(figure(1), ['Parameter_space_Di_Gamma_dataset_',num2str(jj),'_org.png'], "Resolution",300);
%exportgraphics(figure(1), 'output.pdf', 'ContentType', 'vector');
%set(gcf,'PaperOrientation','landscape'); print(figure(1),['Parameter_space_Linear_dataset_',num2str(jj),'_org.pdf'],'-dpdf','-bestfit')
%close all

end
%}














