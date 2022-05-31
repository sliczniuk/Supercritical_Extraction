clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

DATA = {'Data_40_200.mat', 'Data_50_200.mat', 'Data_50_300.mat', 'Data_40_300.mat'};

%% Parameters
nstages = 5;                                                   %

C0fluid = 0;                                            % Extractor initial concentration of extract
% Fluid phase kg / m^3
C0solid = 32;                                           % Solid phase kg / m^3
V       = 0.01/4;                                       % Volume of the extractor  [m3] - Vargas
epsi    = 0.70;                                         % Porosity [-] - Vargas
dp      = 0.001;                                        % Diameter of the particle [m] - Vargas
L       = 0.70;                                         % Length of the extractor [m] - Vargas
rho_s   = 1087.2;                                       % Densisty of the solid phase [kg / m^3] - Vargas
km      = 0.29;                                         % Partition coefficient (?)
mi      = 1/5;                                          % Geometric shape coefficient (?) - Vargas
Ti      = [313.15, 323.15, 333.15, 343.15];
Di      = [ 3.03E-11, 3.66E-11, 1.97E-10, 2.73E-10];
T_kp    = [313.15, 323.15, 333.15, 343.15];
dkp     = [8.13E-01 6.94E-01 1.04E-01 6.67E-02];


Tc      = 304.1;                                        % Critical temperature [K]
Pc      = 73.8;                                         % Critical pressure [bar]
omega   = 0.228;                                        % Acentric factor [-]
R       = 8.314e-5;                                     % Universal gas constant, [m3-bar/K-mol]

kappa   = 0.37464 + 1.54226 * omega - 0.26992 * omega^2;

MW      = 44e-3;                                        % Molar weight of CO2 [Kg/mol]

[EA_Di, betah_Di] = Arrhenius_Di(Ti,Di,R);
[EA_km, betah_km] = Arrhenius_km(T_kp,dkp,R,rho_s);

CP_0 =   4.18686;
CP_A =  19.80;        %4.5980;
CP_B =   7.344E-2;    %0.0125;
CP_C = - 5.602E-5;    %2.86E-06;
CP_D =   1.715E-8;    %-2.70E-09;

cpSolid = 1.5E3;      % J / K / Kg

% axial diffusion correlation
a_axial = 0.1;
b_axial = 0.011;
c_axial = 0.48;

% Heat Conductivity
A1_cond = -105.161;
A2_cond =  0.9007;
A3_cond =  0.0007;
A4_cond =  3.50E-15;
A5_cond =  3.76E-10;
A6_cond =  0.7500;
A7_cond =  0.0017;

% Viscosity
A1_visc = -1.146067E-1;
A2_visc =  6.978380E-7;
A3_visc =  3.976765E-10;
A4_visc =  6.336120E-2;
A5_visc = -1.166119E-2;
A6_visc =  7.142596E-4;
A7_visc =  6.519333E-6;
A8_visc = -3.567559E-1;
A9_visc =  3.180473E-2;

sigma   = 1;
sigma_km = 2;
mu_km    = 1;
Di = 7e-12;
Dx = 0.0001;

%                 1        2    3    4    5  6     7   8   9   10  11  12   13   14   15      16       17    18    19    20    21    22      23        24      25        26       27      28       29        30        31       32       33       34      35       36        37      38       39       40      41        42       43      44   45   46
parameters = {nstages, C0solid, V, epsi, dp, L, rho_s, km, mi, Tc, Pc, R, kappa, MW, EA_Di, betah_Di, CP_0, CP_A, CP_B, CP_C, CP_D, EA_km, betah_km, cpSolid, a_axial, b_axial, c_axial, A1_cond, A2_cond, A3_cond, A4_cond, A5_cond, A6_cond, A7_cond, A1_visc, A2_visc, A3_visc, A4_visc, A5_visc, A6_visc, A7_visc, A8_visc, A9_visc, Di, Dx ,sigma };
which_parameter= {15, 16, 22, 23, 25, 26, 27, 46};
theta = parameters;

%% Time
simulationTime       = 150;                                                 % Minutes
delayTime            = 2;                                                   % Minutes
timeStep             = 1;                                                 % Minutes
timeStep_in_sec      = timeStep *60;                                        % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;          % Seconds
Time                 = [0 Time_in_sec/60];                                   % Minutes
SamplingTime         = 5;                                                  % Minutes
N_Sample             = find(abs(Time-SamplingTime) < timeStep/10)-1;
N_Delay              = delayTime / timeStep;

%% Casadi variables and Model

% Create symbolic variables
u  = MX.sym('u', 3);
x  = MX.sym('x', 3*nstages+1);
k  = MX.sym('k', length(which_parameter));

%Variables
Nx = 3*nstages+1;
Nu = 3;
Nk = length(which_parameter);
Ny = 1;

%% Model
% Model Equations
f = @(x, u, k) modelSFE_SS(x, u, k, parameters, which_parameter);
%g = @(x) modelSFE_out(x, parameters);
g = @(x, u, y_old) modelSFE_out2(x, u, y_old, parameters, timeStep_in_sec);

% Integrator
F = buildIntegrator_ParameterEstimation(f, [Nx,Nu,Nk] , timeStep_in_sec);

%% dummy parameters
RHO = [];
k0 = [1;50;1;1;1;1;1;1];

%%
OCP_solver = casadi.Opti();

% http://www.diva-portal.se/smash/get/diva2:956377/FULLTEXT01.pdf
nlp_opts = struct;
nlp_opts.ipopt.max_iter = 100;
%nlp_opts.ipopt.acceptable_iter = 50;
%nlp_opts.ipopt.acceptable_tol = 1e-6;
%nlp_opts.ipopt.tol = 1e-7;

ocp_opts = {'nlp_opts', nlp_opts};
OCP_solver.solver('ipopt',nlp_opts)

K = OCP_solver.variable(numel(which_parameter));

Cost = 0;

%%   
for l = 1:1%numel(DATA)
    l
    load(DATA{l});

    feedTemp  = T0homog      * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress(1) * ones(1,length(Time_in_sec));  % Bars

    Z = Compressibility(T0homog, feedPress(1), parameters);
    rho = full(rhoPB_Comp(T0homog, feedPress(1), Z, parameters));
    RHO = [RHO, rho];

    feedFlow  = 0.4 * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec

    %% Inital values
    x0 = [C0fluid*ones(nstages,1);
        C0solid*ones(nstages,1);
        T0homog*ones(nstages,1);
        0];

    uu = [feedTemp', feedPress', feedFlow'];

    data2 = diff(data);
    

    %%
    OCP = struct('Nx', Nx, 'Nu', Nu, 'Nk', Nk, 'Ny', Ny, 'N', length(Time_in_sec), ...
        'x_lu', [],  ...                % 0*ones(Nx,1) inf*ones(Nx,1)
        'k_lu', [],  ...                % 0*ones(Nk,1) [1.5;inf;inf]; 0*ones(Nk,1) inf*ones(Nk,1)
        'x_eq', [], ...
        'k_eq', [], ...
        'N_Sample', N_Sample,...
        'N_Delay', N_Delay,...
        'F'   , F, ...
        'L'   , [], ...       %
        'LS'  , @(x) sum( (x-data2) .^2 ),...
        'MSE', @(x,sigma_MSE) -(-length(data2)/2*log(2*pi*sigma_MSE^2) - 1/(2*sigma_MSE^2) * sum( (x-data2) .^2 )), ...                 %'Lf'  , @(x,yd) sum( (g(x)-yd).^2 )   );   % 0.5*( sum(g(x)-yd) )'*( sum(g(x)-yd) ) %'Lf'  , @(x,yd) sum( (g(x)-yd).^2 )   );   % 0.5*( sum(g(x)-yd) )'*( sum(g(x)-yd) )
        'MAP', @(x,km,sigma_MSE) -( -length(data2)/2*log(2*pi*sigma_MSE^2) - 1/(2*sigma_MSE^2) * sum( (x-data2) .^2 ) - 1/2*log(2*pi*sigma_km^2) - 1/(2*sigma_km^2) * sum( (km-mu_km) .^2 )  ) );


    %% Solve LSE
    [J] = singleShooting_ParameterEstimation(OCP, x0, uu',  K,  parameters);
    
    Cost = Cost + J;

end

% state constraints
if ~isempty(OCP.x_lu)
    for nx=1:OCP.Nx
        OCP_solver.subject_to(OCP.x_lu(nx,1)<=X(nx,:)<= OCP.x_lu(nx,2));
    end
end

% parameter constraints
if ~isempty(OCP.k_lu)
    for nk=1:OCP.Nk
        OCP_solver.subject_to(OCP.k_lu(nk,1)<=K(nk,:)<= OCP.k_lu(nk,2));
    end
end

OCP_solver.minimize(Cost);
OCP_solver.set_initial(K, k0);

try
    sol = OCP_solver.solve();
    kout = sol.value(K);
catch
    kout = OCP_solver.debug.value(K);
end

%% Substitute the initial parameters with their estimates
    if ~isempty(which_parameter)
        for i=1:length(which_parameter)
            theta{which_parameter{i}} = kout(i);
        end
    end

    %
%     fprintf('\n-------------------------------------------------\n')
%     fprintf('           | EA_Di | beath_Di |  Di   |  km  \n')
%     fprintf('-------------------------------------------------\n')
%     fprintf(' Vargas    | %-2.0d |  %-2.0d   | %-2.0d | %-2.0d\n', EA_Di, betah_Di,  Di_of_T(T0homog, parameters), km )
%     fprintf(' Estiamted | %-2.0d |  %-2.0d   | %-2.0d | %-2.0d\n', kout(2), kout(3) , Di_of_T(T0homog, theta), kout(1) )
%     fprintf('-------------------------------------------------\n')

%%

for l = 1:1%numel(DATA)
    %% Load data
    load(DATA{l});

    feedTemp  = T0homog      * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress(1) * ones(1,length(Time_in_sec));  % Bars

    Z = Compressibility(T0homog, feedPress(1), parameters);
    rho = full(rhoPB_Comp(T0homog, feedPress(1), Z, parameters));
    RHO = [RHO, rho];

    feedFlow  = 0.4 * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec

    %% Inital values
    x0 = [C0fluid*ones(nstages,1);
        C0solid*ones(nstages,1);
        T0homog*ones(nstages,1);
        0];

    uu = [feedTemp', feedPress', feedFlow'];

    data2 = diff(data);

    %% Simulate system with estiamted parameter
    [yy_out, tt_out, xx_out] = simulateSystem_ParameterEstimation(F, g, x0, uu, kout);

    %% Confidence interval
    % Confidence 80   | 85    | 90    | 95    | 99    | 99.5  | 99.9
    % Interval  1,282 | 1,440 | 1,645 | 1,960 | 2,576 | 2,807 | 3,291

    %CI = [1.282, 1.440, 1.645, 1.960, 2.576, 2.807, 3.291];
    %figure()

    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)], ...
        'renderer','OpenGL');

    yy_out_D = [ zeros(1,N_Delay) yy_out(1:end-N_Delay)];

    CI = 0;

    hold on;
    while CI <= 1.96
        e = CI * kout(end) / sqrt(1);
        lo = yy_out_D - e;
        hi = yy_out_D + e;

        hp = patch([Time, Time(end:-1:1), Time(1)], [lo, hi(end:-1:1), lo(1)], 'b','FaceAlpha',0.01);
        set(hp, 'facecolor', [0, 0, 1], 'edgecolor', 'none');
        CI = CI + 1.96/50;
    end
    plot(Time, yy_out_D, 'color', 'k', 'LineWidth',1);
    plot(Time(1:N_Sample:end), data, 'o','Color','k','LineWidth',1)
    hold off
    xlim([0,Time(end)])
    ylim([-10,110])
    %set(hl, 'color', 'r', 'marker', 'x');
    %set(gcf,'renderer','OpenGL');

    %caption1 = sprintf('T = %g [K], P = %g [bar]\n EA_{Di} = %g, D^0_{Di} = %g, EA_{km} = %g, D^0_{km} = %g\na = %g, b = %g, c = %g, \\sigma = %g', T0homog, feedPress(1), kout');
    %caption1 = sprintf('T = %g [K], P = %g [bar]\n km = %g, EA_{Di} = %g, D^0_{Di} = %g\na = %g, b = %g, c = %g, \\sigma = %g', T0homog, feedPress(1), kout');
    caption1 = sprintf('T = %g [K], P = %g [bar]\n km = %g, EA_{Di} = %g, D^0_{Di} = %g, \\sigma = %g', T0homog, feedPress(1), kout');
    title(caption1, 'Interpreter','tex');

    xlabel('Time [min]')
    ylabel('Yield [%]')

    Name = strcat('MegaFitting_T_',string(T0homog-273.15),'_P_',string(feedPress(1)),'.pdf');

    %print(gcf, Name,'-dpdf','-r700');

end
%%
%}