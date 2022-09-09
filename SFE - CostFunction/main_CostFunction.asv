clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*
%load cost.mat
DATA = {'LUKE_T40_P200.xlsx', 'LUKE_T50_P200.xlsx', 'LUKE_T40_P300.xlsx', 'LUKE_T50_P300.xlsx'};

%% Parameters
nstages = 100;                                           %

m_ref   = 80;                                           % g of product obtained from a kg of biomass

C0fluid = 0;                                            % Extractor initial concentration of extract
% Fluid phase kg / m^3

V       = 0.00165;                                       % Volume of the extractor  [m3] - Vargas
L       = 0.095;                                         % Length of the extractor [m] - Vargas
epsi    = 0.75;                                          % Porosity [-] - Vargas

C0solid = m_ref *1e-3 / (V*(1-epsi));                   % Solid phase kg / m^3

dp      = 0.0010;                                        % Diameter of the particle [m] - Vargas
rho_s   = 1300.0;                                       % Densisty of the solid phase [kg / m^3] - FC
km      = 0.29;                                         % Partition coefficient (?)

mi      = 1/2;                                          % Geometric shape coefficient (?) - Vargas

Ti      = [313.15, 323.15, 333.15, 343.15];
Di      = [3.03E-11, 3.66E-11, 1.97E-10, 2.73E-10];
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

% Axials diffusions correlations
a_axial = 0.1;
b_axial = 0.011;
c_axial = 0.48;

% Heats Conductivitys
A1_cond = -105.161;
A2_cond =  0.9007;
A3_cond =  0.0007;
A4_cond =  3.50E-15;
A5_cond =  3.76E-10;
A6_cond =  0.7500;
A7_cond =  0.0017;

% Viscositys
A1_visc = -1.146067E-1;
A2_visc =  6.978380E-7;
A3_visc =  3.976765E-10;
A4_visc =  6.336120E-2;
A5_visc = -1.166119E-2;
A6_visc =  7.142596E-4;
A7_visc =  6.519333E-6;
A8_visc = -3.567559E-1;
A9_visc =  3.180473E-2;

sigma    = 1;
sigma_km = 0.1;
mu_km    = 0.1;
Di_slack = 1;

%                 1        2    3    4    5  6     7   8   9   10  11  12   13   14   15      16       17    18    19    20    21    22      23        24      25        26       27      28       29        30        31       32       33       34      35       36        37      38       39       40      41        42       43      44        45
parameters = {nstages, C0solid, V, epsi, dp, L, rho_s, km, mi, Tc, Pc, R, kappa, MW, EA_Di, betah_Di, CP_0, CP_A, CP_B, CP_C, CP_D, EA_km, betah_km, cpSolid, a_axial, b_axial, c_axial, A1_cond, A2_cond, A3_cond, A4_cond, A5_cond, A6_cond, A7_cond, A1_visc, A2_visc, A3_visc, A4_visc, A5_visc, A6_visc, A7_visc, A8_visc, A9_visc, Di_slack, sigma };


%% Times
simulationTime       = 150;                                                 % Minutes
delayTime            = 0;                                                   % Minutes
timeStep             = 1;                                                   % Minutes
timeStep_in_sec      = timeStep * 60;                                       % Seconds
Time_in_sec          = (timeStep:timeStep:simulationTime)*60;               % Seconds
Time                 = [0 Time_in_sec/60];                                  % Minutes
SamplingTime         = 5;                                                   % Minutes
N_Sample             = find(abs(Time-SamplingTime) < timeStep/10)-1;
N_Delay              = delayTime / timeStep;

%% Casadis variabless ands Models

% Creates symbolics variabless
u  = MX.sym('u', 3);
x  = MX.sym('x', 3*nstages+1);
k  = MX.sym('k', 3);

%Variabless
Nx = numel(x);
Nu = numel(u);
Nk = numel(k);
Ny = 1;

%% Models
f_r = @(x, u, k) modelSFE_Regression(x, u, k, parameters);
g   = @(x, u, y_old) modelSFE_out2(x, u, y_old, parameters, timeStep_in_sec);
% Integrator
F_r = buildIntegrator_ParameterEstimation(f_r, [Nx,Nu,Nk] , timeStep_in_sec);

%%

km_check =  logspace(-02, +00, 3);
Di_check =  logspace(-13, -10, 4);
Dx_check =  logspace(-05, -03, 3);


%%

LS = zeros(numel(km_check), numel(Di_check),numel(Dx_check),4);

RHO = [];

V_Flow = 0.30;

id = [1 2 3 4];

for i = id% 1:numel(DATA)
    %load(DATA{i});

    LabResults = xlsread(DATA{i});

    T0homog   = LabResults(1,1)+273.15;
    feedPress = LabResults(1,2);
    data      = LabResults(:,4)';
 
    feedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars

    Z         = Compressibility(T0homog(1), feedPress(1), parameters);
    rho       = full(rhoPB_Comp(T0homog(1), feedPress(1), Z, parameters));

    RHO       = [RHO; rho];

    feedFlow  = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec


    %% Inital values
    x0 = [C0fluid*ones(nstages,1);
        C0solid*ones(nstages,1);
        T0homog*ones(nstages,1);
        0];

    uu = [feedTemp', feedPress', feedFlow'];

    for j=1:numel(km_check)
        [i, j]
        tic
        for k=1:numel(Di_check)

            parfor l=1:numel(Dx_check)
            %[i,j,k]
            
                [yy_out, tt_out, xx_out] = simulateSystem_ParameterEstimation(F_r, g, x0, uu, [km_check(j); Di_check(k); Dx_check(l)] );
                
    
                Yield_estimate = xx_out(3*nstages+1,:);
                %Yield_estim\te = yy_out;
                Yield_estimate = [zeros(1,N_Delay) Yield_estimate(1:end-N_Delay)];
                Yield_estimate = Yield_estimate(1:N_Sample:end);
    
                LS(j,k,l,i) = (data-Yield_estimate) * diag(1:1:1)  * (data-Yield_estimate)';

            end

        end
        toc
    end
    clc
end
%%
%save cost4.mat

%% Plot LS
%{
figure(1)

for i=id%1:numel(DATA)
        
        LabResults = xlsread(DATA{i});

        T0homog   = LabResults(1,1)+273.15;
        feedPress = LabResults(1,2);
       
       
        subplot(2,2,i)
        contourf(log10( Di_check ), log10( km_check) ,log10( squeeze( LS(:,:,i) ) ),50, 'Edgecolor','none'); colormap cool;
        [r c] = find(squeeze( LS(:,:,i) == min(min( squeeze( LS(:,:,i) ) ) )) );
        hold on
        s = scatter( log10( Di_check(c) ), log10( km_check(r) ),log10( squeeze( LS(r,c,i) )  ) ,'ko');
        s.SizeData = 100;
        %yline(log10( km_check(r) ),'k');
        %xline(log10( Di_check(c) ),'k');
        %datatip(s,Di_check(c),km_check(r));
        hold off
        xlabel('Di');
        ylabel('km');
        caption = sprintf('T = %g, P = %g, \\rho = %g', T0homog, feedPress, RHO(i) );
        title(caption)
       
end


%% Yield curves with the minimum values of LS
clc
%YIELD = nan(4,numel(Time),2);
figure(2);
%}

for i = id%1:1%numel(DATA)
    LabResults = xlsread(DATA{i});

    T0homog   = LabResults(1,1)+273.15;
    feedPress = LabResults(1,2);
    data      = LabResults(:,4);

    eedTemp  = T0homog   * ones(1,length(Time_in_sec));  % Kelvin
    feedPress = feedPress * ones(1,length(Time_in_sec));  % Bars

    feedFlow  = V_Flow * RHO(i) * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/min -> Kg / sec
    
    [m inj] = min(LS(:,:,:,i),[],"all");
    [a1 a2 a3] = ind2sub(size(squeeze(LS(:,:,:,i))),inj);

    
    %% Inital values
    x0 = [C0fluid*ones(nstages,1);
        C0solid*ones(nstages,1);
        T0homog*ones(nstages,1);
        0];

    uu = [feedTemp', feedPress', feedFlow'];

    subplot(2,2,i)
    hold on;

    [yy_out, tt_out, xx_out] = simulateSystem_ParameterEstimation(F_r, g, x0, uu, [km_check(a1), Di_check(a2), Dx_check(a3)] );

    yy_out_D = [zeros(1,N_Delay) yy_out(1:end-N_Delay)];

    plot(Time,xx_out(end,:)); %plot(Time,yy_out); plot(Time,yy_out_D);
    %plot(Time,yy_out_D);

    
    plot(Time(1:N_Sample:end),data,'o'); hold off
    ylabel('Mass of product [g]')
    xlabel('Time [min]')
    caption = sprintf('T = %g, P = %g, \\rho = %g', T0homog, feedPress(1), RHO(i) );
    title(caption);
    

end
%figure()
%hold on; plot(Time,xx_out(end,:)); plot(Time,yy_out); plot(Time,yy_out_D); plot(Time(1:50:end),data,'o'); hold off
