clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.5');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

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


Tc      = 304.2;                                        % Critical temperature [K]
Pc      = 73.765;                                       % Critical pressure [bar]
omega   = 0.225;                                        % Acentric factor [-]
R       = 83.14;                                        % Universal gas constant, [cm3-bar/mol-K]

kappa   = 0.37464 + 1.5422 * omega - 0.2699 * omega^2;

MW      = 44.01;                                        % Molar weight of CO2 [Kg/mol]

[EA_Di, betah_Di] = Arrhenius_Di(Ti,Di,R);
[EA_km, betah_km] = Arrhenius_km(T_kp,dkp,R,rho_s);

CP_0 =  4.1868;
CP_A =  4.7280;        %4.5980;
CP_B =  0.0175;    %0.0125;
CP_C = -1.34E-05;    %2.86E-06;
CP_D =  4.10E-09;    %-2.70E-09;

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

%                 1        2    3    4    5  6     7   8   9   10  11  12   13   14   15      16       17    18    19    20    21    22      23        24      25        26       27      28       29        30        31       32       33       34      35       36        37      38       39       40      41        42       43
parameters = {nstages, C0solid, V, epsi, dp, L, rho_s, km, mi, Tc, Pc, R, kappa, MW, EA_Di, betah_Di, CP_0, CP_A, CP_B, CP_C, CP_D, EA_km, betah_km, cpSolid, a_axial, b_axial, c_axial, A1_cond, A2_cond, A3_cond, A4_cond, A5_cond, A6_cond, A7_cond, A1_visc, A2_visc, A3_visc, A4_visc, A5_visc, A6_visc, A7_visc, A8_visc, A9_visc};

%% load viscosityt data

rho_c = 467.6;

[dataset, title]  = xlsread("Viscosity.xlsx");
T_train           = dataset(:,1);
rho_train         = dataset(:,2);
Data              = dataset(:,3)*1e3;

%% Test RBF - optimze mu, weights and the location

RBF = 10;

rho_0   = rho_train(1)/rho_c + (rho_train(end)/rho_c-rho_train(1)/rho_c)*rand(RBF,1); 
T_0     = T_train(1)/Tc + (T_train(end)/Tc-T_train(1)/Tc)*rand(RBF,1);

rho_min = min(rho_train)*ones(RBF,1)/rho_c;
T_min   = min(T_train)*ones(RBF,1)/Tc;

rho_max = max(rho_train)*ones(RBF,1)/rho_c;
T_max   = max(T_train)*ones(RBF,1)/Tc;

mu      = MX.sym('mu'   , RBF );
rho_RBF = MX.sym('P_RBF', RBF );
T_RBF   = MX.sym('T_RBF', RBF );

w       = MX.sym('w', RBF );

APP     = MX(numel(T_train),1);

%[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [rho_RBF, T_RBF];

%% Evaluate the function at all the point

disp(['Create the matrix of data'])
for i = 1:numel(T_train)
    rho      = rho_train(i)/rho_c;
    T        = T_train(i)/Tc;
    distance = sqrt(( centers(:,1)-rho).^2 + (centers(:,2)-T).^2);
    APP(i)   = sum( w.*exp(- distance ./ mu) );
end

%% Set the optimzation problem to find mu
disp(['Set the structure of the NLP'])
nlp = struct;            % NLP declaration
nlp.x = [mu; w; rho_RBF; T_RBF];              % decision vars

D = (Data-APP).^2;
MSE = sum(D(:))/numel(Data);
%MSE = norm(Data-APP,'fro')^2/numel(Data);
nlp.f = MSE;               % objective - mean squared error

disp(['Create the solver'])
% Create solver instance
F = nlpsol('F','ipopt',nlp);


disp(['Solve NLP'])
% Solve the problem 

tic
res = F('x0',[ones(numel(mu),1); ones(numel(w),1); rho_0; T_0 ] ,'ubx',[inf*ones(numel(mu),1); inf*ones(numel(w),1); rho_max; T_max] ,'lbx',[ zeros(numel(mu),1); -inf*ones(numel(w),1); rho_min; T_min ] );
toc

%% Extract the solution
solution = full(res.x);

mu_opt        = solution(0*RBF+1:1*RBF);
Weights_opt   = solution(1*RBF+1:2*RBF);
rho_RBF_opt   = solution(2*RBF+1:3*RBF);
T_RBF_opt     = solution(3*RBF+1:4*RBF);

centers_opt = [rho_RBF_opt , T_RBF_opt ];

%%

APP_opt = nan(numel(T_train), 1);
disp(['Create the matrix of data'])
for i = 1:numel(T_train)
    rho = rho_train(i)/rho_c;
    T   = T_train(i)/Tc;
    distance_opt = sqrt(( centers_opt(:,1)-rho).^2 + (centers_opt(:,2)-T).^2);
    APP_opt(i) = sum(Weights_opt.*exp(- distance_opt ./ mu_opt));
end
%%
figure()
plot(APP_opt - Data)
ylabel('Difference between estiamted and the datapoint')
xlabel('Number of points')

%% 

T_check   = linspace(250,350,400);
rho_check = linspace(100,8000,500);

MU_test      = nan(numel(rho_check), numel(T_check));

for i = 1:numel(rho_check)

    if mod(i,round(numel(rho_check)/100)) == 0
        clc
        fprintf('%f %%\n',i/numel(rho_check)*100)
    end

    for j = 1:numel(T_check)

        rho = rho_check(i)/rho_c;
        T = T_check(j)/Tc;

        distance_opt = sqrt(( centers_opt(:,1)-rho/rho_c).^2 + (centers_opt(:,2)-T/Tc).^2);
        MU_test(i,j,1)   = sum(Weights_opt.*exp(- distance_opt ./ mu_opt));
    end
end

%%
figure()
hold on
contourf(T_check,rho_check,MU_test, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(T_RBF_opt*Tc, rho_RBF_opt*rho_c, 'k')
scatter(T_0*Tc, rho_0*rho_c, '*', 'k')
plot( [T_RBF_opt*Tc T_0*Tc ]' , [rho_RBF_opt*rho_c rho_0*rho_c]' ,'--k' )

hold off
ylabel('Density [kg/m3]')
xlabel('Temperature [K]')

%% Approximate the dataset

correlation_viscosity    = {'Amooey', 'Fenghour', 'Laesecke'};

T_check = linspace(250,350,400);
P_check = linspace(1,200,500);

Z       = nan(numel(P_check), numel(T_check), 3 );
RHO     = nan(numel(P_check), numel(T_check), 3 );
MU      = nan(numel(P_check), numel(T_check), 3, numel(correlation_viscosity) );
MU_RBF  = nan(numel(P_check), numel(T_check), 3 );

for i = 1:numel(P_check)

    if mod(i,round(numel(P_check)/100)) == 0
        clc
        fprintf('%f %%\n',i/numel(P_check)*100)
    end

    for j = 1:numel(T_check)

        T = T_check(j);
        P = P_check(i);

        Tr      = T ./ Tc;
        a       = 0.4572350 .* R.^2 .* Tc.^2 ./ Pc;
        b       = 0.0777961 .* R    .* Tc    ./ Pc;

        alpha = (1 + kappa .* (1 - sqrt(Tr))).^2;

        A = a .* alpha .* P ./ R.^2 ./ T.^2;
        B = b .* P ./ R ./ T;

        z = roots([1, - (1 - B), (A - 2.*B - 3.*B.^2), - ( A .* B - B.^2 - B.^3) ]);

        z = z(real(z)>0&imag(z)==0);

        if numel(z) == 1
            Z(i,j,1) = z;
        else
            Z(i,j,2) = max(z);
            Z(i,j,3) = min(z);
        end

        rho = full(rhoPB_Comp(T, P, z, parameters)) ;

        if numel(z) == 1
            RHO(i,j,1) = rho;
        else
            RHO(i,j,2) = max(rho);
            RHO(i,j,3) = min(rho);
        end

        for k=1:numel(correlation_viscosity)
            correlation_mass  = correlation_viscosity{k};
            if numel(z) == 1
                MU(i,j,1,k) = Viscosity(T,P,rho,parameters,correlation_mass);
            else
                MU(i,j,2,k) = Viscosity(T,P,max(rho),parameters,correlation_mass);
                MU(i,j,3,k) = Viscosity(T,P,min(rho),parameters,correlation_mass);
            end
        end

        if numel(z) == 1
            distance_opt = sqrt(( centers_opt(:,1)-rho/rho_c).^2 + (centers_opt(:,2)-T/Tc).^2);
            MU_RBF(i,j,1)   = sum(Weights_opt.*exp(- distance_opt ./ mu_opt));
        else
            distance_opt = sqrt(( centers_opt(:,1)-max(rho)/rho_c).^2 + (centers_opt(:,2)-T/Tc).^2);
            MU_RBF(i,j,2)   = sum(Weights_opt.*exp(- distance_opt ./ mu_opt));

            distance_opt = sqrt(( centers_opt(:,1)-min(rho)/rho_c).^2 + (centers_opt(:,2)-T/Tc).^2);
            MU_RBF(i,j,3)   = sum(Weights_opt.*exp(- distance_opt ./ mu_opt));
        end
    end
end


%%

root = 3;

clc
figure()
subplot(1,3,1)
hold on
contourf(T_check,P_check, MU(:,:,root,end), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
hold off
ylabel('Pressure [bar]')
xlabel('Temperature [K]')
title Correlation 

subplot(1,3,2)
contourf(T_check,P_check, MU_RBF(:,:,root), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
ylabel('Pressure [bar]')
xlabel('Temperature [K]')
title RBF

subplot(1,3,3)
contourf(T_check,P_check, MU_RBF(:,:,root)-MU(:,:,root,end), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
ylabel('Pressure [bar]')
xlabel('Temperature [K]')
title Difference

%% Evaluate opt mu and corresponding weights
%{
APP_opt = nan(numel(P_check), numel(T_check));

%mu = full(res.x);
%disp(['The optimal value of mu is ',num2str(mu)])

%% Approximate the dataset
for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;

        distance_opt = sqrt(( centers_opt(:,1)-P).^2 + (centers_opt(:,2)-T).^2);
        APP_opt(i,j) = sum(Weights_opt.*exp(- distance_opt ./ mu_opt));
    end
end

%% Plot the results with opt mu and compare to the original data

figure(1)
subplot(1,3,1);
contourf(T_check,P_check, Data, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
title('Original function')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2);
hold on
contourf(T_check,P_check,APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(T_RBF_opt*Tc, P_RBF_opt*Pc, 'k')
scatter(T_0*Tc, P_0*Pc, '*', 'k')
plot( [T_RBF_opt*Tc T_0*Tc ]' , [P_RBF_opt*Pc P_0*Pc]' ,'--k' )

%quiver( T_0*Tc, P_0*Pc, T_RBF_opt*Tc-T_0*Tc, P_RBF_opt*Pc-P_0*Pc, 0, 'k')

hold off
title('RBF')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3);
hold on
contourf(T_check,P_check,Data-APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(T_RBF_opt*Tc, P_RBF_opt*Pc, 'k')
scatter(T_0*Tc, P_0*Pc, '*', 'k')

hold off
title('Difference')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')
%}