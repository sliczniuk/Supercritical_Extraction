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


%% dummy parameters
clc

correlation_viscosity    = {'Amooey', 'Fenghour', 'Laesecke'};
correlation_conductivity = {'Amooey', 'Bahadori', 'Jarrahian', 'Rostami', 'Rostamian', 'Huber'};

T_check = linspace(Tc,350,50);
P_check = linspace(Pc,100,50);
%T_check = Tc;
%P_check = Pc;
%T_check = [0.85*Tc, 0.95*Tc, 0.99*Tc,     Tc,  1.01*Tc, 1.05*Tc, ];
%P_check = [0.01*Pc, 0.10*Pc, 0.25*Pc, 0.50*Pc, 0.75*Pc,      Pc, 1.10*Pc, 1.5*Pc];

N_polynomial = linspace(0,1,100);
pol = nan(numel(P_check), numel(T_check), numel(N_polynomial) );
Z   = nan(numel(P_check), numel(T_check), 3 );
RHO = nan(numel(P_check), numel(T_check), 3 );
CP  = nan(numel(P_check), numel(T_check), 3 );
MU  = nan(numel(P_check), numel(T_check), 3, numel(correlation_viscosity) );
KT  = nan(numel(P_check), numel(T_check), 3, numel(correlation_conductivity) );

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

        f = @(z) z.^3 - (1 - B)*z.^2 + (A - 2.*B - 3.*B.^2).*z + - ( A .* B - B.^2 - B.^3);
        pol(i,j,:) = f(N_polynomial);

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

        if numel(z) == 1
            CP(i,j,1) = SpecificHeatComp(T, P, z, rho, parameters);
        else
            CP(i,j,2) = SpecificHeatComp(T, P, max(z), max(rho), parameters);
            CP(i,j,3) = SpecificHeatComp(T, P, min(z), min(rho), parameters);
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

        for k=1:numel(correlation_conductivity)
            correlation_heat = correlation_conductivity{k};
            if numel(z) == 1
                KT(i,j,1, k) = HeatConductivity_Comp(T, P, z, rho, MU(i,j,1,end), parameters, correlation_heat);
            else
                KT(i,j,2, k) = HeatConductivity_Comp(T, P, max(z), max(rho), MU(i,j,2,end), parameters, correlation_heat);
                KT(i,j,3, k) = HeatConductivity_Comp(T, P, min(z), min(rho), MU(i,j,3,end), parameters, correlation_heat);
            end
        end

    end
end


%% Plotting
clc
%{
rho_check = 1:1200;
T_check = 303:0.01:306;

KT  = nan(numel(rho_check), numel(T_check),4);

for i = 1:numel(rho_check)
    if mod(i,round(numel(rho_check)/100)) == 0
        clc
        fprintf('%f %%\n',i/numel(rho_check)*100)
    end
    for j = 1:numel(T_check)

        T = T_check(j);
        rho = rho_check(i);

        KT(i,j,:) = HeatConductivity_Comp(T, Pc, 0.53, rho, [], parameters, correlation_conductivity{6});
    end
end
%plot(rho_check, KT)

%%
subplot(4,1,1)
plot(rho_check,squeeze(KT(:,:,1)))
subplot(4,1,2)
plot(rho_check,squeeze(KT(:,:,2)))
subplot(4,1,3)
plot(rho_check,squeeze(KT(:,:,3)))
subplot(4,1,4)
plot(rho_check,squeeze(KT(:,:,4)))
%}

%%

%{

% Compressibility

figure(1);
set(gcf,'PaperOrientation','landscape', 'visible','off')
subplot(1,3,1); contourf(T_check,P_check,squeeze(Z(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'Z');
title('Compressibility 1-phase region',' ')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2); contourf(T_check,P_check,squeeze(Z(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'Z');
title('Compressibility 2-phase region:','the biggest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3); contourf(T_check,P_check,squeeze(Z(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'Z');
title('Compressibility 2-phase region:','the smallest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'Compressibility.pdf'); close all

% Rho

figure(2)
set(gcf,'PaperOrientation','landscape', 'visible','off')
subplot(1,3,1); contourf(T_check,P_check,squeeze(RHO(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'\rho_{CO_2} [kg/m^3]');
title('Density 1-phase region',' ')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2); contourf(T_check,P_check,squeeze(RHO(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'\rho_{CO_2} [kg/m^3]');
title('Density 2-phase region:','the biggest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3); contourf(T_check,P_check,squeeze(RHO(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'\rho_{CO_2} [kg/m^3]');
title('Density 2-phase region:','the smallest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'RHO.pdf'); close all

% CP

figure(3)
mCP = log10(CP);
set(gcf,'PaperOrientation','landscape', 'visible','off')
subplot(1,3,1); contourf(T_check,P_check,squeeze(mCP(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat 1-phase region:',' ')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2); contourf(T_check,P_check,squeeze(mCP(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat 2-phase region:','the biggest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3); contourf(T_check,P_check,squeeze(mCP(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat 2-phase region:','the smallest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'CP.pdf'); close all

% mu

figure(4)
set(gcf, 'visible','off')

for i=0:numel(correlation_viscosity)-1
    subplot(numel(correlation_viscosity),3,3*i+1); contourf(T_check,P_check,squeeze(MU(:,:,1,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'\mu[Pa \times s]');
    title('Viscosity 1-phase region:',' ')
    
    ylabel({correlation_viscosity{i+1};'Pressure [bar]'})
    xlabel('Temperature [K]')
    
    subplot(numel(correlation_viscosity),3,3*i+2); contourf(T_check,P_check,squeeze(MU(:,:,2,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'\mu[Pa \times s]');
    title('Viscosity 2-phase region:','the biggest root')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
    
    subplot(numel(correlation_viscosity),3,3*i+3); contourf(T_check,P_check,squeeze(MU(:,:,3,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'\mu[Pa \times s]');
    title('Viscosity 2-phase region:','the smallest root')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
end
print(gcf, '-dpdf', '-fillpage', 'MU.pdf'); close all

%}

%{

% kt

figure(5)
set(gcf, 'visible','on')

for i=0:numel(correlation_conductivity)-1
    subplot(numel(correlation_conductivity),3,3*i+1); contourf(T_check,P_check,squeeze(KT(:,:,1,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'kt');
    title('Conductivity 1-phase region:',' ')
    ylabel({correlation_conductivity{i+1};'Pressure [bar]'})
    xlabel('Temperature [K]')
    
    subplot(numel(correlation_conductivity),3,3*i+2); contourf(T_check,P_check,squeeze(KT(:,:,2,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'kt');
    title('Conductivity 2-phase region:','the biggest root')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
    
    subplot(numel(correlation_conductivity),3,3*i+3); contourf(T_check,P_check,squeeze(KT(:,:,3,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'kt');
    title('Conductivity 2-phase region:','the smallest root')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
end

print(gcf, '-dpdf', '-fillpage', 'KT.pdf'); close all

%}

%% Plot polynomial

%{
set(gcf,'PaperOrientation','landscape', 'visible','off')
tiledlayout(numel(P_check), numel(T_check))
count = numel(P_check)-1;

for i = numel(P_check):-1:1

    for j = 1:numel(T_check)
        
        T = T_check(j);
        P = P_check(i);
        nexttile
        title([num2str(T-273.15),'[C],',num2str(P),'[bar]'],'FontSize',5)
        hold on
        scatter(Z(i,j,1),0,'k','o')
        scatter(Z(i,j,2),0,'k','o')
        scatter(Z(i,j,3),0,'k','o')
        yline(0)
        plot(N_polynomial, squeeze( pol(i,j,:) )); 
        ylim([-0.15, 0.45])
        hold off
        box off; axis square tight

    end
    count = count - 1;
end
print(gcf, '-dpdf', '-fillpage', 'Polynomials.pdf'); close all
%}

%% Test RBF - one value

%{
Data = Z(:,:,1);

testpoint_P = 5;
testpoint_T = 5;

P_RBF = P_check(1:testpoint_P:end)/Pc;
T_RBF = T_check(1:testpoint_T:end)/Tc;

MU_RBF = [1e-3, 1e-2, 2e-2, 5e-2, 1e-1, 5e-1, 1e0];
APP = nan(numel(P_check), numel(T_check));

[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [X(:), Y(:)];
distance = pdist2(centers,centers);

% Kernel Function
mu = MU_RBF(k);

GramMatrix = exp( - distance / mu );

% Weights
observations = Data(1:testpoint_P:end, 1:testpoint_T:end);
observations = reshape(observations.',[],1);

% Evaluate weights
%Weights = (GramMatrix'*GramMatrix)^(-1) * GramMatrix' * observations ;
Weights = GramMatrix \ observations ;

%% Evaluate the function at all the point

for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        APP(i,j) = sum(Weights.*exp(- pdist2([P T], centers) / mu)');
    end
end
%}

%% Test RBF - loop over multiple values of mu

%{
Data = Z(:,:,1);

testpoint_P = 5;
testpoint_T = 5;

P_RBF = P_check(1:testpoint_P:end)/Pc;
T_RBF = T_check(1:testpoint_T:end)/Tc;

MU_RBF = [1e-3, 1e-2, 2e-2, 5e-2, 1e-1, 5e-1, 1e0];
APP = nan(numel(P_check), numel(T_check), numel(MU_RBF));

[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [X(:), Y(:)];
distance = pdist2(centers,centers);

for k=1:numel(MU_RBF)

    if mod(k,round(numel(MU_RBF)/100)) == 0
        clc
        fprintf('%f %%\n',k/numel(MU_RBF)*100)
    end

    % Kernel Function
    mu = MU_RBF(k);

    GramMatrix = exp( - distance / mu );

    % Weights
    observations = Data(1:testpoint_P:end, 1:testpoint_T:end);
    observations = reshape(observations.',[],1);

    % Evaluate weights
    %Weights = (GramMatrix'*GramMatrix)^(-1) * GramMatrix' * observations ;
    Weights = GramMatrix \ observations ;

    %% Evaluate the function at all the point

    for i = 1:numel(P_check)

        for j = 1:numel(T_check)

            P = P_check(i)/Pc;
            T = T_check(j)/Tc;

            APP(i,j,k) = sum(Weights.*exp(- pdist2([P T], centers) / mu)');

        end

    end

    %max(max(abs(Data-APP)))

end

%%
for i=0:numel(MU_RBF)-1
    subplot(numel(MU_RBF),3,3*i+1);
    contourf(T_check,P_check,Data, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
    title('Original function')
    ylabel({['mu=',mat2str(MU_RBF(i+1))];'Pressure [bar]'})
    xlabel('Temperature [K]')

    subplot(numel(MU_RBF),3,3*i+2);
    hold on
    contourf(T_check,P_check,APP(:,:,i+1), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
    scatter(Y*Tc, X*Pc, 'k','SizeData',0.1)
    hold off
    title('RBF')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')

    subplot(numel(MU_RBF),3,3*i+3);
    hold on
    contourf(T_check,P_check,Data-APP(:,:,i+1), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
    scatter(Y*Tc, X*Pc, 'k','SizeData',0.1)
    hold off
    title('Difference')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
end
%}

%% Test RBF - optimze mu

%{
mu = MX.sym('mu');

Data = RHO(:,:,1);

testpoint_P = 10;
testpoint_T = 10;

P_RBF = P_check(1:testpoint_P:end)/Pc;
T_RBF = T_check(1:testpoint_T:end)/Tc;

APP = MX(numel(P_check), numel(T_check));

[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [X(:), Y(:)];
distance = pdist2(centers,centers);

% Kernel Function
GramMatrix = exp( - distance / mu );

% Weights
observations = Data(1:testpoint_P:end, 1:testpoint_T:end);
observations = reshape(observations.',[],1);

% Evaluate weights
disp(['Find the symbolic weights'])
%Weights = (GramMatrix'*GramMatrix)^(-1) * GramMatrix' * observations ;
tic
Weights = GramMatrix \ observations ;
toc

%% Evaluate the function at all the point

disp(['Create the matrix of data'])
for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        APP(i,j) = sum(Weights.*exp(- pdist2([P T], centers) / mu)');
    end
end

%% Set the optimzation problem to find mu
disp(['Set the structure of the NLP'])
nlp = struct;            % NLP declaration
nlp.x = mu;              % decision vars

D = (Data-APP).^2;
MSE = sum(D(:))/numel(Data);
%MSE = norm(Data-APP,'fro')^2/numel(Data);
nlp.f = MSE;               % objective - mean squared error

disp(['Create the solver'])
% Create solver instance
F = nlpsol('F','ipopt',nlp);

disp(['Solve NLP'])
% Solve the problem using a guess
tic
res = F('x0',[1e-4]);
toc

%% Evaluate opt mu and corresponding weights
APP_opt = nan(numel(P_check), numel(T_check));

mu = full(res.x);
disp(['The optimal value of mu is ',num2str(mu)])

GramMatrix = exp( - distance / mu );
tic
Weights = GramMatrix \ observations ;
toc

%% Approximate the dataset
for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        APP_opt(i,j) = sum(Weights.*exp(- pdist2([P T], centers) / mu)');
    end
end

%% Plot the results with opt mu and compare to the original data
subplot(1,3,1);
contourf(T_check,P_check, Data, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
title('Original function')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2);
hold on
contourf(T_check,P_check,APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(Y*Tc, X*Pc, 'k')
hold off
title('RBF')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3);
hold on
contourf(T_check,P_check,Data-APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(Y*Tc, X*Pc, 'k')
hold off
title('Difference')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

%}

%% Test RBF - optimze mu and weights

%{
mu = MX.sym('mu');

Data = RHO(:,:,1);

trainpoint_P = 5;
trainpoint_T = 5;

P_RBF = P_check(1:trainpoint_P:end)/Pc;
T_RBF = T_check(1:trainpoint_T:end)/Tc;

w = MX.sym('w', numel(P_RBF)*numel(T_RBF) );

APP = MX(numel(P_check), numel(T_check));

[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [X(:), Y(:)];
%distance = pdist2(centers,centers);

%% Evaluate the function at all the point

disp(['Create the matrix of data'])
for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        APP(i,j) = sum(w.*exp(- pdist2([P T], centers) / mu)');
    end
end

%% Set the optimzation problem to find mu
disp(['Set the structure of the NLP'])
nlp = struct;            % NLP declaration
nlp.x = [mu; w];              % decision vars

D = (Data-APP).^2;
MSE = sum(D(:))/numel(Data);
%MSE = norm(Data-APP,'fro')^2/numel(Data);
nlp.f = MSE;               % objective - mean squared error

disp(['Create the solver'])
% Create solver instance
F = nlpsol('F','ipopt',nlp);

disp(['Solve NLP'])
% Solve the problem using a guess
tic
res = F('x0',[1; ones(numel(w),1)],'ubx',inf ,'lbx',[ 0; -inf*ones(numel(w),1) ] );
toc

%% Extract the solution
solution = full(res.x);

mu_opt  = solution(1);
Weights_opt = solution(2:numel(w)+1);

%% Evaluate opt mu and corresponding weights

APP_opt = nan(numel(P_check), numel(T_check));

%mu = full(res.x);
%disp(['The optimal value of mu is ',num2str(mu)])

%% Approximate the dataset
for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        APP_opt(i,j) = sum(Weights_opt.*exp(- pdist2([P T], centers) / mu_opt)');
    end
end

%% Plot the results with opt mu and compare to the original data
subplot(1,3,1);
contourf(T_check,P_check, Data, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
title('Original function')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2);
hold on
contourf(T_check,P_check,APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(Y*Tc, X*Pc, 'k')
hold off
title('RBF')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3);
hold on
contourf(T_check,P_check,Data-APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(Y*Tc, X*Pc, 'k')
hold off
title('Difference')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

%}

%% Test RBF - optimze mu, weights and the location

Data = Z(:,:,1);

RBF = 15;

P_0 = P_check(1)/Pc + (P_check(end)/Pc-P_check(1)/Pc)*rand(RBF,1); 
T_0 = T_check(1)/Tc + (T_check(end)/Tc-T_check(1)/Tc)*rand(RBF,1);
[X_0,Y_0] = meshgrid(P_0,T_0);

P_max = P_check(end)*ones(RBF,1)/Pc;
T_max = T_check(end)*ones(RBF,1)/Tc;

P_min = P_check(1)*ones(RBF,1)/Pc;
T_min = T_check(1)*ones(RBF,1)/Tc;

mu    = MX.sym('mu', RBF);
P_RBF = MX.sym('P_RBF', RBF );
T_RBF = MX.sym('T_RBF', RBF );

w = MX.sym('w', RBF );

APP = MX(numel(P_check), numel(T_check));

%[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [P_RBF, T_RBF];

%% Evaluate the function at all the point

disp(['Create the matrix of data'])
for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        distance = sqrt(( centers(:,1)-P).^2 + (centers(:,2)-T).^2);
        APP(i,j) = sum( w.*exp(- distance ./ mu) );
    end
end

%% Set the optimzation problem to find mu

nlp = struct;            % NLP declaration
nlp.x = [mu; w; P_RBF; T_RBF];              % decision vars

D = (Data-APP).^2;
MSE = sum(D(:))/numel(Data);
%MSE = norm(Data-APP,'fro')^2/numel(Data);
nlp.f = MSE;               % objective - mean squared error

% Create solver instance
F = nlpsol('F','ipopt',nlp);


% Solve the problem 
tic
res = F('x0',[ones(numel(mu),1); ones(numel(w),1); P_0; T_0 ] ,'ubx',[inf*ones(numel(mu),1); inf*ones(numel(w),1); 1.05*P_max; 1.05*T_max] ,'lbx',[ zeros(numel(mu),1); -inf*ones(numel(w),1); 0.95*P_min; 0.95*T_min ] );
toc

%% Extract the solution
solution = full(res.x);

mu_opt      = solution(0*RBF+1:1*RBF);
Weights_opt = solution(1*RBF+1:2*RBF);
P_RBF_opt   = solution(2*RBF+1:3*RBF);
T_RBF_opt   = solution(3*RBF+1:4*RBF);

%[X_opt,Y_opt] = meshgrid(P_RBF_opt,T_RBF_opt);
centers_opt = [P_RBF_opt , T_RBF_opt ];

%% Evaluate opt mu and corresponding weights

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

hold off
title('RBF')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3);
hold on
contourf(T_check,P_check,Data-APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(T_RBF_opt*Tc, P_RBF_opt*Pc, 'k')
scatter(T_0*Tc, P_0*Pc, '*', 'k')
plot( [T_RBF_opt*Tc T_0*Tc ]' , [P_RBF_opt*Pc P_0*Pc]' ,'--k' )
hold off
title('Difference')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')


%%
%{
figure(2)
axis square tight

hold on
for i = 1:RBF
    mu = centers_opt(i,:);
    Sigma =  [mu_opt(1) 1 ; 1 mu_opt(1)];
    %[X,Y] = meshgrid( T_min(1):.1:T_max(1) , P_min(1):.1:P_max(1) );
    [X,Y] = meshgrid( -10:.1:10 , -10:.1:10 );
    p = mvnpdf([X(:) Y(:)],mu,Sigma);
    p = reshape(p,size(X));
    contour(X,Y,p)
    scatter(T_RBF_opt(i), P_RBF_opt(i), 'k')
end
yline(1)
xline(1)
hold off

%}