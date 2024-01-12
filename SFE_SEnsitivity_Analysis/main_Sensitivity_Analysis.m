startup;
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

%% Load paramters
m_total                 = 80;

V_Flow                  = 0.4;                                              % Volumetric flow rate l/min

% Bed geometry
before                  = 0.1;                                              % Precentage of length before which is empty
bed                     = 0.165;                                            % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 2000;
timeStep                = 0.5;                                              % Minutes
SamplingTime            = 5;                                                % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);

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
T0homog                 = 45+273.15;
feedPress               = 250;

Z                       = Compressibility( T0homog, feedPress,         Parameters );
rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );

enthalpy_rho            = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );

feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;  % Kelvin

feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;  % Bars
%feedPress(100:200)     = 300;

feedFlow                = V_Flow * rho * 1e-3 / 60 * ones(1,length(Time_in_sec));  % l/min -> kg/s

uu                      = [feedTemp', feedPress', feedFlow'];

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

% Initial conditions
x0                      = [ C0fluid'                        ;
                            C0solid  * bed_mask             ;
                            enthalpy_rho * ones(nstages,1)  ;
                            feedPress(1)                    ;
                            0                               ;
                            ];

%% Run plain simulation
Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];
[xx_0]                 = simulateSystem(F, [], x0, Parameters_init_time );

%% Set sensitivity analysis

name_v = {'T_{in}', 'P', 'F'};

My_Font = 24;
My_Font_cont = 8;
num_levels = 100;

for ii = 1

        %% Sensitivities calculations
        Parameters{8} = ii;
        Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];
        [S,p,Sdot]              = Sensitivity(x, xdot, u, ii );

        % Initial conditions
        x0_SA                   = [ x0; zeros(length(S)-length(xdot),1) ];

        f_SA = @(S, p) Sdot(S, p, bed_mask);
        Results = Integrator_SS(Time*60, x0_SA, S, p, Sdot, Parameters_init_time);
        Res = Results(Nx+1:end,:);

        %% Sensitivities plot - normal scale
        %subplot(3,2,1)
        imagesc(Time, L_nstages, Res(1*nstages+1:2*nstages,:)); cb = colorbar;
        ax=gca; ax.FontSize = My_Font;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'k--');
        plot([ExtractionTime-40, ExtractionTime-40], [L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))], 'k--');
        text(ExtractionTime-45,[mean([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))])],'fixed bed', 'FontSize', My_Font-8, 'Interpreter', 'latex', 'Color', 'black', 'HorizontalAlignment','center','VerticalAlignment','middle', 'Rotation', 90);
        hold off
        cb.Label.String = ['$\frac{dc_s}{d',name_v{ii},'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 1.5*My_Font;
        cb.Label.Rotation = 0; % to rotate the text
        xlabel('Time [min]'); ylabel('Length [m]'); 
        %title(sprintf('$\\rho_f$=%4.2f',rho))
        pbaspect([1 1 1])
        
        set(gcf,'PaperOrientation','landscape'); print(figure(1),['CS_',name_v{ii},'.pdf'],'-dpdf','-bestfit')
        close all;

        %subplot(3,2,2)
        imagesc(Time, L_nstages, Res(2*nstages+1:3*nstages,:) ); cb = colorbar;
        ax=gca; ax.FontSize = My_Font;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'k--');
        plot([ExtractionTime-40, ExtractionTime-40], [L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))], 'k--');
        text(ExtractionTime-45,[mean([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))])],'fixed bed', 'FontSize', My_Font-8, 'Interpreter', 'latex', 'Color', 'black', 'HorizontalAlignment','center','VerticalAlignment','middle', 'Rotation', 90);
        hold off
        cb.Label.String = ['$\frac{d(h\times\rho)}{d',name_v{ii},'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 1.5*My_Font;
        cb.Label.Rotation = 0; % to rotate the text
        xlabel('Time [min]'); ylabel('Length [m]')
        pbaspect([1 1 1])

        set(gcf,'PaperOrientation','landscape'); print(figure(1),['H_',name_v{ii},'.pdf'],'-dpdf','-bestfit')
        close all;

        %subplot(3,2,3)
        imagesc(Time, L_nstages, Res(0*nstages+1:1*nstages,:)); cb = colorbar;
        ax=gca; ax.FontSize = My_Font;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'k--');
        plot([ExtractionTime-40, ExtractionTime-40], [L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))], 'k--');
        text(ExtractionTime-45,[mean([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))])],'fixed bed', 'FontSize', My_Font-8, 'Interpreter', 'latex', 'Color', 'black', 'HorizontalAlignment','center','VerticalAlignment','middle', 'Rotation', 90);
        hold off
        cb.Label.String = ['$\frac{dc_f}{d',name_v{ii},'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 1.5*My_Font;
        cb.Label.Rotation = 0; % to rotate the text
        xlabel('Time [min]'); ylabel('Length [m]')
        pbaspect([1 1 1])

        set(gcf,'PaperOrientation','landscape'); print(figure(1),['CF_',name_v{ii},'.pdf'],'-dpdf','-bestfit')
        close all;

        %% Sensitivities plot - log scale

        [C,h] = contourf(log10(Time+1), L_nstages, Res(1*nstages+1:2*nstages,:), num_levels, 'EdgeColor', 'none'); view(0,-90); cb = colorbar;
        clabel(C,h,'FontSize',My_Font_cont)
        ax=gca; ax.FontSize = My_Font;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'k--');
        plot(log10([ExtractionTime-40, ExtractionTime-40]), [L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))], 'k--');
        text(log10(ExtractionTime-45),[mean([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))])],'fixed bed', 'FontSize', My_Font-8, 'Interpreter', 'latex', 'Color', 'black', 'HorizontalAlignment','center','VerticalAlignment','middle', 'Rotation', 90);
        hold off
        cb.Label.String = ['$\frac{dc_s}{d',name_v{ii},'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 1.5*My_Font;
        cb.Label.Rotation = 0; % to rotate the text
        xlabel('log_{10}(Time) [min]'); ylabel('Length [m]'); 
        %title(sprintf('$\\rho_f$=%4.2f',rho))
        pbaspect([1 1 1])
        axis tight
        
        set(gcf,'PaperOrientation','landscape'); print(figure(1),['CS_',name_v{ii},'_log.pdf'],'-dpdf','-bestfit')
        close all;

        %subplot(3,2,2)
        [C,h] = contourf(log10(Time+1), L_nstages, Res(2*nstages+1:3*nstages,:), num_levels, 'EdgeColor', 'none'); view(0,-90); cb = colorbar;
        clabel(C,h,'FontSize',My_Font_cont)
        ax=gca; ax.FontSize = My_Font;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'k--');
        plot(log10([ExtractionTime-40, ExtractionTime-40]), [L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))], 'k--');
        text(log10(ExtractionTime-45),[mean([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))])],'fixed bed', 'FontSize', My_Font-8, 'Interpreter', 'latex', 'Color', 'black', 'HorizontalAlignment','center','VerticalAlignment','middle', 'Rotation', 90);
        hold off
        cb.Label.String = ['$\frac{d(h\times\rho)}{d',name_v{ii},'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 1.5*My_Font;
        cb.Label.Rotation = 0; % to rotate the text
        xlabel('log_{10}(Time) [min]'); ylabel('Length [m]')
        pbaspect([1 1 1])
        axis tight

        set(gcf,'PaperOrientation','landscape'); print(figure(1),['H_',name_v{ii},'_log.pdf'],'-dpdf','-bestfit')
        close all;

        %subplot(3,2,3)
        [C,h] = contourf(log10(Time+1), L_nstages, Res(0*nstages+1:1*nstages,:), num_levels, 'EdgeColor', 'none'); view(0,-90); cb = colorbar;
        clabel(C,h,'FontSize',My_Font_cont);
        ax=gca; ax.FontSize = My_Font;
        hold on
        yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'k--');
        plot(log10([ExtractionTime-40, ExtractionTime-40]), [L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))], 'k--');
        text(log10(ExtractionTime-45),[mean([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))])],'fixed bed', 'FontSize', My_Font-8, 'Interpreter', 'latex', 'Color', 'black', 'HorizontalAlignment','center','VerticalAlignment','middle', 'Rotation', 90);
        hold off
        cb.Label.String = ['$\frac{dc_f}{d',name_v{ii},'}$']; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 1.5*My_Font;
        cb.Label.Rotation = 0; % to rotate the text
        xlabel('log_{10}(Time) [min]'); ylabel('Length [m]');    
        pbaspect([1 1 1])
        axis tight

        set(gcf,'PaperOrientation','landscape'); print(figure(1),['CF_',name_v{ii},'_log.pdf'],'-dpdf','-bestfit')
        close all;

        %% Sensitivities plot - P and y

        %subplot(3,2,4)
        plot(Time, Res(end-1,:), LineWidth=2); cb = colorbar;
        ax=gca; ax.FontSize = My_Font;
        xlabel('Time [min]'); ylabel(['$\frac{d P_{t-1}}{d',name_v{ii},'}$'])
        pbaspect([1 1 1])

        set(gcf,'PaperOrientation','landscape'); print(figure(1),['P_',name_v{ii},'.pdf'],'-dpdf','-bestfit')
        close all;

        %subplot(3,2,5)
        hold on
        plot(Time, Res(end,:), LineWidth=2); colorbar;
        yline(0, LineWidth=2)
        xlabel('Time [min]'); ylabel(['$\frac{d y}{d',name_v{ii},'}$'])
        hold off
        ax=gca; ax.FontSize = My_Font;
        colormap jet
        pbaspect([1 1 1])

        set(gcf,'PaperOrientation','landscape'); print(figure(1),['Y_',name_v{ii},'.pdf'],'-dpdf','-bestfit')
        close all;

        %U = [U; Res(end,:)];
end
%}