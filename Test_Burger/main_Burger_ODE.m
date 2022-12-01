clc; clear all;

%% Selection of numerical parameters (time step, grid spacing, etc...).
xend = 2;   % x-axis size.
tend = 0.5;   % t-axis size.

N = 1000;
dx = xend/N; % Grid spacing

dt = 0.001;
x  = 0:dx:xend;

nt = floor(tend/dt);
dt = tend / nt;

%% Set up the initial solution values.

xshift = round(N/2);

uL = 1;
uR = 2;

u0 = ones(numel(x),1);
u0(1:xshift) = uL;
u0(xshift+1:end) = uR;

%u0 = exp(-2*(x-2).^2)';

%{

%% Selection of equation and method.

%'Up-wind nonconservative','Up-wind conservative','Lax-Friedrichs','Lax-Wendroff','MacCormack','Godunov'
method = 2;

%'Piecewise constant (shock)','Piecewise constant (expansion)','Gaussian','Piecewise continuous'
ictype = 1;

u0  = uinit(x,ictype);  %Call to the function "uinit".


u   = u0;
unew = 0*u;

%% Implementation of the numerical methods.
    for i = 1:nt
        switch method
            case 1  %Up-wind nonconservative
                
                unew = u - dt/dx * u .* (u - [u(1), u(1:end-1) ]);
                
            case 2  %Up-wind conservstive
                unew = u - dt/dx * (f(u) - f( [u(1), u(1:end-1)] ) );
                
            case 3  %Lax-Friedrichs
                unew(2:end-1) = 0.5*(u(3:end)+u(1:end-2)) - 0.5*dt/dx * ...
                    (f(u(3:end)) - f(u(1:end-2)));
                unew(1)   = u(1);
                unew(end) = u(end);

            case 4  %Lax-Wendroff
                unew(2:end-1) = u(2:end-1) ...
                    - 0.5*dt/dx * (f(u(3:end)) - f(u(1:end-2))) ...
                    + 0.5*(dt/dx)^2 * ...
                    ( df(0.5*(u(3:end) + u(2:end-1))) .* (f(u(3:end)) - f(u(2:end-1))) - ...
                    df(0.5*(u(2:end-1) + u(1:end-2))) .* (f(u(2:end-1)) - f(u(1:end-2))) );
                unew(1)   = u(1);
                unew(end) = u(end);    

            case 5  %MacCormack
                us    = u(1:end-1) - dt/dx * (f(u(2:end)) - f(u(1:end-1)));
                unew(2:end-1)= 0.5*(u(2:end-1) + us(2:end)) - ...
                    0.5*dt/dx* (f(us(2:end)) - f(us(1:end-1)));
                unew(1)   = u(1);
                unew(end) = u(end);

            case 6  %Godunov
                unew(2:end-1) =u(2:end-1)- dt/dx*(nf(u(2:end-1),u(3:end)) - nf(u(1:end-2),u(2:end-1)));
                unew(1)   = u(1);
                unew(end) = u(end);  

        end
        u = unew;
        U(i,:) = u(:);
    end

U=[u0;U];
%}

T=0:dt:tend;

%%
[tNConODE45,yNConODE45] = ode45(@(t,y) Burger_NonConservative(t,y,dx), T, u0);
[tConODE45,yConODE45]   = ode45(@(t,y) Burger_Conservative(t,y,dx), T, u0);

yConE  = Burger_E_Conservative(T,u0,dx,dt);
yNConE = Burger_E_NonConservative(T,u0,dx,dt);

[tNConODE45_CS,yNConODE45_CS] = ode45(@(t,y) Burger_NonConservative_CS(t,y,dx), T, u0);

yConLF = Burger_E_LF(T,u0,dx,dt);
[tConLF_ODE45,yConLF_ODE45] = ode45(@(t,y) Burger_NonConservative_CS(t,y,dx), T, u0);

%% Plot of the solutions.
f = figure(1);
for i =1:10:numel(T)
    
    sgtitle(num2str(tNConODE45(i),'%.2f'));

    subplot(2,3,1)
    plot(x,yNConODE45(i,:))
    title('Non-conservative, ODE45')
    ylim([0, max(yNConODE45(:))])

    subplot(2,3,2)
    plot(x,yConODE45(i,:))
    title('Conservative, ODE45')
    ylim([0, max(yConODE45(:))])

    subplot(2,3,4)
    plot(x,yNConE(i,:))
    title('Non-conservative, EULER')
    ylim([0, max(yNConE(:))])

    subplot(2,3,5)
    plot(x,yConE(i,:))
    title('Conservative, EULER')
    ylim([0, max(yConE(:))])

    subplot(2,3,3)
    plot(x,yConLF_ODE45(i,:))
    title('Non-Conservative, Lax-Friedrichs ODE45')
    ylim([0, max(yConLF_ODE45(:))])   

    subplot(2,3,6)
    plot(x,yConLF(i,:))
    title('Non-Conservative, Lax-Friedrichs EULER')
    ylim([0, max(yConLF(:))])

    drawnow
end