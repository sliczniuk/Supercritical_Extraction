function [xout, kout] = singleShooting_ParameterEstimation(OCP, x0, u, k0,  parameters)

%addpath('C:\Users\slicz\Desktop\Aalto\Matlab\SFE\casadi-windows-matlabR2016a-v3.5.1');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

nstages = parameters{1};

OCP_solver = casadi.Opti();

% http://www.diva-portal.se/smash/get/diva2:956377/FULLTEXT01.pdf
nlp_opts = struct;
nlp_opts.ipopt.max_iter = 100;
%nlp_opts.ipopt.acceptable_iter = 50;
%nlp_opts.ipopt.acceptable_tol = 1e-6;
%nlp_opts.ipopt.tol = 1e-7;

ocp_opts = {'nlp_opts', nlp_opts};
OCP_solver.solver('ipopt',nlp_opts)

K = OCP_solver.variable(OCP.Nk);
X = [MX(x0) zeros(OCP.Nx,OCP.N)];
J = 0;

yout = MX(zeros(OCP.Ny,OCP.N+1));   

%g = @(x, u, y_old) modelSFE_out2(x, u, y_old, parameters, timeStep_in_sec);

for j=1:OCP.N
    %J=J+OCP.L(U(:,j));
    X(:,j+1)=OCP.F(X(:,j),u(:,j),K);
    %yout(:,j+1) = g(X(:,j+1), u(:,j), yout(:,j));
end

%data = yy(1:N_Sample:end);
data = X(3*nstages+1,:);
data = [zeros(1,OCP.N_Delay) data(1:end-OCP.N_Delay)];
data = data(1:OCP.N_Sample:end);

%J = J + OCP.LS(data);
J = J + OCP.MSE(data,K(end));
%J = J + OCP.MAP(data,K(1),K(end));

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

OCP_solver.minimize(J);
OCP_solver.set_initial(K, k0);

try
    sol = OCP_solver.solve();
    kout = sol.value(K);
catch
    kout = OCP_solver.debug.value(K);
end

xout = 1;
end