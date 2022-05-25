function [yout,tout,xout] = simulateSystem_ParameterEstimation(varargin)

    % [y, t, x] = Simulate(F, g, x0, u); 
    
    import casadi.*
    
    F = varargin{1}; 
    g = varargin{2}; 
    x0 = varargin{3};
    
    u = varargin{4};
    N = size(u,1);
        
    k = varargin{5};

    f = @(x,u,k) full(F(x,u,k));
    
    % aux
    Nx = size(x0,1); 
    Ny = size(g(x0),1);
    
    tout = 1:N;
    
    xout = zeros(Nx,N+1);  
    xout(:,1) = x0;
    
    yout = zeros(Ny,N+1);   
    yout(:,1) = g(xout(:,1));
    
    % sim
    for i = 1:N
        xout(:,i+1) = f(xout(:,i), u(i,:),k);
        yout(:,i+1) = g(xout(:,i+1));
        
    end
    
    % extract the initial state from the vector that is returned
    %xout = xout(:,2:end); yout = yout(:,2:end);
end
