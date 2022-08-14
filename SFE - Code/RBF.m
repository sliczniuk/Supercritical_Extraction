function [] = RBF()
    %% import libraries
    %addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
    addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
    import casadi.*

    %% define data - toy example
    x = linspace(0,1,14);
    func = @(x) exp(x .* cos(3*pi*x) );
    y = func(x);

    %% choose a kernel
    phi = @(r) exp(-(3*r).^2);
    
    %% define phi matrix and fill it with data
    phi_matrix = nan(numel(x));

    for i=1:numel(x)
        for j=1:numel(x)
            % find the distance between each point
            r = pdist([x(i); x(j)],'euclidean');
            % fill the matrix with phi
            phi_matrix(i,j) = phi(r);
        end
    end
    phi_matrix
    %% solve the matrix equation Ax=b to find weights
    w = phi_matrix / y

    %% evaluate the function
    s = phi_matrix * w;
    
    %%
    hold on
    plot(x,y)
    plot(x,s,'o')
    hold off

end