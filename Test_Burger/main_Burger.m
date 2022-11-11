clc, clear all

Nx = 100;
L = linspace(0,1,Nx);

dz = L(2) - L(1);

u0 = zeros(Nx,1);
u0(45:55) = 1;

tspan = [0 1];

a = -0.1;

[t1,u1] = ode45(@(t,u) burger_Gudanov(t, u, Nx, dz, a), tspan, u0);
[t2,u2] = ode45(@(t,u) burger_Backward(t, u, Nx, dz, a), tspan, u0);

subplot(2,1,1)
imagesc(t1,L,u1); colorbar; colormap jet
subplot(2,1,2)
imagesc(t2,L,u2); colorbar; colormap jet