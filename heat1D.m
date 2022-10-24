% Finite Difference Solver for the 1D Heat Equation
% Unknown: u(x,t)
% PDE: d_t u = k d_x^2 u, 0 < x < L, t > 0
% BC: u(0,t) = 0, u(1,t) = 0 (Dirichlet)
% IC: u(x,0) = f(x)

% setup
clear;
showplot = 1;

L = 1; % length of domain
Nt = 201; % number of time steps
Nx = 31; % number of spatial steps
Tf = 0.5; % stopping time
dt = Tf/Nt; % time step
dx = L/Nx; % spatial step
f = @(x) (sin(pi*x/L)).^8; % initial condition
k = 1/2*(dx)^2/dt; % thermal diffusivity (chosen for stability)
alpha = k*dt/dx^2;

xgrid = linspace(0,L,Nx);
tgrid = linspace(0,Tf,Nt);

u = nan(Nx,Nt);
u(:,1) = f(xgrid); % initial condition

% PDE solver
for t = 2:Nt
    u(2:end-1,t) = u(2:end-1,t-1) + alpha*( u(3:end,t-1)...
        -2*u(2:end-1,t-1)+u(1:end-2,t-1) ); % forward difference
    u(1,t) = 0; % BC1
    u(end,t) = 0; % BC2
end

% plot
if showplot
    figure(1)
    for t = 1:Nt
        plot(xgrid,u(:,t));
        xlim([0 L])
        ylim([0 1])
        pause(0.025);
    end
end
