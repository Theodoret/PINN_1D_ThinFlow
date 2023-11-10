clear all; clc; close all;
%% Initialization

% Hyperparameters
delta = 1;
theta = 15;
alpha = sind(theta);
beta = -delta/3*cosd(theta);

% Time
dt = 0.1;
time = 120; % time max
epoch_time = floor(time/dt);
time_step = 100; % graph changing every mod(i,time_step)==0

% Grid 1D
Nx = 101; % N-many elements in x vector
xmin = 0; xmax = 100; % boundary x vector
x = linspace(xmin,xmax,Nx); % Grid
dx = x(2) - x(1);

% Initial values
U = zeros(1,Nx);
for i = 1:Nx
    if x(i) >= 10 && x(i) <= 20
        U(i) = sin(0.1*pi*(x(i)-10));
    end
end

% History
hUfd = zeros(epoch_time,Nx);



%% Other solutions
% Numerical solution - Finite Difference Method
Ufd = U;
sa = dt*alpha/(2*dx);
sb = dt*beta/dx^2;

% A = eye([Nx,Nx]);
% for i = 2:Nx-1
%     for j = 2:Nx-1
%         if A(i,j) == 1
%             A(i,j) = 1 - 2*sb;
%             A(i,j-1) = -sa + sb;
%             A(i,j+1) = sa + sb;
%         end
%     end
% end

middle = (1 - 2*sb)*ones(1,Nx-2);
left = (-sa + sb)*ones(1,Nx-3);
right = (sa + sb)*ones(1,Nx-3);
A = diag(middle) + diag(left,-1) + diag(right,1);

t0fd = tic; % stopwatch start
for i = 1:epoch_time
    Ufd(1,2:Nx-1) = A\Ufd(1,2:Nx-1)';
%     Ufd(1) = 0; Ufd(end) = 0; % boundary conditions
    hUfd(i,:) = Ufd;
end
tfd = toc(t0fd); % stopwatch end



%% Results and Visualisations
% Visualisations
% Final result
figure(1);
plot(x,Ufd,'g--'); % Finite Difference solution
title("Results") ; xlabel ("x"); ylabel ("y","Rotation",0);

figure(2);
nexttile; surf(x,1:epoch_time,hUfd); shading interp; colorbar;
% view(0,90);