clear all; clc; close all;
%% Initialization

% Hyperparameters
alpha = 1.0; % unused
delta = 1;
theta = 15;
a = sind(theta);
b = -delta/3*cosd(theta);

% Time
dt = 0.1;
time = 120; % time max
epoch_time = floor(time/dt);
time_step = 2; % graph changing every mod(i,time_step)==0

% Grid 1D
Nx = 1001; % N-many elements in x vector
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
% 1
sa = dt*a/(2*dx)/2;
% sa = dt*a/(2*dx); % implicit du/dx
sb = dt*b/dx^2/2;

% 2
% sa = dt*a/(2*dx);
% sb = dt*b/dx^2;

A = eye([Nx,Nx]);
A2 = A; % An+1
A1 = A; % An
for i = 2:Nx-1
    for j = 2:Nx-1
        if A(i,j) == 1
            % 1
            A1(i,j) = -1 - 2*sb;
            A1(i,j-1) = -sa + sb;
            A1(i,j+1) = sa + sb;

            A2(i,j) = -1 + 2*sb;
            A2(i,j-1) = sa - sb;
            A2(i,j+1) = -sa - sb;

            % 2
%             A1(i,j) = -2 - 2*sb;
%             A1(i,j-1) = -sa + sb;
%             A1(i,j+1) = sa + sb;
% 
%             A2(i,j) = -2 + 2*sb;
%             A2(i,j-1) = sa - sb;
%             A2(i,j+1) = -sa - sb;
            
            % is equal to 1
%             A1(i,j) = 1 + 2*sb;
%             A1(i,j-1) = sa - sb;
%             A1(i,j+1) = -sa - sb;
% 
%             A2(i,j) = 1 - 2*sb;
%             A2(i,j-1) = -sa + sb;
%             A2(i,j+1) = sa + sb;
            
            % Implicit du/dx
%             A1(i,j) = 1 + 2*sb;
%             A1(i,j-1) = sa - sb;
%             A1(i,j+1) = -sa - sb;
% 
%             A2(i,j) = 1 - 2*sb;
%             A2(i,j-1) = sb;
%             A2(i,j+1) = sb;
        end
    end
end

% t0fd = tic; % stopwatch start
% for i = 1:epoch_time
%     Ufd = (A1*Ufd.').'/A2;
%     Ufd(1) = 0; Ufd(end) = 0; % boundary conditions
%     hUfd(i,:) = Ufd;
% end
% tfd = toc(t0fd); % stopwatch end

t0fd = tic; % stopwatch start
for i = 1:epoch_time
    Ufd(1,2:Nx-1) = A2(2:Nx-1,2:Nx-1)\(A1(2:Nx-1,2:Nx-1)*Ufd(1,2:Nx-1)');
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