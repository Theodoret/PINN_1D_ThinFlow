clear all; clc; close all;
theta =  15;
delta = 1;
alpha = sind(theta);
beta = -delta/3 * cosd(theta);

Nx = 1;
Ny = 101;

x = linspace(0,100,Ny);
dy = x(2)-x(1); % space step

T = 120; % total time
dt = 0.1; % time step
n = floor(T/dt); % number of iterations

% Analytical solution
% P1D = sin(pi.*x).*exp(1).^(-pi^2.*T);

U = zeros(Nx,Ny);
for i = 1:Ny
    if x(i) >= 10 && x(i) <= 20
        U(1,i) = sin(0.1*pi*(x(i)-10));
    end
end
Ud = U;

% Initialiaze K1 and K2
K1U = zeros(Nx,Ny);
K2U = zeros(Nx,Ny);

history = zeros(n,Ny);

step_plot = 100;
t = 0;

t0iter = tic;

for i = 1:n
    
    [K1U] = RK2(U,dy,Nx,Ny,alpha,beta);
    Ud = U + dt *K1U;
 
    % Neumann conditions: derivatives at edges are null.
    [Ud] = neumann(Ud,Ny);

    
    [K2U] = RK2(Ud,dy,Nx,Ny,alpha,beta);
    U = U + (dt/2) * (K1U+K2U);
   
    % Neumann conditions: derivatives at edges are null.
    [U] = neumann(U,Ny);
    t = t +dt;
    
    history(i,:) = U;
 
    % We plot the state of the system at 9 different times.
    if mod(i,step_plot)==0 %&& i < 9 * step_plot
        minU = min(U(:));  maxU = max(U(:));
         fprintf(2, 'i = %5d ; time = %7.4f ; min U = %7.4f ; max U = %7.4f\n', i, t, minU, maxU);
        
%         figure(i)
        figure(gcf)
        colormap jet;
        plot(x,U);
    end
end

titer = toc(t0iter);

figure(1);plot(x,U,'r-.');
% hold on;
% plot(x,P1D) % Exact solution
% lh = legend('PINN prediction','Exact solution'); 
% set(lh,'interpreter','latex','fontsize',16,'location','northwest');
% set(gca,'fontsize',16)
% hold off;

figure(2);
s = surf(x,(1:n)*dt,history); s.EdgeColor = "interp";
xlabel("x"); ylabel("t"); zlabel("eta");

% error = mean((P1D-U).^2);



function [KU] = RK2(U,dx,Nx,Ny,alpha,beta)
    % We compute the Laplacian of u and v.
    deltaU2 = laplacian2(U,dx,Nx,Ny);
    deltaU1 = laplacian1(U,dx,Nx,Ny);
    % We look for K1 & K2 value
    KU = -(alpha.*(deltaU1) + beta.*(deltaU2));
end

function [U] = neumann(U,Ny)
    U(1,1)  = 0;
    U(1,Ny) = 0;
end

function L = laplacian1(Z,h,Nx,Ny)
    L = Z;
    for j = 2:Ny-1
        L(1,j) = (Z(Nx,j+1)-Z(Nx,j-1))/(2*h);
    end
end


function L = laplacian2(Z,h,Nx,Ny)
    L = Z;
    for j = 2:Ny-1
        L(1,j) = (Z(Nx,j+1)-2*Z(Nx,j)+Z(Nx,j-1))/h^2;
    end
end