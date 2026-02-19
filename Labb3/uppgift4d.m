N = 200; 
T = 2;
dx = 1/N;
dt = dx/2.0; 
M = round(T/dt);
c = 1; 

u = zeros(N+1, M+1);
p = zeros(N+1, M+1);

main_diag = -2 * ones(N+1, 1);
off_diag = ones(N, 1);
A_N = (diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1)) / dx^2;

A_N(1, 2) = 2 / dx^2;   
A_N(N+1, N) = 2 / dx^2; 

x = dx * (0:N)';

% Begynnelsepuls fr fr
g = @(x) exp(-200*(x-0.5).^2);
u(:,1) = g(x);
p(:,1) = zeros(N+1, 1); 

E = zeros(1, M+1);

weights = ones(N+1, 1);
weights(1) = 0.5;
weights(end) = 0.5;
W = diag(weights);

E(1) = 0.5 * (p(:,1)' * W * p(:,1)) - 0.5 * c^2 * (u(:,1)' * W * A_N * u(:,1));

figure('Name', 'Vågpropagering (Neumann)');
hPlot = plot(x, u(:,1), 'r', 'Linewidth', 1.5);
axis([0 1 -1 1]);
title('Tidsutveckling av vågekvationen (Fria ändar)');
xlabel('Position x'); 
ylabel('Amplitud u');
drawnow;

% Symplektisk Euler 
for m = 1:M
    p(:,m+1) = p(:,m) + dt * c^2 * A_N * u(:,m);
    
    u(:,m+1) = u(:,m) + dt * p(:,m+1);
    
    E(m+1) = 0.5 * (p(:,m+1)' * W * p(:,m+1)) - 0.5 * c^2 * (u(:,m+1)' * W * A_N * u(:,m+1));

    set(hPlot, 'YData', u(:,m+1));
    drawnow;
end

figure('Name', 'Energins tidsutveckling');

t_vec = linspace(0, T, M+1); 

plot(t_vec, E, 'k', 'LineWidth', 1.5);
grid on;

title('Total energi över tid');
xlabel('Tid t');
ylabel('Energi E');
