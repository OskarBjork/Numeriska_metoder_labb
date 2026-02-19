N = 200; 
T = 2;
dx = 1/N;
dt = dx/2.0; 
M = round(T/dt);
c = 1; 

u = zeros(N-1, M+1);
p = zeros(N-1, M+1);

main_diag = -2 * ones(N-1, 1);
off_diag = ones(N-2, 1);
A = (diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1)) / dx^2;

x = dx * (1:N-1)';

% Begynnelsepuls fr fr
g = @(x) exp(-200*(x-0.5).^2);
u(:,1) = g(x);
p(:,1) = zeros(N-1, 1); 

E = zeros(1, M+1);
E(1) = 0.5 * norm(p(:,1))^2 - 0.5 * c^2 * (u(:,1)' * A * u(:,1));

figure('Name', 'Vågpropagering (Dirichlet)');
X_full = [0; x; 1]; 
U_full = [0; u(:,1); 0]; 
hPlot = plot(X_full, U_full, 'b', 'Linewidth', 1.5);
axis([0 1 -1 1]);
title('Tidsutveckling av vågekvationen');
xlabel('Position x'); 
ylabel('Amplitud u');
drawnow;


% Symplektisk Euler
for m = 1:M
    p(:,m+1) = p(:,m) + dt * c^2 * A * u(:,m);
    
    u(:,m+1) = u(:,m) + dt * p(:,m+1);
    
    E(m+1) = 0.5 * norm(p(:,m+1))^2 - 0.5 * c^2 * (u(:,m+1)' * A * u(:,m+1));
    
    U_full = [0; u(:,m+1); 0];
    set(hPlot, 'YData', U_full);
    drawnow;
end

figure('Name', 'Energins tidsutveckling');
t_vec = linspace(0, T, M+1);
plot(t_vec, E, 'k', 'LineWidth', 1.5);
