% Vågekvationen: Implementation med Dirichlet-randvillkor (Uppgift 4c)
N = 200; 
T = 2;
dx = 1/N;
dt = dx/2.0; % CFL-stabilitetsvillkoret uppfyllt (0.5 <= 1)
M = round(T/dt);
c = 1; % Vågens fortplantningshastighet

% Allokering av lösningsmatriser u (rum, tid) och hastigheter p
u = zeros(N-1, M+1);
p = zeros(N-1, M+1);

% Uppbyggnad av systemmatrisen A (Andraderivatan med finita differenser)
main_diag = -2 * ones(N-1, 1);
off_diag = ones(N-2, 1);
A = (diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1)) / dx^2;

% Positionsvektor (enbart inre punkter evalueras)
x = dx * (1:N-1)';

% Applicering av begynnelsedata (en gaussisk puls enligt Uppgift 4c)
g = @(x) exp(-200*(x-0.5).^2);
u(:,1) = g(x);
p(:,1) = zeros(N-1, 1); % Strängen hålls i stillhet innan den släpps

% Energivektor för diagnostik
E = zeros(1, M+1);
E(1) = 0.5 * norm(p(:,1))^2 - 0.5 * c^2 * dx * (u(:,1)' * (A * u(:,1)));

% --- Förberedelse av visualisering ---
figure('Name', 'Vågpropagering (Dirichlet)');
X_full = [0; x; 1]; % Addera rändernas positioner
U_full = [0; u(:,1); 0]; % Addera nollorna från ränderna för plotten
hPlot = plot(X_full, U_full, 'b', 'Linewidth', 1.5);
axis([0 1 -1 1]);
title('Tidsutveckling av vågekvationen');
xlabel('Position x'); 
ylabel('Amplitud u');
drawnow;

% --- Tidsstegning: Symplektisk Euler ---
for m = 1:M
    % Steg 1: Uppdatera hastigheten (p) explicit med hjälp av aktuella positionen (u)
    p(:,m+1) = p(:,m) + dt * c^2 * A * u(:,m);
    
    % Steg 2: Uppdatera positionen (u) explicit via den NYA hastigheten
    u(:,m+1) = u(:,m) + dt * p(:,m+1);
    
    % Beräkna momentan energi i varje steg
    E(m+1) = 0.5 * norm(p(:,m+1))^2 - 0.5 * c^2 * dx * (u(:,m+1)' * (A * u(:,m+1)));
    
    % Uppdatera den grafiska animationen
    U_full = [0; u(:,m+1); 0];
    set(hPlot, 'YData', U_full);
    drawnow;
end