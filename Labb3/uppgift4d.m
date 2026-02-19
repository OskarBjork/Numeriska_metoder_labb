% Vågekvationen: Implementation med Neumann-randvillkor (Uppgift 4d)
N = 200; 
T = 2;
dx = 1/N;
dt = dx/2.0; % CFL-stabilitetsvillkoret
M = round(T/dt);
c = 1; % Vågens fortplantningshastighet

% 1. Allokering av lösningsmatriser u och p
% Dimensionen är nu (N+1) eftersom ränderna vid j=0 och j=N är frihetsgrader
u = zeros(N+1, M+1);
p = zeros(N+1, M+1);

% 2. Konstruktion av utökad matris A_N för Neumann-randvillkor, dim (N+1)x(N+1)
main_diag = -2 * ones(N+1, 1);
off_diag = ones(N, 1);
A_N = (diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1)) / dx^2;

% Modifiering av koefficienterna på grund av reflektion av spökpunkter (Eq 13)
A_N(1, 2) = 2 / dx^2;   % j=0 ser u_1 två gånger
A_N(N+1, N) = 2 / dx^2; % j=N ser u_{N-1} två gånger

% 3. Positionsvektor (inkluderar nu hela strängen från x=0 till x=1)
x = dx * (0:N)';

% Applicering av begynnelsedata (gaussisk puls)
g = @(x) exp(-200*(x-0.5).^2);
u(:,1) = g(x);
p(:,1) = zeros(N+1, 1); % Startar från vila

% FÖRE LOOPEN:
% Skapa energivektorn
E = zeros(1, M+1);

% Skapa viktmatris W för trapetsregeln (halv vikt på randnoderna)
weights = ones(N+1, 1);
weights(1) = 0.5;
weights(end) = 0.5;
W = diag(weights);

% Initial energi (före loopen) - Inget dx!
E(1) = 0.5 * (p(:,1)' * W * p(:,1)) - 0.5 * c^2 * (u(:,1)' * W * A_N * u(:,1));
% --- Förberedelse av visualisering ---
figure('Name', 'Vågpropagering (Neumann)');
% Till skillnad från 4c behöver vi inte manuellt lägga till nollor i kanterna
hPlot = plot(x, u(:,1), 'r', 'Linewidth', 1.5);
axis([0 1 -1 1]);
title('Tidsutveckling av vågekvationen (Fria ändar)');
xlabel('Position x'); 
ylabel('Amplitud u');
drawnow;

% --- Tidsstegning: Symplektisk Euler ---
for m = 1:M
    % Uppdatera hastigheten via den modifierade Neumann-matrisen
    p(:,m+1) = p(:,m) + dt * c^2 * A_N * u(:,m);
    
    % Uppdatera positionen explicit
    u(:,m+1) = u(:,m) + dt * p(:,m+1);
    
    % ...och inuti loopen:
    E(m+1) = 0.5 * (p(:,m+1)' * W * p(:,m+1)) - 0.5 * c^2 * (u(:,m+1)' * W * A_N * u(:,m+1));

    set(hPlot, 'YData', u(:,m+1));
    drawnow;
end

E
% --- Visualisering av energins utveckling över tid ---
figure('Name', 'Energins tidsutveckling (Neumann)');

% Skapa en tidsvektor från 0 till T, med M+1 punkter
t_vec = linspace(0, T, M+1); 

% Rita upp energin
plot(t_vec, E, 'k', 'LineWidth', 1.5);
grid on;

% Snygga till grafen med titlar och axel-etiketter
title('Total energi över tid (Symplektisk Euler, Neumann)');
xlabel('Tid t');
ylabel('Energi E');

% Frivilligt: Anpassa y-axeln för att zooma in på de små oscillationerna
% e_mean = mean(E);
% ylim([e_mean - 1e-4, e_mean + 1e-4]);
