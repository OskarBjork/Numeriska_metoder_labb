% Skapande av referenslösning med ode45 (adaptiv Runge-Kutta)

% Fysikaliska och numeriska parametrar
a = 0.5; % Excentricitet som styr banans form
t_end = 100; % Total simuleringstid
% Begynnelsevillkor: q1, q2, p1, p2 (kroppen befinner sig vid perihelium)
y0 = [1-a; 0; 0; sqrt((1+a)/(1-a))];

f = @(y) [y(3); % p1 = q1_p
          y(4); % p2 = q2_p
          -y(1)/norm(y(1:2))^3; % p1_p
          -y(2)/norm(y(1:2))^3]; % p2_p

Jf = @(y) [0, 0, 1, 0; 
           0, 0, 0, 1; 
           (2*y(1)^2 - y(2)^2)/norm(y(1:2))^5, 3*y(1)*y(2)/norm(y(1:2))^5, 0, 0; 
           3*y(1)*y(2)/norm(y(1:2))^5, (2*y(2)^2 - y(1)^2)/norm(y(1:2))^5, 0, 0];
energy = @(y) 0.5*(y(3)^2 + y(4)^2) - 1/norm(y(1:2));
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
tic;
[~, y_ref] = ode45(@(t,y) f(y), [0 t_end], y0, options);
time_ode45 = toc;
y_true = y_ref(end, :)'; % Den "sanna" lösningsvektorn vid t=100

% Parametervektor för konvergensstudie
h_vals = [0.1, 0.05, 0.025, 0.0125];
err_se_y = zeros(size(h_vals));
err_se_energy = zeros(size(h_vals));
err_mp_y = zeros(size(h_vals));
err_mp_energy = zeros(size(h_vals));
time_se = zeros(size(h_vals));
time_mp = zeros(size(h_vals));

for i = 1:length(h_vals)
    h = h_vals(i);
    N = round(t_end/h);
    
    % Test av Symplektisk Euler
    tic;
    y = y0;
    for n = 1:N
        % p uppdateras explicit från q
        p_next = y(3:4) - h * (y(1:2) / norm(y(1:2))^3);
        % q uppdateras explicit från det nya p
        y(1:2) = y(1:2) + h * p_next;
        y(3:4) = p_next;
    end
    time_se(i) = toc;
    err_se_y(i) = norm(y - y_true); % Euklidiskt avstånd (L2-norm)
    err_se_energy(i) = norm(energy(y)- energy(y_true));
    
    % Test av Implicit Mittpunktsmetod
    tic;
    y = y0;
    for n = 1:N
        yg = y + h * f(y); % Initialgissning med framåt Euler
        % Newton iteration
        for iter = 1:7
            ym = 0.5 * (y + yg); % Mittpunkt
            J = eye(4) - 0.5 * h * Jf(ym); % Analytisk Jacobian
            F_val = yg - y - h * f(ym);
            delta = J \ F_val;
            yg = yg - delta;
            if norm(delta) < 1e-8, break; end
        end
        y = yg;
    end
    time_mp(i) = toc;
    err_mp_y(i) = norm(y - y_true);
    err_mp_energy(i) = norm(energy(y)-energy(y_true));
end

err_se_y;
err_mp_y;
err_se_energy;
err_mp_energy;

time_se
time_mp
time_ode45


% Beräkning av den empiriska konvergensordningen
ord_se = log2(err_se_y(1:end-1)./ err_se_y(2:end));
ord_mp = log2(err_mp_y(1:end-1)./ err_mp_y(2:end));