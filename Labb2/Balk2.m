function Balk()
    % Parametrar
    A = 9; alpha = 0.5; omega = 4;
    
    % Funktioner för y(t) och y'(t)
    y = @(t) A * exp(-alpha*t) .* cos(omega*t);
    y_prim = @(t) A * (-alpha * exp(-alpha*t) .* cos(omega*t) - omega * exp(-alpha*t) .* sin(omega*t));
    
    tol = 1e-8;
    
    % --- Uppgift 2a: Newtons metod (H=0.5) ---
    H1 = 0.5;
    f1 = @(t) y(t) - H1;
    df1 = y_prim;
    
    fprintf('--- Uppgift 2a (Newton H=0.5) ---\n');
    t = 1.5; % Startgissning (behöver vara rimlig, titta på graf om osäker)
    for k = 1:100
        val = f1(t);
        fprintf('Iter %d: t = %.10f, fel = %.2e\n', k, t, abs(val));
        if abs(val) < tol, break; end
        t = t - val / df1(t);
    end
    t_H1 = t;
    
    % --- Uppgift 2b: Sekantmetoden (H=0.5) ---
    fprintf('\n--- Uppgift 2b (Sekant H=0.5) ---\n');
    t0 = 1.4; t1 = 1.6; % Två startgissningar
    for k = 1:100
        f0_val = f1(t0);
        f1_val = f1(t1);
        if abs(f1_val) < tol, break; end
        
        t_new = t1 - f1_val * (t1 - t0) / (f1_val - f0_val);
        t0 = t1;
        t1 = t_new;
        fprintf('Iter %d: t = %.10f, fel = %.2e\n', k, t1, abs(f1_val));
    end
    
    % --- Uppgift 2c: Newton (H = 4.13...) ---
    H2 = 4.1355432362;
    f2 = @(t) y(t) - H2;
    fprintf('\n--- Uppgift 2c (Newton H=4.13...) ---\n');
    t = 0.5; % Startgissning
    for k = 1:100
        val = f2(t);
        deriv = df1(t); 
        fprintf('Iter %d: t = %.10f, derivata = %.4e\n', k, t, deriv);
        if abs(val) < tol, break; end
        t = t - val / deriv;
    end
    
    % 2d Diskussion:
    % Vid H2 är derivatan mycket nära noll (vi är på en topp).
    % Det gör att t ändras kraftigt även för små fel i H.
end