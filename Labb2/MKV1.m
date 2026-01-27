function MKV()
    % --- Uppgift 1a: Linjär modell ---
    % Ladda data (ersätt med din fil)
    load dollarkurs.mat; % Antas innehålla vektorn X (730x1) [cite: 14]
    
    X = USDSEK
    N = length(X);
    t = (1:N)'; % Tidsvektor
    
    % Modell: f(t) = c0 + c1*t
    % System matris A * c = X
    A = [ones(N,1), t]; 
    
    % Lösning med backslash
    c = A \ X; 
    c0 = c(1); c1 = c(2);
    
    % Beräkna modell och fel
    f_lin = A * c;
    resid_lin = X - f_lin;
    E_lin = (1/N) * sum(resid_lin.^2); % Medelkvadratfel 
    
    fprintf('Uppgift 1a - Linjär modell:\n');
    fprintf('c0 = %.4f, c1 = %.4e, Medelkvadratfel E = %.4e\n', c0, c1, E_lin);
    
    figure(1);
    subplot(2,1,1); plot(t, X, 'k.', t, f_lin, 'r-'); title('Data och Linjär anpassning');
    subplot(2,1,2); plot(t, resid_lin); title('Residualer (Fel)');
    
    % --- Uppgift 1b: Normalekvationer vs Backslash ---
    c_norm = (A'*A) \ (A'*X); % Via normalekvationer
    diff_norm = norm(c - c_norm);
    cond_ATA = cond(A'*A);
    
    fprintf('\nUppgift 1b:\n');
    fprintf('Skillnad i koefficienter (Backslash vs Normal): %.4e\n', diff_norm);
    fprintf('Konditionstal cond(A''A): %.4e\n', cond_ATA);
    
    % --- Uppgift 1c: Periodisk modell (L fixerad) ---
    % Uppskatta L från plotten i 1a (t.ex. räkna toppar). 
    % Låt oss gissa L baserat på en okulärbesiktning (exempelvärde, justera efter plot)
    L_guess = 460; 
    
    % Modell: f(t) = d0 + d1*t + d2*sin(2pi*t/L) + d3*cos(2pi*t/L) [cite: 32]
    A_per = [ones(N,1), t, sin(2*pi*t/L_guess), cos(2*pi*t/L_guess)];
    d = A_per \ X;
    
    f_per = A_per * d;
    resid_per = X - f_per;
    E_per = (1/N) * sum((X - f_per).^2);
    

    fprintf('\nUppgift 1c - Periodisk modell (L=%d):\n', L_guess);
    fprintf('E = %.4e\n', E_per);
    
    % --- Uppgift 1d: Gauss-Newton (L är variabel) ---
    % Parametervektor p = [d0; d1; d2; d3; L]
    p = [d; L_guess]; 
    tol = 1e-6; max_iter = 20;
    
    for k = 1:max_iter
        d0=p(1); d1=p(2); d2=p(3); d3=p(4); L=p(5);
        
        arg = 2*pi*t/L;
        % Residualer r = X - f(t)
        f_val = d0 + d1*t + d2*sin(arg) + d3*cos(arg);
        r = X - f_val; 
        
        % Jacobian J (N x 5) - Derivator map parametrarna
        % df/dd0 = 1, df/dd1 = t, etc.
        % df/dL = d2*cos(arg)*(-2pi*t/L^2) + d3*(-sin(arg))*(-2pi*t/L^2)
        
        d_inner = -2*pi*t ./ (L^2); % Inre derivata map L
        J_L = d2*cos(arg).*d_inner - d3*sin(arg).*d_inner;
        
        J = [ones(N,1), t, sin(arg), cos(arg), J_L];
        
        % Gauss-Newton steg: h = (J'J)^(-1) J'r
        h = J \ r; 
        p = p + h;
        
        if norm(h) < tol
            break;
        end
    end
    
    f_gn = d0 + d1*t + d2*sin(arg) + d3*cos(arg); % Uppdaterad modell
    E_gn = (1/N) * sum((X - f_gn).^2);
    resid_gn = X - f_gn;
    
    fprintf('\nUppgift 1d - Gauss-Newton:\n');
    fprintf('Optimerat L = %.4f\n', p(5));
    fprintf('Medelkvadratfel E = %.4e\n', E_gn);
    
    figure(2);
    plot(t, X, 'k.', t, f_lin, 'b--', t, f_per, 'g--', t, f_gn, 'r-', 'LineWidth', 1.5);
    %subplot(2,1,1); plot(t, X, 'k.', t, f_lin, 'r-'); title('Data och Linjär anpassning');
    legend('Data', 'Linjär', 'Periodisk (fast L)', 'Gauss-Newton (opt L)');
    title('Jämförelse av modeller');
    figure(3);
    plot(t, resid_per, t, resid_gn); title('test');
    E_lin, E_per, E_gn    
    
end

% linjär modell minskade 0,001715



