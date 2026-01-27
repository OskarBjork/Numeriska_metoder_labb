function TrapzMain()
    % Huvudprogram för Uppgift 4
    
    % --- Uppgift 4a: Fel vs N ---
    % --- Uppgift 4a: Fel vs N (KORRIGERAD) ---
    x_vals = [0.11, 0.32, 1.14];
    N_vals = round(logspace(1.3, 2.7, 10)); % Logaritmiskt fördelade N
    
    figure(5); clf;
    colors = {'b', 'r', 'y'}; % Färger för de olika x-värdena [cite: 144]
    
    for i = 1:length(x_vals)
        x = x_vals(i);
        col = colors{i};
        errors = [];
        
        for N = N_vals
            integral_val = my_trapz(x, N);
            exact_val = erf(x); 
            errors = [errors, abs(integral_val - exact_val)];
        end
        
        % Plotta uppmätt fel (heldragen linje)
        loglog(N_vals, errors, ['-o' col], 'LineWidth', 1.5, ...
               'DisplayName', sprintf('x=%.2f', x)); hold on;
        
        % --- TEORETISK FELGRÄNS (JUSTERAD) ---
        % Maxvärde av g''(t) = 4/sqrt(pi) * e^(-t^2) * (2t^2-1) är vid t=0.
        % Maxvärdet är 4/sqrt(pi) ≈ 2.257
        % max_g_pp = 4/sqrt(pi); 
        
        % Formel: (b-a)^3 / (12*N^2) * max|g''|
        % Vid symmetri [-x, x] är b-a = 2x. Vi tar halva felet (0.5 * ...).
        % Teori = 0.5 * ( (2*x)^3 / (12*N_vals.^2) ) * max_g_pp
        
        % C_theor = 0.5 * ( (2*x)^3 / 12 ) * max_g_pp;
        
        % Plotta teoretisk gräns (streckad linje)
        % loglog(N_vals, C_theor ./ N_vals.^2, ['--' col], 'LineWidth', 1, ...
        %        'DisplayName', ['Teori 1/N^2 (x=' num2str(x) ')']);
    end
    
    title('Fel i trapetsregeln vs N');
    xlabel('N'); ylabel('Absolut fel');
    legend('Location', 'best'); grid on;

    % --- Uppgift 4b: Fel vs x (Fixed N) ---
    N_fixed = [50, 120, 400];
    x_range = linspace(0.1, 6, 100);
    
    figure(6); clf;
    for N = N_fixed
        errors_x = [];
        for x = x_range
            val = my_trapz(x, N);
            errors_x = [errors_x, abs(val - erf(x))];
        end
        semilogy(x_range, errors_x, '-', 'DisplayName', sprintf('N=%d', N)); hold on;
    end
    title('Fel i trapetsregeln vs x (Fix N)');
    xlabel('x'); ylabel('Absolut fel');
    legend; grid on;
    
    % Svar 4b: Felet sjunker drastiskt för stora x eftersom funktionen
    % och dess derivator går mot noll vid ränderna, vilket ger 
    % exponentiell konvergens för trapetsregeln.
end

function val = my_trapz(x, N)
    % Numerisk integration av Error function 2/sqrt(pi) * integral(e^-t^2)
    % Intervall [0, x], dock använder vi symmetrin [-x, x]/2 enligt tips [cite: 139]
    % för att utnyttja trapetsregelns egenskaper bättre.
    
    a = 0; b = x;
    h = (b-a)/N;
    t = linspace(a, b, N+1);
    
    g = @(v) (2/sqrt(pi)) * exp(-v.^2);
    y = g(t);
    
    % Trapetsformeln: h * (0.5*y0 + y1 + ... + 0.5*yN)
    T = sum(y) - 0.5*y(1) - 0.5*y(end);
    val = h * T; % Faktorn 0.5 pga symmetritipset i labben
end
