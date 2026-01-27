function TrapzMain()
    % Huvudprogram för Uppgift 4
    
    % --- Uppgift 4a: Fel vs N ---
    x_vals = [0.11, 0.32, 1.14];
    N_vals = round(logspace(1, 3, 10)); % Logaritmiskt fördelade N
    
    figure(5); clf;
    for x = x_vals
        errors = [];
        for N = N_vals
            N
            integral_val = my_trapz(x, N);
            exact_val = erf(x); % Matlab's inbyggda 
            errors = [errors, abs(integral_val - exact_val)];
        end
        loglog(N_vals, errors, '-o', 'DisplayName', sprintf('x=%.2f', x)); hold on;
        
        % Teoretisk felgräns: C * x^3 / N^2
        % Max andradrivata av g(t) = 2/sqrt(pi) * e^(-t^2)
        % g'(t) = -2t * g(t)
        % g''(t) = (-2 + 4t^2) * g(t). Maximeras nära t=0 eller t=sqrt(1.5)?
        % För små x är max g''(0) = 2/sqrt(pi) * (-2) (absolutvärde ca 2.25)
        %C_theor = (x^3 / 12) * 2.25; % Approximativ konstant
        C = 1/(3*sqrt(pi)); 
        loglog(N_vals, C ./ N_vals.^2, '--', 'DisplayName', 'Teori 1/N^2');   
    end
    title('Fel i trapetsregeln vs N');
    xlabel('N'); ylabel('Absolut fel');
    legend; grid on;

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
    
    a = -x; b = x;
    h = (b-a)/N;
    t = linspace(a, b, N+1);
    
    g = @(v) (2/sqrt(pi)) * exp(-v.^2);
    y = g(t);
    
    % Trapetsformeln: h * (0.5*y0 + y1 + ... + 0.5*yN)
    T = sum(y) - 0.5*y(1) - 0.5*y(end);
    val = 0.5 * h * T; % Faktorn 0.5 pga symmetritipset i labben
end
