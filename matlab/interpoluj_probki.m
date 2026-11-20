function yq = interpoluj_probki(x, y, xq, metoda)
% INTERPOLUJ_PROBKI Interpolacja liniowa, kwadratowa, PCHIP lub spline.
% EN: Linear, local quadratic, PCHIP, or cubic spline interpolation.
% EN: The quadratic stencil uses the preceding node for interior intervals
% EN: (three successive samples); this intentionally uses the same three-node rule.
% yq = interpoluj_probki(x, y, xq, metoda)
% x musi byc rosnacym wektorem co najmniej 3 roznych wezlow.
% xq musi nalezec do przedzialu [x(1), x(end)].
% Metoda 'quadratic' stosuje lokalny wielomian Lagrange'a stopnia 2,
% oparty o trzy kolejne wezly; interpolant jest ciagly, ale jego
% pierwsza pochodna nie musi byc ciagla w wezlach.
    x = x(:);
    y = y(:);
    if numel(x) ~= numel(y) || numel(x) < 3 || any(diff(x) <= 0)
        error('Potrzebne sa co najmniej 3 rozne, rosnace wezly x i zgodne y.');
    end
    if ~ischar(metoda) || size(metoda,1) ~= 1
        error('Metoda musi byc napisem: linear, quadratic, pchip lub spline.');
    end
    metoda = lower(strtrim(metoda));
    if ~any(strcmp(metoda, {'linear','quadratic','pchip','spline'}))
        error('Nieznana metoda interpolacji: %s.', metoda);
    end
    if ~strcmp(metoda, 'quadratic')
        yq = interp1(x, y, xq, metoda);
        return;
    end
    rozmiar = size(xq);
    q = xq(:);
    if any(q < x(1)) || any(q > x(end))
        error('Dla quadratic punkty xq musza nalezec do przedzialu wezlow.');
    end
    n = numel(x);
    % Przedzial [x(i),x(i+1)] otrzymuje trojke najblizszych wezlow.
    % Dla i=1 sa to wezly 1:3; dla ostatniego n-2:n.
    yq_vec = zeros(size(q));
    for i = 1:n-1
        if i == n-1
            maska = (q >= x(i)) & (q <= x(i+1));
        else
            maska = (q >= x(i)) & (q < x(i+1));
        end
        if ~any(maska)
            continue;
        end
        k0 = max(1, min(i-1, n-2));
        k = k0:(k0+2);
        t = q(maska);
        % Lagrange: trzy skladniki, bez polyfit (lepsza stabilnosc
        % dla wezlow blisko siebie; skalujemy lokalna wspolrzedna).
        s = (x(k)-x(k(1))) / (x(k(3))-x(k(1)));
        z = (t-x(k(1))) / (x(k(3))-x(k(1)));
        wart = zeros(size(z));
        for j = 1:3
            p = mod(j,3)+1;
            r = mod(j+1,3)+1;
            wart = wart + y(k(j)) * ((z-s(p)).*(z-s(r))) / ...
                ((s(j)-s(p))*(s(j)-s(r)));
        end
        yq_vec(maska) = wart;
    end
    yq = reshape(yq_vec, rozmiar);
end
