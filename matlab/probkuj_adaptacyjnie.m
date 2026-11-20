function [x, y, info] = probkuj_adaptacyjnie(f, xmin, xmax, a, opcje)
% PROBKUJ_ADAPTACYJNIE Oszczedne probkowanie rzeczywistej funkcji f(x,a).
% EN: Adaptive sampling of a real-valued, computationally expensive f(x,a).
% EN: Original v7 interface and sampling logic have been preserved.
% EN: Disagreement between interpolants is only a heuristic, NOT an error bound.
%
% [x,y,info] = probkuj_adaptacyjnie(f,xmin,xmax,a,opcje)
%
% f       - uchwyt do funkcji przyjmujacej SKALARNE x oraz ustalone a;
%           wynik f(x,a) musi byc pojedyncza, skonczona liczba rzeczywista.
% xmin    - lewy koniec przedzialu.
% xmax    - prawy koniec przedzialu.
% a       - parametry przekazywane bez zmian do f(x,a) (moga byc []).
% opcje   - opcjonalna struktura z polami:
%   x0, y0          - wczesniej obliczone probki; kazda bedzie uzyta ponownie.
%   maks_nowych     - limit NOWYCH wywolan f w tym uruchomieniu (domyslnie 100).
%   maks_punktow    - limit lacznej liczby probek (domyslnie 200).
%   tol_bezwzgl     - bezwzgledna tolerancja roznicy interpolantow (1e-3).
%   tol_wzgl        - tolerancja wzgledna do rozpietosci probek (1e-2).
%   maks_krok       - maksymalny dopuszczalny odstep x (Inf = bez wymuszenia).
%   interpolacja1  - pierwsza metoda (domyslnie 'linear').
%   interpolacja2  - druga metoda (domyslnie 'pchip').
%                   Dostepne: 'linear', 'quadratic', 'pchip', 'spline'.
%   losowosc        - maks. losowe przesuniecie x jako czesc dlugosci
%                      sasiedniego przedzialu (domyslnie 0.04 = 4%%).
%                      Zakres 0..0.25; 0 wylacza losowosc.
%                      Dotyczy tez srodkowej probki startowej.
%   tryb_siatki    - 'lokalna' (domyslnie) lub 'globalna'.
%   liczba_testow   - liczba kandydatow WEWNATRZ kazdego przedzialu
%                    w trybie 'lokalna' (domyslnie 20).
%   liczba_testow_globalnych - liczba punktow wewnetrznych na calym
%                    [xmin,xmax] w trybie 'globalna' (domyslnie 400).
%   udzial_eksploracji - prawdopodobienstwo (domyslnie 0.25) wyboru
%                     srodka najdluzszego przedzialu.
%                     EN: Probability of longest-interval midpoint selection.
%   ziarno_losowe  - opcjonalne ziarno RNG (puste: pozostaw stan RNG).
%                     EN: Optional reproducible random seed.
%   priorytet_wyboru - 'rozbieznosc' (domyslnie) lub 'krok'.
%                     'rozbieznosc': najpierw adaptacja do krzywizny,
%                       a maks_krok dopiero po spelnieniu tolerancji.
%                     'krok': stara strategia, najpierw najdluzszy
%                       przedzial przekraczajacy maks_krok.
%   wybor_punktu    - 'roznica' (domyslnie) lub 'wierzcholek'.
%                     Wierzcholek lokalnej paraboli w wybranym przedziale,
%                     gdy istnieje wewnatrz niego; inaczej maksimum roznicy.
%   rysuj           - czy rysowac postep (false).
%   odswiezaj_co    - odswiezaj wykres co tyle nowych punktow (1).
%   opoznienie      - opoznienie animacji w sekundach (0).
%   numer_rysunku   - numer okna figure; puste oznacza nowe okno ([]).
%
% ALGORYTM:
% 1) Zaczyna od koncow i srodka przedzialu (jesli nie ma danych wstepnych).
% 2) Porownuje interpolanty w 20 punktach WEWNATRZ kazdego
%    przedzialu; opcjonalnie korzysta ze starej siatki globalnej.
% 3) Wybiera miejsce maksymalnej roznicy miedzy interpolantami.
% 4) Najpierw lokalizuje rozbieznosci modeli; gdy mieszcza sie w
%    tolerancji, opcjonalnie dogeszcza zbyt dlugie przedzialy.
%    Tryb zgodnosci 'krok' odtwarza stare pierwszenstwo maks_krok.
% 5) Dodaje niewielkie losowe przesuniecie polozenia punktu wzgledem
%    odstepu miedzy sasiadami (opcje.losowosc).
% 6) Przerywa po uzyskaniu tolerancji i wymaganej gestosci lub limitu.
%
% UWAGA: roznica interpolantow jest HEURYSTYKA bledu, nie scislym
% oszacowaniem bledu wzgledem nieznanej funkcji. Bez ograniczenia maks_krok
% mozna przeoczyc waskie rezonanse lub oscylacje miedzy probkami.
%
% Funkcja wykorzystuje tylko standardowe funkcje MATLAB i GNU Octave.

    if nargin < 4
        error('Nalezy podac f, xmin, xmax oraz parametry a (ewentualnie []).');
    end
    if nargin < 5 || isempty(opcje)
        opcje = struct();
    end
    if ~isa(f, 'function_handle')
        error('Argument f musi byc uchwytem do funkcji, np. @(x,a) ...');
    end
    if ~isscalar(xmin) || ~isscalar(xmax) || ~isreal(xmin) || ...
            ~isreal(xmax) || ~isfinite(xmin) || ~isfinite(xmax) || xmin >= xmax
        error('Wymagane sa skonczone konce przedzialu: xmin < xmax.');
    end
    if ~isstruct(opcje) || numel(opcje) ~= 1
        error('Opcje musza byc pojedyncza struktura.');
    end

    o.maks_nowych = pobierz(opcje, 'maks_nowych', 100);
    o.maks_punktow = pobierz(opcje, 'maks_punktow', 200);
    o.tol_bezwzgl = pobierz(opcje, 'tol_bezwzgl', 1e-3);
    o.tol_wzgl = pobierz(opcje, 'tol_wzgl', 1e-2);
    o.maks_krok = pobierz(opcje, 'maks_krok', Inf);
    o.liczba_testow = pobierz(opcje, 'liczba_testow', 20);
    o.tryb_siatki = pobierz(opcje, 'tryb_siatki', 'lokalna');
    o.liczba_testow_globalnych = pobierz(opcje, 'liczba_testow_globalnych', 400);
    o.losowosc = pobierz(opcje, 'losowosc', 0.04);
    o.interpolacja1 = pobierz(opcje, 'interpolacja1', 'linear');
    o.interpolacja2 = pobierz(opcje, 'interpolacja2', 'pchip');
    o.wybor_punktu = pobierz(opcje, 'wybor_punktu', 'roznica');
    o.priorytet_wyboru = pobierz(opcje, 'priorytet_wyboru', 'rozbieznosc');
    o.udzial_eksploracji = pobierz(opcje, 'udzial_eksploracji', 0.25);
    o.ziarno_losowe = pobierz(opcje, 'ziarno_losowe', []);
    o.rysuj = pobierz(opcje, 'rysuj', false);
    o.odswiezaj_co = pobierz(opcje, 'odswiezaj_co', 1);
    o.opoznienie = pobierz(opcje, 'opoznienie', 0);
    o.numer_rysunku = pobierz(opcje, 'numer_rysunku', []);

    if ~ischar(o.wybor_punktu) || size(o.wybor_punktu,1) ~= 1
        error('wybor_punktu musi byc napisem: roznica lub wierzcholek.');
    end
    o.wybor_punktu = lower(strtrim(o.wybor_punktu));
    if ~any(strcmp(o.wybor_punktu, {'roznica','wierzcholek'}))
        error('wybor_punktu: dostepne wartosci to roznica, wierzcholek.');
    end
    if ~ischar(o.priorytet_wyboru) || size(o.priorytet_wyboru,1) ~= 1
        error('priorytet_wyboru musi byc napisem: rozbieznosc lub krok.');
    end
    o.priorytet_wyboru = lower(strtrim(o.priorytet_wyboru));
    if ~any(strcmp(o.priorytet_wyboru, {'rozbieznosc','krok'}))
        error('priorytet_wyboru: dostepne wartosci to rozbieznosc, krok.');
    end
    o.interpolacja1 = sprawdz_metode(o.interpolacja1, 'interpolacja1');
    o.interpolacja2 = sprawdz_metode(o.interpolacja2, 'interpolacja2');
    if strcmp(o.interpolacja1, o.interpolacja2)
        error('interpolacja1 i interpolacja2 musza byc rozne.');
    end

    if ~isnumeric(o.udzial_eksploracji) || ~isscalar(o.udzial_eksploracji) || ...
            ~isreal(o.udzial_eksploracji) || ~isfinite(o.udzial_eksploracji) || ...
            o.udzial_eksploracji < 0 || o.udzial_eksploracji > 1
        error('udzial_eksploracji musi nalezec do [0,1].');
    end
    if ~isempty(o.ziarno_losowe)
        sprawdz_nieujemna_calkowita(o.ziarno_losowe, 'ziarno_losowe');
        rng(o.ziarno_losowe); % EN: Reproducible random choices.
    end
    sprawdz_nieujemna_calkowita(o.maks_nowych, 'maks_nowych');
    sprawdz_dodatnia_calkowita(o.maks_punktow, 'maks_punktow');
    sprawdz_dodatnia_calkowita(o.liczba_testow, 'liczba_testow');
    sprawdz_dodatnia_calkowita(o.liczba_testow_globalnych, 'liczba_testow_globalnych');
    if ~ischar(o.tryb_siatki) || size(o.tryb_siatki,1) ~= 1
        error('tryb_siatki musi byc napisem: lokalna lub globalna.');
    end
    o.tryb_siatki = lower(strtrim(o.tryb_siatki));
    if ~any(strcmp(o.tryb_siatki, {'lokalna','globalna'}))
        error('tryb_siatki musi byc rowny lokalna lub globalna.');
    end
    sprawdz_dodatnia_calkowita(o.odswiezaj_co, 'odswiezaj_co');
    if ~isnumeric(o.losowosc) || ~isscalar(o.losowosc) || ...
            ~isreal(o.losowosc) || ~isfinite(o.losowosc) || ...
            o.losowosc < 0 || o.losowosc > 0.25
        error('losowosc musi nalezec do przedzialu [0, 0.25].');
    end
    if o.maks_punktow < 3
        error('maks_punktow musi wynosic co najmniej 3.');
    end
    if ~isscalar(o.tol_bezwzgl) || ~isfinite(o.tol_bezwzgl) || o.tol_bezwzgl < 0 || ...
       ~isscalar(o.tol_wzgl) || ~isfinite(o.tol_wzgl) || o.tol_wzgl < 0
        error('Tolerancje musza byc skonczone i nieujemne.');
    end
    if ~isscalar(o.maks_krok) || ~isreal(o.maks_krok) || ...
            isnan(o.maks_krok) || o.maks_krok <= 0
        error('maks_krok musi byc dodatni (dopuszczalne Inf).');
    end
    if ~isscalar(o.opoznienie) || ~isfinite(o.opoznienie) || o.opoznienie < 0
        error('opoznienie musi byc nieujemne.');
    end
    if ~isscalar(o.rysuj) || ~(isnumeric(o.rysuj) || islogical(o.rysuj))
        error('rysuj powinno byc wartoscia logiczna.');
    end
    if ~isempty(o.numer_rysunku) && (~isscalar(o.numer_rysunku) || ...
       ~isfinite(o.numer_rysunku) || o.numer_rysunku < 1 || ...
       o.numer_rysunku ~= floor(o.numer_rysunku))
        error('numer_rysunku musi byc dodatnia liczba calkowita.');
    end

    jest_x0 = isfield(opcje, 'x0');
    jest_y0 = isfield(opcje, 'y0');
    if jest_x0 ~= jest_y0
        error('Dane wstepne wymagaja jednoczesnego podania x0 i y0.');
    end
    if jest_x0
        x = opcje.x0(:);
        y = opcje.y0(:);
    else
        x = zeros(0,1);
        y = zeros(0,1);
    end
    if numel(x) ~= numel(y) || ~isnumeric(x) || ~isnumeric(y) || ...
       ~isreal(x) || ~isreal(y) || any(~isfinite(x)) || any(~isfinite(y))
        error('x0 i y0 musza byc rownej dlugosci, skonczone i rzeczywiste.');
    end
    % Konwersja chroni przed zaokraglaniem nowych x przy calkowitych x0.
    x = double(x);
    y = double(y);
    if any(x < xmin) || any(x > xmax)
        error('Wszystkie punkty x0 musza nalezec do [xmin,xmax].');
    end
    [x, indeksy] = sort(x);
    y = y(indeksy);
    if any(diff(x) <= 0)
        error('Punkty x0 musza miec rozne wspolrzedne.');
    end

    liczba_starych = numel(x);
    liczba_nowych = 0;
    historia_x = zeros(0,1);
    historia_roznic = zeros(0,1);
    historia_trybow = cell(0,1); % EN: Chosen interval selection methods.
    liczba_wierzcholkow = 0;
    poprzednio_wybrany_x = NaN;
    poprzednio_wybrany_y = NaN;

    % Zachowujemy dostarczone pomiary. Uzupelniamy konce przedzialu,
    % a jesli bylo mniej niz 3 punkty - rowniez jego srodek.
    if numel(x) < 3
        wymagane = [xmin; xmin + 0.5*(xmax-xmin); xmax];
    else
        wymagane = [xmin; xmax];
    end
    braki = wymagane(~ismember(wymagane, x));
    if numel(x) + numel(braki) > o.maks_punktow || numel(braki) > o.maks_nowych
        error(['Limity probkowania nie pozwalaja uzupelnic punktow ' ...
               'startowych (koncow przedzialu i ewentualnie srodka).']);
    end
    for k = 1:numel(braki)
        x_now = braki(k);
        % Zlam symetrie takze przy pierwszym probkowaniu: konce
        % sa stale, losujemy jedynie brakujacy punkt srodkowy.
        srodek = xmin + 0.5*(xmax-xmin);
        if x_now == srodek && o.losowosc > 0
            x_now = przesun_losowo(x_now, xmin, xmax, o.losowosc);
        end
        y_now = wywolaj_f(f, x_now, a);
        [x, y] = wstaw_punkt(x, y, x_now, y_now);
        liczba_nowych = liczba_nowych + 1;
        historia_x(end+1,1) = x_now;
        historia_roznic(end+1,1) = NaN;
        poprzednio_wybrany_x = x_now;
        poprzednio_wybrany_y = y_now;
    end
    if numel(x) < 3
        error('Do porownania interpolantow potrzeba co najmniej 3 probek.');
    end

    wykres = [];
    if logical(o.rysuj)
        wykres = przygotuj_wykres(o.numer_rysunku, o.interpolacja1, o.interpolacja2);
        aktualizuj_wykres(wykres, x, y, liczba_nowych, ...
                          poprzednio_wybrany_x, poprzednio_wybrany_y, ...
                          o.interpolacja1, o.interpolacja2);
        if o.opoznienie > 0, pause(o.opoznienie); end
    end

    powod = 'limit_punktow';
    while true
        [rozbieznosc, kandydaci] = porownaj_interpolacje( ...
            x, y, o.liczba_testow, o.interpolacja1, o.interpolacja2, ...
            o.tryb_siatki, o.liczba_testow_globalnych);
        dlugosci = diff(x);
        prog = o.tol_bezwzgl + o.tol_wzgl * (max(y)-min(y));
        maks_roznica = max(rozbieznosc);
        maks_dl = max(dlugosci);

        % Wymuszone dogeszczenie jest niezalezne od wskaznika krzywizny.
        dlugie = find(dlugosci > o.maks_krok);
        if isempty(dlugie) && maks_roznica <= prog
            powod = 'tolerancja';
            break;
        end
        if liczba_nowych >= o.maks_nowych || numel(x) >= o.maks_punktow
            powod = 'limit_punktow';
            break;
        end

        % Unikamy identycznych liczb zmiennoprzecinkowych x po zaokragleniu.
        krok_maszynowy = 16 * eps(max([1, abs(xmin), abs(xmax)]));
        dopuszczalne = find(dlugosci > 2*krok_maszynowy);
        if isempty(dopuszczalne)
            powod = 'rozdzielczosc_numeryczna';
            break;
        end

        % EN: 25% exploratory steps bisect the longest eligible interval.
        % PL: W 25% iteracji wybieramy dokladny srodek najdluzszego
        % dopuszczalnego przedzialu, niezaleznie od roznicy interpolantow.
        % Pozostale 75% uzywa dotychczasowego wyboru adaptacyjnego.
        eksploracja = rand() < o.udzial_eksploracji;
        czy_wymuszony_krok = false;
        if strcmp(o.priorytet_wyboru, 'krok') && ~isempty(dlugie)
            czy_wymuszony_krok = true;
        elseif maks_roznica <= prog && ~isempty(dlugie)
            czy_wymuszony_krok = true;
        end
        if eksploracja
            [~, j] = max(dlugosci(dopuszczalne));
            i = dopuszczalne(j);
            x_now = x(i) + 0.5*dlugosci(i);
            tryb_wyboru = 'eksploracja';
        elseif czy_wymuszony_krok
            dlugie_dopuszczalne = intersect(dlugie, dopuszczalne);
            if isempty(dlugie_dopuszczalne)
                powod = 'rozdzielczosc_numeryczna';
                break;
            end
            [~, j] = max(dlugosci(dlugie_dopuszczalne));
            i = dlugie_dopuszczalne(j);
            x_now = x(i) + 0.5*dlugosci(i);
            tryb_wyboru = 'krok_awaryjny';
        else
            [~, j] = max(rozbieznosc(dopuszczalne));
            i = dopuszczalne(j);
            skala_maszynowa = 64 * eps(max([1; abs(y)]));
            if rozbieznosc(i) <= skala_maszynowa
                x_now = x(i) + 0.5*dlugosci(i);
            else
                x_now = kandydaci(i);
            end
            if strcmp(o.wybor_punktu, 'wierzcholek')
                [x_wierz, znaleziono] = wierzcholek_lokalnej_paraboli(x, y, i);
                if znaleziono
                    x_now = x_wierz;
                    liczba_wierzcholkow = liczba_wierzcholkow + 1;
                end
            end
            if o.losowosc > 0
                x_now = przesun_losowo(x_now, x(i), x(i+1), o.losowosc);
            end
            tryb_wyboru = 'rozbieznosc';
        end
        % EN: Midpoints are not jittered and never shifted to a vertex.
        % PL: Punktow eksploracyjnych nie przesuwamy losowo.
        if ~(x_now > x(i) && x_now < x(i+1))
            powod = 'rozdzielczosc_numeryczna';
            break;
        end

        % EN: This is the only expensive function evaluation in this iteration.
        % Tylko tutaj wykonujemy kosztowne wywolanie f(x,a).
        y_now = wywolaj_f(f, x_now, a);
        [x, y] = wstaw_punkt(x, y, x_now, y_now);
        liczba_nowych = liczba_nowych + 1;
        historia_x(end+1,1) = x_now;
        historia_roznic(end+1,1) = maks_roznica;
        historia_trybow{end+1,1} = tryb_wyboru;
        poprzednio_wybrany_x = x_now;
        poprzednio_wybrany_y = y_now;

        if logical(o.rysuj) && mod(liczba_nowych, o.odswiezaj_co) == 0
            aktualizuj_wykres(wykres, x, y, liczba_nowych, x_now, y_now, ...
                               o.interpolacja1, o.interpolacja2);
            if o.opoznienie > 0, pause(o.opoznienie); end
        end
    end

    % Stan po ostatnim wywolaniu, takze w przypadku limitu punktow.
    [rozbieznosc, ~] = porownaj_interpolacje(x, y, o.liczba_testow, ...
        o.interpolacja1, o.interpolacja2, o.tryb_siatki, ...
        o.liczba_testow_globalnych);
    info.interpolacja1 = o.interpolacja1;
    info.tryb_siatki = o.tryb_siatki;
    info.liczba_testow = o.liczba_testow;
    info.liczba_testow_globalnych = o.liczba_testow_globalnych;
    info.losowosc = o.losowosc;
    info.wybor_punktu = o.wybor_punktu;
    info.priorytet_wyboru = o.priorytet_wyboru;
    info.udzial_eksploracji = o.udzial_eksploracji;
    info.historia_trybow = historia_trybow;
    info.liczba_eksploracji = sum(strcmp(historia_trybow,'eksploracja'));
    info.liczba_adaptacyjnych = sum(strcmp(historia_trybow,'rozbieznosc'));
    info.liczba_krokow_awaryjnych = sum(strcmp(historia_trybow,'krok_awaryjny'));
    info.liczba_wybranych_wierzcholkow = liczba_wierzcholkow;
    info.interpolacja2 = o.interpolacja2;
    info.liczba_nowych_wywolan = liczba_nowych;
    info.liczba_wczesniejszych_probek = liczba_starych;
    info.liczba_punktow = numel(x);
    info.historia_nowych_x = historia_x;
    info.historia_maks_roznicy = historia_roznic;
    info.maks_roznica_interpolantow = max(rozbieznosc);
    info.maks_odstep = max(diff(x));
    info.min_odstep = min(diff(x));
    info.sredni_odstep = (x(end)-x(1))/(numel(x)-1);
    info.iloraz_odstepow = info.sredni_odstep/info.min_odstep;
    info.prog = o.tol_bezwzgl + o.tol_wzgl * (max(y)-min(y));
    info.powod_zatrzymania = powod;
    if logical(o.rysuj)
        aktualizuj_wykres(wykres, x, y, liczba_nowych, ...
                          poprzednio_wybrany_x, poprzednio_wybrany_y, ...
                          o.interpolacja1, o.interpolacja2);
    end
end

function [xw, ok] = wierzcholek_lokalnej_paraboli(x, y, i)
    % Ta sama trojka wezlow, co w interpolacji 'quadratic'.
    % Normalizacja x ogranicza utrate dokladnosci przy bardzo malych h.
    n = numel(x);
    k0 = max(1, min(i-1, n-2));
    k = k0:(k0+2);
    roz = x(k(3))-x(k(1));
    xw = NaN;
    ok = false;
    if ~(roz > 0), return; end
    z = (x(k)-x(k(1)))/roz;
    % Ilorazy roznicowe paraboli y = A*z^2+B*z+C.
    d01 = (y(k(2))-y(k(1)))/(z(2)-z(1));
    d12 = (y(k(3))-y(k(2)))/(z(3)-z(2));
    A = (d12-d01)/(z(3)-z(1));
    B = d01-A*(z(1)+z(2));
    if ~isfinite(A) || ~isfinite(B) || ...
            abs(A) <= 128*eps(max([1;abs(y(k))]))
        return;
    end
    zw = -B/(2*A);
    xw = x(k(1))+zw*roz;
    % Dodatkowy margines zapobiega wyborowi niemal istniejacego wezla.
    h = x(i+1)-x(i);
    marg = max(32*eps(max([1,abs(x(i)),abs(x(i+1))])),1e-6*h);
    ok = isfinite(xw) && xw > x(i)+marg && xw < x(i+1)-marg;
end

function x_nowe = przesun_losowo(x, lewy, prawy, losowosc)
    % Symetryczny szum o amplitudzie losowosc * (prawy-lewy).
    % Unika punktow przy samych wezlach (zageszczanie bez postepu).
    h = prawy-lewy;
    x_nowe = x + (2*rand()-1)*losowosc*h;
    margines = max(16*eps(max([1,abs(lewy),abs(prawy)])), 0.01*h);
    if 2*margines < h
        x_nowe = min(max(x_nowe, lewy+margines), prawy-margines);
    else
        x_nowe = lewy + 0.5*h;
    end
end

function metoda = sprawdz_metode(metoda, nazwa)
    if ~ischar(metoda) || size(metoda,1) ~= 1
        error('Opcja %s musi byc napisem (char).', nazwa);
    end
    metoda = lower(strtrim(metoda));
    if ~any(strcmp(metoda, {'linear','quadratic','pchip','spline'}))
        error('Nieznana metoda dla %s: %s.', nazwa, metoda);
    end
end

function wartosc = pobierz(opcje, pole, domyslnie)
    if isfield(opcje, pole) && ~isempty(opcje.(pole))
        wartosc = opcje.(pole);
    else
        wartosc = domyslnie;
    end
end

function sprawdz_dodatnia_calkowita(wartosc, nazwa)
    if ~isscalar(wartosc) || ~isreal(wartosc) || ~isfinite(wartosc) || ...
            wartosc < 1 || wartosc ~= floor(wartosc)
        error('Opcja %s musi byc dodatnia liczba calkowita.', nazwa);
    end
end

function sprawdz_nieujemna_calkowita(wartosc, nazwa)
    if ~isscalar(wartosc) || ~isreal(wartosc) || ~isfinite(wartosc) || ...
            wartosc < 0 || wartosc ~= floor(wartosc)
        error('Opcja %s musi byc nieujemna liczba calkowita.', nazwa);
    end
end

function y_now = wywolaj_f(f, x_now, a)
    y_now = f(x_now, a);
    if ~isnumeric(y_now) || ~isscalar(y_now) || ~isreal(y_now) || ...
            ~isfinite(y_now)
        error('f(x,a) musi zwracac skonczona, rzeczywista liczbe skalarna.');
    end
    y_now = double(y_now);
end

function [x, y] = wstaw_punkt(x, y, x_now, y_now)
    x(end+1,1) = x_now;
    y(end+1,1) = y_now;
    [x, indeksy] = sort(x);
    y = y(indeksy);
end

function [maksima, kandydaci] = porownaj_interpolacje(x, y, n_testow, m1, m2, tryb, n_globalnych)
    % EN: Candidate comparisons reuse existing samples; no calls to f(x,a).
    % Wszystkie interpolacje sa tanie: nie wywolujemy tutaj f(x,a).
    h = diff(x);
    liczba_przedzialow = numel(h);

    if strcmp(tryb, 'lokalna')
        % Dokladnie n_testow punktow wewnatrz KAZDEGO przedzialu.
        % Nawet bardzo krotkie przedzialy maja wlasna siatke.
        t = (1:n_testow)/(n_testow+1);
        XQ = bsxfun(@plus, x(1:end-1), h*t);
        q = XQ(:);
        y1 = interpoluj_probki(x, y, q, m1);
        y2 = interpoluj_probki(x, y, q, m2);
        R = reshape(abs(y2-y1), size(XQ));
        [maksima, indeksy] = max(R, [], 2);
        wiersze = (1:liczba_przedzialow)';
        kandydaci = XQ(sub2ind(size(XQ), wiersze, indeksy));
    else
        % Opcjonalna siatka jednorodna na calym przedziale.
        % Srodki lokalnych przedzialow gwarantuja co najmniej
        % jednego kandydata, nawet gdy siatka globalna omija przedzial.
        q_globalne = x(1) + ((1:n_globalnych)'/(n_globalnych+1)) * ...
            (x(end)-x(1));
        srodki = x(1:end-1) + 0.5*h;
        q = [q_globalne; srodki];
        roznice = abs(interpoluj_probki(x,y,q,m2) - ...
                       interpoluj_probki(x,y,q,m1));
        maksima = roznice(n_globalnych+1:end);
        kandydaci = srodki;
        for i = 1:liczba_przedzialow
            indeksy = find(q_globalne > x(i) & q_globalne < x(i+1));
            if ~isempty(indeksy)
                [r, j] = max(roznice(indeksy));
                if r > maksima(i)
                    maksima(i) = r;
                    kandydaci(i) = q_globalne(indeksy(j));
                end
            end
        end
    end
end

function w = przygotuj_wykres(numer, m1, m2)
    if isempty(numer)
        figure('Name', 'Adaptacyjne probkowanie funkcji', 'NumberTitle', 'off');
    else
        figure(numer);
    end
    clf;
    w.gorny = subplot(2,1,1);
    w.pierwszy = plot(NaN, NaN, 'k--', 'LineWidth', 1.0); hold on;
    w.drugi = plot(NaN, NaN, 'b-', 'LineWidth', 1.7);
    w.probki = plot(NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    w.ostatni = plot(NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 7);
    hold off;
    grid on;
    xlabel('x'); ylabel('f(x)');
    legend(['Interpolacja: ' m1], ['Interpolacja: ' m2], ...
           'Obliczone probki', 'Ostatni nowy punkt', 'Location', 'best');
    w.dolny = subplot(2,1,2);
    w.roznica = plot(NaN, NaN, 'b-', 'LineWidth', 1.2);
    grid on;
    xlabel('x'); ylabel(['|' m2 ' - ' m1 '|']);
    title('Rozbieznosc interpolantow (nie jest granica bledu)');
end

function aktualizuj_wykres(w, x, y, ile_nowych, nowy_x, nowy_y, m1, m2)
    % Jednorodna siatka ORAZ wszystkie wezly: nie omijamy ich na rysunku.
    q = unique([linspace(x(1), x(end), max(2000, 5*numel(x)))'; x(:)]);
    y1 = interpoluj_probki(x, y, q, m1);
    y2 = interpoluj_probki(x, y, q, m2);
    set(w.pierwszy, 'XData', q, 'YData', y1);
    set(w.drugi, 'XData', q, 'YData', y2);
    set(w.probki, 'XData', x, 'YData', y);
    if isfinite(nowy_x)
        set(w.ostatni, 'XData', nowy_x, 'YData', nowy_y);
    else
        set(w.ostatni, 'XData', NaN, 'YData', NaN);
    end
    set(w.roznica, 'XData', q, 'YData', abs(y2-y1));
    h = diff(x);
    wsp = ((x(end)-x(1))/(numel(x)-1))/min(h);
    title(w.gorny, sprintf('%d probek, %d nowych wywolan, sredni/min = %.3g', ...
                           numel(x), ile_nowych, wsp));
    drawnow;
end
