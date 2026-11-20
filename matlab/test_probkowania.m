% TEST_PROBKOWANIA Proste testy bez okna graficznego i bez toolboxow.
% EN: Regression tests of the MATLAB/Octave implementation.
clear;

% Test 1: funkcja liniowa, 3 probki wystarczaja bez maks_krok.
f = @(x,a) a.x0 + a.s*x;
a.x0 = -0.6;
a.s = 2.0;
opcje.maks_nowych = 20;
opcje.maks_punktow = 50;
opcje.tol_bezwzgl = 1e-9;
opcje.tol_wzgl = 0;
opcje.rysuj = false;
[x,y,info] = probkuj_adaptacyjnie(f, -2, 3, a, opcje);
assert(numel(x) == 3);
assert(info.liczba_nowych_wywolan == 3);
assert(max(abs(y-(a.x0+a.s*x))) < 1e-12);
assert(strcmp(info.powod_zatrzymania, 'tolerancja'));

% Test 2: wymuszona gestosc (mozliwe wykrywanie ukrytych pikow).
opcje.maks_nowych = 40;
opcje.maks_krok = 0.25;
[x1,y1,info1] = probkuj_adaptacyjnie(f, -2, 3, a, opcje);
assert(info1.maks_odstep <= 0.25 + 1e-12);

% Test 3: wznowienie nie zmienia wczesniejszych probek.
opcje.x0 = x1;
opcje.y0 = y1;
opcje.maks_krok = 0.10;
opcje.maks_nowych = 100;
opcje.maks_punktow = 200;
[x2,y2,info2] = probkuj_adaptacyjnie(f, -2, 3, a, opcje);
assert(all(ismember(x1, x2)));
for k = 1:numel(x1)
    j = find(x2 == x1(k), 1);
    assert(~isempty(j) && y2(j) == y1(k));
end
assert(info2.liczba_nowych_wywolan == numel(x2)-numel(x1));
assert(info2.maks_odstep <= 0.10 + 1e-12);

% Test 4: kwadratowa krzywa jest dogeszczana wedlug roznicy modeli.
a = [];
f2 = @(x,a) x.^2;
opcje2.maks_nowych = 12;
opcje2.maks_punktow = 15;
opcje2.tol_bezwzgl = 0;
opcje2.tol_wzgl = 0;
[x3,y3,info3] = probkuj_adaptacyjnie(f2, 0, 1, a, opcje2);
assert(numel(x3) > 3);
assert(numel(x3) == numel(y3));
assert(all(diff(x3) > 0));
assert(info3.liczba_nowych_wywolan == numel(x3));

fprintf('Wszystkie testy probkowania przeszly poprawnie.\n');

% Test 5: wszystkie metody sa dostepne i interpoluja wezly.
xw = [0; 0.2; 0.6; 1.1; 2.0];
yw = [0; 0.5; -0.1; 0.8; 0.2];
metody = {'linear','quadratic','pchip','spline'};
for j = 1:numel(metody)
    wyniki = interpoluj_probki(xw, yw, xw, metody{j});
    assert(max(abs(wyniki(:)-yw)) < 1e-10);
end

% Test 6: lokalna interpolacja kwadratowa dokladnie odtwarza parabole.
q = linspace(0,2,81);
yq = interpoluj_probki(xw, xw.^2-3*xw+1, q, 'quadratic');
assert(max(abs(yq(:)-(q(:).^2-3*q(:)+1))) < 1e-10);

% Test 7: wybrana para jest zapisywana; wznowienie bez ponownych wywolan.
opcje3.maks_nowych = 12;
opcje3.maks_punktow = 30;
opcje3.rysuj = false;
opcje3.interpolacja1 = 'linear';
opcje3.interpolacja2 = 'quadratic';
[x4,y4,info4] = probkuj_adaptacyjnie(@(x,a) sin(3*x),0,2,[],opcje3);
assert(strcmp(info4.interpolacja2,'quadratic'));
opcje3.x0 = x4; opcje3.y0 = y4;
opcje3.interpolacja1 = 'quadratic';
opcje3.interpolacja2 = 'spline';
[x5,y5,info5] = probkuj_adaptacyjnie(@(x,a) sin(3*x),0,2,[],opcje3);
assert(strcmp(info5.interpolacja1,'quadratic'));
assert(strcmp(info5.interpolacja2,'spline'));
assert(all(ismember(x4,x5)));
assert(info5.liczba_nowych_wywolan == numel(x5)-numel(x4));

fprintf('Testy dodatkowych interpolacji przeszly poprawnie.\n');

% Test 8: probka srodkowa jest lekko losowo przesuwana.
% Punkt startowy pozostaje wewnatrz ograniczenia +/-4% zakresu.
opcje4.maks_nowych = 3;
opcje4.maks_punktow = 3;
opcje4.losowosc = 0.04;
[xx,yy,ii] = probkuj_adaptacyjnie(@(x,a) x.^2, 0, 10, [], opcje4);
assert(numel(xx) == 3);
assert(xx(2) >= 4.6 && xx(2) <= 5.4);
assert(ii.losowosc == 0.04);

% Test 9: losowosc = 0 odtwarza dokladnie srodek.
opcje4.losowosc = 0;
[xx0,~,~] = probkuj_adaptacyjnie(@(x,a) x.^2, 0, 10, [], opcje4);
assert(xx0(2) == 5);

% Test 10: wszystkie stare punkty zachowane po wznowieniu z losowoscia.
opcje4.x0 = xx; opcje4.y0 = yy;
opcje4.losowosc = 0.04;
opcje4.maks_nowych = 12; opcje4.maks_punktow = 15;
opcje4.maks_krok = 1;
[xx2,yy2,ii2] = probkuj_adaptacyjnie(@(x,a) x.^2, 0, 10, [], opcje4);
assert(all(ismember(xx, xx2)));
assert(all(diff(xx2) > 0));
assert(ii2.liczba_nowych_wywolan == numel(xx2)-numel(xx));
for k = 1:numel(xx)
    assert(yy2(xx2 == xx(k)) == yy(k));
end
fprintf('Testy losowego przesuniecia przeszly poprawnie.\n');

% Test 11: domyslnie 20 testow lokalnych w kazdym przedziale.
opcje5 = struct();
opcje5.maks_nowych = 3;
opcje5.maks_punktow = 3;
opcje5.losowosc = 0;
opcje5.rysuj = false;
[~,~,info6] = probkuj_adaptacyjnie(@(x,a) sin(2*x),0,1,[],opcje5);
assert(strcmp(info6.tryb_siatki, 'lokalna'));
assert(info6.liczba_testow == 20);

% Test 12: oba tryby dzialaja na mocno niejednorodnym podziale,
% obejmujacym bardzo krotki przedzial.
opcje6 = struct();
opcje6.x0 = [0; 0.1; 0.100001; 0.7; 1];
opcje6.y0 = sin(3*opcje6.x0);
opcje6.maks_nowych = 1;
opcje6.maks_punktow = 10;
opcje6.losowosc = 0;
opcje6.tol_bezwzgl = 0;
opcje6.tol_wzgl = 0;
opcje6.liczba_testow = 20;
for j = 1:2
    if j == 1
        opcje6.tryb_siatki = 'lokalna';
    else
        opcje6.tryb_siatki = 'globalna';
        opcje6.liczba_testow_globalnych = 30;
    end
    [xt,yt,it] = probkuj_adaptacyjnie(@(x,a) sin(3*x),0,1,[],opcje6);
    assert(numel(xt) == numel(opcje6.x0)+1);
    assert(all(diff(xt) > 0));
    assert(all(isfinite(yt)));
    assert(strcmp(it.tryb_siatki, opcje6.tryb_siatki));
end
fprintf('Testy lokalnej i globalnej siatki przeszly poprawnie.\n');

% Test 13: strategia wierzcholka paraboli dla gladkiego maksimum.
op = struct(); op.maks_nowych = 7; op.maks_punktow = 10;
op.losowosc = 0; op.tol_bezwzgl = 0; op.tol_wzgl = 0;
op.wybor_punktu = 'wierzcholek';
op.udzial_eksploracji = 0; % EN: Disable exploration to test vertex logic.
[xv,yv,iv] = probkuj_adaptacyjnie(@(x,a) -(x-0.37).^2,0,1,[],op);
assert(any(abs(xv-0.37)<1e-10));
assert(iv.liczba_wybranych_wierzcholkow >= 1);
assert(abs(iv.iloraz_odstepow - mean(diff(xv))/min(diff(xv))) < 1e-10);

% Test 14: krzywe interpolowane dokladnie przechodza przez wezly.
for j=1:numel(metody)
    w=interpoluj_probki(xv,yv,unique([linspace(0,1,40)'; xv]),metody{j});
    assert(all(isfinite(w)));
    assert(max(abs(interpoluj_probki(xv,yv,xv,metody{j})-yv))<1e-10);
end
fprintf('Testy wyboru wierzcholka i zgodnosci wezlow przeszly.\n');

% Test 15: regresja bledu v6. Dlugi, prawie liniowy fragment nie moze
% przeslonic istotnej roznicy interpolantow w krotszym przedziale.
% Poczatkowe wezly sa zachowane. Przy maks_nowych=1 porownujemy wybor.
op = struct();
op.x0 = [0; 0.3; 0.5; 0.7; 1.0; 2.0];
op.y0 = [0; 0; 1; 0; 0; 0];
op.maks_nowych = 1;
op.maks_punktow = 7;
op.maks_krok = 0.4;
op.tol_bezwzgl = 1e-8;
op.tol_wzgl = 0;
op.losowosc = 0;
op.interpolacja1 = 'linear';
op.interpolacja2 = 'quadratic';
op.priorytet_wyboru = 'rozbieznosc';
op.udzial_eksploracji = 0; % EN: Priority regression is deterministic.
[~,~,ir] = probkuj_adaptacyjnie(@(t,a) exp(-((t-0.5)/0.12).^2),0,2,[],op);
assert(ir.historia_nowych_x(1) < 1.0);
op.priorytet_wyboru = 'krok';
[~,~,ik] = probkuj_adaptacyjnie(@(t,a) exp(-((t-0.5)/0.12).^2),0,2,[],op);
assert(ik.historia_nowych_x(1) > 1.0);
fprintf('Regresja priorytetu rozbieznosci przeszla poprawnie.\n');

% Test 16: przy 100%% eksploracji nowy punkt to srodek najdluzszego przedzialu.
% EN: A forced exploratory draw bisects the longest eligible interval exactly.
op = struct();
op.x0=[0;0.5;1;4]; op.y0=[0;1;0;0];
op.maks_nowych=1; op.maks_punktow=5;
op.udzial_eksploracji=1; op.losowosc=0.20;
op.wybor_punktu='wierzcholek'; op.ziarno_losowe=2026;
op.tol_bezwzgl=0; op.tol_wzgl=0;
[~,~,ie] = probkuj_adaptacyjnie(@(x,a) x.^2,0,4,[],op);
assert(abs(ie.historia_nowych_x(1)-2.5)<1e-12);
assert(ie.liczba_eksploracji==1);
assert(strcmp(ie.historia_trybow{1},'eksploracja'));
fprintf('Test eksploracji najdluzszego przedzialu przeszedl.\n');
