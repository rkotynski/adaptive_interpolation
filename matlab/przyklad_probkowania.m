% PRZYKLAD_PROBKOWANIA / EXAMPLE OF ADAPTIVE SAMPLING
% PL: Funkcja z ekstremami o roznych szerokosciach; zapis GIF.
% EN: Wide peak, medium dip and a much narrower resonance; save a GIF.
% Uruchom / Run: przyklad_probkowania
clear; close all;

% EN: Reuse the same benchmark function as in examples/demo_python.py.
% PL: Ta sama funkcja co w przykladzie Python.
f = @(x,a) 0.04 + 0.035*x + ...
    1.00*exp(-((x-2.45)/1.10).^2) - ...
    0.72*exp(-((x-5.30)/0.36).^2) + ...
    1.18*exp(-((x-7.43)/0.105).^2);
% EN: Uniform start, without special points close to the narrow peak.
% PL: Rownomierny start, bez specjalnych probek przy waskim maksimum.
x0 = linspace(0,10,6)';
y0 = arrayfun(@(z) f(z,[]),x0);
opcje = struct();
opcje.x0=x0; opcje.y0=y0;
opcje.maks_nowych=90; opcje.maks_punktow=180;
opcje.tol_bezwzgl=2e-4; opcje.tol_wzgl=3e-4;
opcje.maks_krok=0.22;
opcje.interpolacja1='linear'; opcje.interpolacja2='quadratic';
opcje.tryb_siatki='lokalna'; opcje.liczba_testow=20;
opcje.priorytet_wyboru='rozbieznosc';
opcje.udzial_eksploracji=0.25; % EN: 25 percent longest-interval bisection.
opcje.ziarno_losowe=2026; % EN: Fixed seed for reproducibility.
opcje.wybor_punktu='roznica';
opcje.losowosc=0; opcje.rysuj=false;
[x,y,info] = probkuj_adaptacyjnie(f,0,10,[],opcje);

% EN: Reconstruct frames WITHOUT repeated calls to the expensive function.
% PL: Odtwarzamy klatki z zapamietanych wartosci, bez ponownego wywolywania f.
outdir=fullfile(fileparts(mfilename('fullpath')),'..','assets');
if ~exist(outdir,'dir'), mkdir(outdir); end
filename=fullfile(outdir,'adaptive_sampling_matlab.gif');
fig=figure('Name','Adaptive sampling / Probkowanie adaptacyjne','Color','w');
set(fig,'Position',[100 100 980 650]);
for k=0:info.liczba_nowych_wywolan
    if k==0
        xi=x0; yi=y0; added=NaN;
    else
        % EN: Supplied initial points are kept unchanged.
        % PL: Zachowujemy wszystkie dostarczone probki poczatkowe.
        xi=[x0;info.historia_nowych_x(1:k)];
        [xi,~]=sort(xi);
        % To reconstruct y, read values from the final output, not f(x).
        % Wartosci y bierzemy z koncowego wyniku, a nie z nowych wywolan f.
        [~,ids]=ismember(xi,x);
        yi=y(ids);
        added=info.historia_nowych_x(k);
    end
    q=unique([linspace(0,10,2000)';xi(:)]);
    yl=interpoluj_probki(xi,yi,q,'linear');
    yq=interpoluj_probki(xi,yi,q,'quadratic');
    subplot(2,1,1);
    plot(q,yl,'b-',q,yq,'r-','LineWidth',1.5); hold on;
    plot(xi,yi,'ko','MarkerFaceColor','k','MarkerSize',3);
    if isfinite(added)
        [~,j]=min(abs(xi-added));
        plot(xi(j),yi(j),'go','MarkerSize',10,'LineWidth',2);
    end
    hold off; grid on;
    xlim([0 10]); ylim([-0.95 1.45]); ylabel('f(x)');
    if k==0, tryb='start'; else, tryb=info.historia_trybow{k}; end
    title(sprintf('75%% adaptive + 25%% exploration: %d samples (+%d), %s',numel(xi),k,tryb));
    legend('Linear / liniowa','Quadratic / kwadratowa','Samples / probki','Location','northeast');
    subplot(2,1,2);
    plot(q,abs(yq-yl),'m-','LineWidth',1.2);
    xlim([0 10]); grid on; xlabel('x'); ylabel('|Q - L|');
    drawnow;
    try
        frame=getframe(fig);
        [im,map]=rgb2ind(frame2im(frame),256);
        if k==0
            imwrite(im,map,filename,'gif','LoopCount',Inf,'DelayTime',0.17);
        else
            imwrite(im,map,filename,'gif','WriteMode','append','DelayTime',0.17);
        end
    catch ME
        warning('GIF export unsupported by graphics backend: %s',ME.message);
        break;
    end
end
fprintf('Added %d points, exploratory choices %d; total %d. Stop: %s\n', ...
    info.liczba_nowych_wywolan,info.liczba_eksploracji,numel(x),info.powod_zatrzymania);
