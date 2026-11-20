# MATLAB / GNU Octave

EN: The implementation is based on the author's v7 archive, **extended** with the 75% adaptive / 25% exploratory interval selection. The additional option `udzial_eksploracji` (default `0.25`) specifies the independent probability that the next point is the exact midpoint of the longest eligible interval. Set `ziarno_losowe` for repeatable pseudo-random choices. `info.historia_trybow` records the choice type for each new selected point. This changes the original v7 sampling logic; set `udzial_eksploracji=0` to disable it.

PL: Implementacja opiera się na dostarczonym archiwum v7 i dodaje **strategię mieszaną**: 75% prób adaptacyjnych i 25% prób eksploracyjnych. Parametr `udzial_eksploracji` (domyślnie `0.25`) jest prawdopodobieństwem wyboru dokładnego środka najdłuższego dopuszczalnego przedziału. `ziarno_losowe` pozwala odtworzyć sekwencję losową, a `info.historia_trybow` zapisuje sposób wyboru każdego punktu. Aby wyłączyć eksplorację, należy ustawić `udzial_eksploracji=0`.

Run / Uruchom: `test_probkowania` and `przyklad_probkowania`. GIF export needs a graphics backend with `getframe` and `imwrite` support. / Eksport GIF wymaga odpowiedniej obsługi graficznej.
