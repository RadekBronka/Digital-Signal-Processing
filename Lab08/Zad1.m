%% Parametry sygnału
fs = 1000;                % częstotliwość próbkowania [Hz] - zgodna z treścią zadania
fc = 200;                 % częstotliwość nośnej [Hz] - zgodna z treścią zadania
t = 0:1/fs:1-1/fs;        % czas trwania 1 sekunda - zgodny z treścią zadania

%% Wczytanie sygnału
load('lab08_am.mat');     % wczytanie zarejestrowanego sygnału AM z pliku (zgodnie z zadaniem)
% Wybór odpowiedniego sygnału na podstawie przedostatniej cyfry nr indeksu
% x = s0; x = s1; x = s2;
 x = s3;                   % <-- ustaw tu swoją realizację (dopasuj do nr indeksu)

%% Projektowanie filtru Hilberta (FIR)
N = 201;                   % długość filtru (nieparzysta liczba próbek, typowe dla filtru FIR)
n = -(N-1)/2:(N-1)/2;     % oś próbkowania odpowiedzi impulsowej

% Tworzenie idealnej odpowiedzi impulsowej filtru Hilberta
h_ideal = (1./(pi*n)).*(1 - cos(pi*n));  % h[n] = 1/(pi*n) * (1 - cos(pi*n))
h_ideal((N+1)/2) = 0;     % ustawienie wartości w n=0 (uniknięcie NaN) zgodnie z definicją

% Zastosowanie okna Hamminga (poprawa charakterystyki filtru - zgodnie z treścią zadania)
h = h_ideal .* hamming(N)';

% Filtracja sygnału x przez filtr Hilberta
delay = (N-1)/2;           % obliczenie opóźnienia filtru FIR
y = filter(h,1,x);         % filtracja sygnału
x_ht = [y(delay+1:end), zeros(1,delay)]; % kompensacja opóźnienia (zgodnie z podręcznikiem TZ, str. 353)

%% Obliczenie obwiedni sygnału
env = sqrt(x.^2 + x_ht.^2);  % Obwiednia sygnału AM: sqrt(x^2 + HT(x)^2) = |sygnał analityczny|

%% Wykresy sygnałów w dziedzinie czasu
% Sekcja wizualizacyjna: prezentacja sygnału x(t), jego transformaty Hilberta oraz obwiedni m(t)
figure; 
subplot(3,1,1);
plot(t, x);
title('Oryginalny sygnał x(t)');
xlabel('Czas [s]'); ylabel('Amplituda');

subplot(3,1,2);
plot(t, x_ht);
title('Transformata Hilberta HT(x)');
xlabel('Czas [s]'); ylabel('Amplituda');

subplot(3,1,3);
plot(t, env);
title('Obwiednia sygnału (m(t))');
xlabel('Czas [s]'); ylabel('Amplituda');

%% Wykres widma obwiedni
% Analiza częstotliwościowa obwiedni m(t)
L = length(env);
Env_fft = fft(env);                   % FFT obwiedni
mag = abs(Env_fft)/L;                % Normalizacja widma amplitudowego
fv_full = (0:L-1)*(fs/L);            % Oś częstotliwości
half = ceil(L/2);                    % Połowa widma (z uwagi na symetrię dla sygnałów rzeczywistych)
figure;
plot(fv_full(1:half), mag(1:half));  % Wykres widma obwiedni
title('Widmo obwiedni');
xlabel('Częstotliwość [Hz]'); ylabel('Amplituda');

%% FFT obwiedni i ekstrakcja parametrów m(t)
% Znalezienie dominujących składowych częstotliwościowych m(t)
L = length(env);
Env_fft = fft(env);
mag = abs(Env_fft)/L;
fv = (0:L-1)*(fs/L);
half = ceil(L/2);
mag_half = mag(2:half);      % pomijamy DC (pierwszy element)
fv_half = fv(2:half);        % odpowiadające częstotliwości

% Wyznaczenie trzech największych składowych
[sorted_mag, idx_desc] = sort(mag_half, 'descend');  % sortujemy malejąco
top_idx = idx_desc(1:3);     % wybieramy 3 największe

% Wyznaczenie amplitud (uwzględniając symetrię) i częstotliwości
A_vals = 2 * sorted_mag(1:3);  % przemnożenie przez 2 (bo FFT tylko połówka widma)
f_vals = fv_half(top_idx);

% Posortowanie częstotliwości rosnąco (f1 < f2 < f3)
[f_sorted, order] = sort(f_vals);
A_sorted = A_vals(order);

%% Wyświetlenie wyników parametrów
for k = 1:3
    fprintf('A%d = %.3f, f%d = %.1f Hz', k, A_sorted(k), k, f_sorted(k));
end


%% REKONSTRUKCJA SYGNAŁU x(t) NA PODSTAWIE WYZNACZONYCH PARAMETRÓW m(t) (opcjonalne +0.25 pkt)

% --- UWAGA ---
% Zakładamy, że znamy już wyznaczone parametry:
% A_sorted = [A1, A2, A3], f_sorted = [f1, f2, f3]
% oraz nośną fc i czas t

% Rekonstrukcja sygnału m(t) = 1 + A1*cos(2πf1t) + A2*cos(2πf2t) + A3*cos(2πf3t)
m_recon = 1 ...
    + A_sorted(1)*cos(2*pi*f_sorted(1)*t) ...
    + A_sorted(2)*cos(2*pi*f_sorted(2)*t) ...
    + A_sorted(3)*cos(2*pi*f_sorted(3)*t);

% Odtworzony sygnał AM: x_recon(t) = m(t) * cos(2π f_c t)
x_recon = m_recon .* cos(2*pi*fc*t);

%% PORÓWNANIE SYGNAŁU ODTWORZONEGO Z ORYGINALNYM

% Tworzymy wykres porównawczy
figure;
subplot(2,1,1);
plot(t, x);
title('Oryginalny sygnał x(t) z pliku');
xlabel('Czas [s]'); ylabel('Amplituda');

subplot(2,1,2);
plot(t, x_recon);
title('Odtworzony sygnał x_{recon}(t) na podstawie m(t)');
xlabel('Czas [s]'); ylabel('Amplituda');

% Obliczamy błąd (np. MSE – średni błąd kwadratowy), żeby ocenić jakość rekonstrukcji
mse = mean((x - x_recon).^2);
fprintf('\nŚredni błąd rekonstrukcji (MSE) = %.6f\n', mse);
