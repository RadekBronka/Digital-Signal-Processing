%% Parametry sygnału
fs = 1000;                % częstotliwość próbkowania [Hz] - zgodna z treścią zadania
fc = 200;                 % częstotliwość nośnej [Hz] - zgodna z treścią zadania
t = 0:1/fs:1-1/fs;        % czas trwania 1 sekunda - zgodny z treścią zadania

%% Wczytanie sygnału
load('lab08_am.mat');     % wczytanie zarejestrowanego sygnału AM z pliku (zgodnie z zadaniem)
% Wybór odpowiedniego sygnału na podstawie przedostatniej cyfry nr indeksu
% x = s0; x = s1; x = s2;
x = s0;                   % <-- ustaw tu swoją realizację (dopasuj do nr indeksu)

%% Projektowanie filtru Hilberta (FIR)
N = 201;                   % długość filtru (nieparzysta liczba próbek, typowe dla filtru FIR)
n = -(N-1)/2:(N-1)/2;     % oś próbkowania odpowiedzi impulsowej (indeksy próbek)

% Tworzenie idealnej odpowiedzi impulsowej filtru Hilberta
h_ideal = (1./(pi*n)).*(1 - cos(pi*n));  % h[n] = 1/(pi*n) * (1 - cos(pi*n))
h_ideal((N+1)/2) = 0;     % ustawienie wartości w n=0 (uniknięcie NaN) zgodnie z definicją

% Zastosowanie okna Hamminga (poprawa charakterystyki filtru - zgodnie z treścią zadania)
h = h_ideal .* hamming(N)';

% Filtracja sygnału x przez filtr Hilberta
delay = (N-1)/2;           % obliczenie opóźnienia filtru FIR
y = filter(h,1,x);         % filtracja sygnału
x_ht = [y(delay+1:end), zeros(1,delay)]; % kompensacja opóźnienia

%% Obliczenie obwiedni sygnału
env = sqrt(x.^2 + x_ht.^2);  % Obwiednia sygnału AM: sqrt(x^2 + HT(x)^2) ze wzoru

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

%% Analiza częstotliwościowa obwiedni i ekstrakcja parametrów m(t)
L = length(env);                         % Długość obwiedni
Env_fft = fft(env);                      % FFT obwiedni
mag = abs(Env_fft) / L;                  % Normalizacja widma amplitudowego
fv_full = (0:L-1)*(fs/L);                % Oś częstotliwości (częstotliwości)

% Pomijamy zerową częstotliwość (DC), bierzemy tylko dodatnie częstotliwości
fv_half = fv_full(2:floor(L/2));         % Oś częstotliwości bez DC
mag_half = mag(2:floor(L/2));            % Widmo obwiedni bez DC

% Wyświetlenie widma obwiedni
figure;
plot(fv_half, mag_half);                 % Wykres widma obwiedni
title('Widmo obwiedni sygnału AM');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda');
grid on;

% Pomijamy zerową częstotliwość (DC), bierzemy tylko dodatnie częstotliwości
[sorted_mag, idx_desc] = sort(mag_half, 'descend');  % Sortowanie malejąco
top_idx = idx_desc(1:3);               % Wybieramy trzy największe składowe

% Częstotliwości (indeksy muszą być przesunięte o 1, bo zaczynamy od 2 elementu)
f_sorted = fv_half(top_idx);            % Częstotliwości składowych
A_sorted = 2 * sorted_mag(1:3);        % Amplitudy (mnożymy przez 2, bo FFT jest symetryczne)

% Wyświetlenie wyników
fprintf('A1 = %.3f, f1 = %.1f Hz\n', A_sorted(1), f_sorted(1));
fprintf('A2 = %.3f, f2 = %.1f Hz\n', A_sorted(2), f_sorted(2));
fprintf('A3 = %.3f, f3 = %.1f Hz\n', A_sorted(3), f_sorted(3));

%% Rekonstrukcja sygnału x(t) na podstawie wyznaczonych parametrów m(t)
% Rekonstrukcja sygnału m(t) = 1 + A1*cos(2πf1t) + A2*cos(2πf2t) + A3*cos(2πf3t)
m_recon = 1 ...
    + A_sorted(1)*cos(2*pi*f_sorted(1)*t) ...
    + A_sorted(2)*cos(2*pi*f_sorted(2)*t) ...
    + A_sorted(3)*cos(2*pi*f_sorted(3)*t);

% Odtworzony sygnał AM: x_recon(t) = m(t) * cos(2π f_c t)
x_recon = m_recon .* cos(2*pi*fc*t);

%% Porównanie sygnału odtworzonego z oryginalnym
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
