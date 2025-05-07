clc; 
clear;

%% Parametry
fs = 1000;        % częstotliwość próbkowania
fc = 200;         % częstotliwość nośna
t = 0:1/fs:1-1/fs; % wektor czasu
N = length(t);    

%% Wczytanie danych
load('lab08_am.mat'); % Zmienna np. x (macierz: próbki x różnych realizacji)

% Ustaw swój numer realizacji (np. 5, jeśli przedostatnia cyfra legitymacji to 5)
x_selected = s0;

%% Filtr Hilberta FIR (projekt filtru)
order = 100; % rząd filtru (parzysty)
h = firpm(order, [0.05 0.95], [1 1], 'hilbert'); % pasmo 5%-95% Nyquista

% Transformacja Hilberta
x_hilbert = conv(x_selected, h, 'same'); % HT(x)

%% Obwiednia sygnału (czyli m(t))
envelope = sqrt(x_selected.^2 + x_hilbert.^2);

%% Analiza częstotliwościowa obwiedni
m = envelope;
m_ac = m - mean(m); % usunięcie składowej stałej (1)

M_fft = abs(fft(m_ac));
f_axis = (0:N-1)*(fs/N);

% Wyznacz składowe częstotliwościowe powyżej np. 1 Hz
[~, locs] = findpeaks(M_fft(2:N/2), 'SortStr', 'descend', 'NPeaks', 3);
locs = locs + 1; % indeksy względem pełnego FFT
f1 = f_axis(locs(1));
f2 = f_axis(locs(2));
f3 = f_axis(locs(3));
A1 = 2*M_fft(locs(1))/N;
A2 = 2*M_fft(locs(2))/N;
A3 = 2*M_fft(locs(3))/N;

%% Posortuj częstotliwości rosnąco (opcjonalnie)
[f_sorted, idx] = sort([f1, f2, f3]);
A_sorted = [A1, A2, A3];
A_sorted = A_sorted(idx);

f1 = f_sorted(1); f2 = f_sorted(2); f3 = f_sorted(3);
A1 = A_sorted(1); A2 = A_sorted(2); A3 = A_sorted(3);

%% Wyświetl wyniki
fprintf('f1 = %.2f Hz, A1 = %.4f\n', f1, A1);
fprintf('f2 = %.2f Hz, A2 = %.4f\n', f2, A2);
fprintf('f3 = %.2f Hz, A3 = %.4f\n', f3, A3);

%% (Opcjonalnie) Rekonstrukcja sygnału x
m_recon = 1 + A1*cos(2*pi*f1*t) + A2*cos(2*pi*f2*t) + A3*cos(2*pi*f3*t);
x_recon = m_recon .* cos(2*pi*fc*t);

% Porównanie oryginału i rekonstrukcji
figure;
subplot(3,1,1);
plot(t, x_selected);
title('Oryginalny sygnał x(t)');
xlabel('Czas [s]'); ylabel('Amplituda');

subplot(3,1,2);
plot(t, envelope);
title('Obwiednia sygnału m(t)');
xlabel('Czas [s]'); ylabel('Amplituda');

subplot(3,1,3);
plot(t, x_recon);
title('Zrekonstruowany sygnał x(t)');
xlabel('Czas [s]'); ylabel('Amplituda');

% Porównanie błędu
mse = mean((x_selected - x_recon).^2);  % oblicz średni błąd rekonstrukcji
fprintf('Średni błąd rekonstrukcji (MSE): %.6f\n', mse);
