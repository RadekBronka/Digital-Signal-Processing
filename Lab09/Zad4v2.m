clear;
close all;

load('ECG100.mat');      % Załaduj plik
ekg = val(1, :);         % Wybierz pierwszy kanał sygnału
fs = 1000;               % Częstotliwość próbkowania
t = (0:length(ekg)-1)/fs;

% Dodaj zakłócenie 50 Hz
f_noise = 50;
noise = 100 * sin(2*pi*f_noise*t);
ekg_noisy = ekg + noise;

% Parametry LMS
M = 12;                      % Długość filtru
mi = 0.01;                   % Współczynnik uczenia
d = ekg_noisy;               % Sygnał wyjściowy (z zakłóceniem)
x = sin(2*pi*f_noise*t);     % Referencja: czysty sinus 50 Hz

% Inicjalizacja
bx = zeros(M, 1);            % Bufor próbek x
h = zeros(M, 1);             % Wektory wag
y = zeros(size(x));          % Wyjście filtru
e = zeros(size(x));          % Błąd

% Implementacja filtru LMS
for n = 1:length(x)
    % Opóźnienie i buforowanie próbek x
    bx = [x(n); bx(1:end-1)];
    
    % Wyjście filtru
    y(n) = h' * bx;
    
    % Oblicz błąd: sygnał odniesienia minus wyjście filtru
    e(n) = d(n) - y(n);
    
    % Aktualizacja wag
    h = h + mi * e(n) * bx;
end

% Sygnał odszumiony
ekg_filtered = e;

% Wizualizacja
figure;
subplot(3,1,1);
plot(t, ekg);
title('Oryginalny sygnał EKG');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(3,1,2);
plot(t, ekg_noisy);
title('EKG z zakłóceniem 50 Hz');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(3,1,3);
plot(t, ekg_filtered);
title('EKG po własnym filtrze LMS');
xlabel('Czas [s]');
ylabel('Amplituda');

% Ocena błędu
mae = mean(abs(ekg - ekg_filtered));
fprintf('MAE = %.6f\n', mae);
