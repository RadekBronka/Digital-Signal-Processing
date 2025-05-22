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




%%dodatkowe
emg_noise = randn(size(t));
fc = 150; % Częstotliwość graniczna
[b_emg, a_emg] = butter(4, fc/(fs/2), 'low'); % Filtr LP
emg_filtered = filter(b_emg, a_emg, emg_noise);
emg_filtered = 30 * emg_filtered; % Skalowanie

ekg_noisy_emg = ekg_noisy + emg_filtered;

% === 4. Filtr LMS - Interference Cancellation (IC) ===
M = 12;              % Długość filtru
mi = 0.0000001;           % Krok uczenia
d = ekg_noisy_emg;   % Sygnał z zakłóceniami
x = sin(2*pi*f_noise*t); % Referencja: znany 50 Hz

h = zeros(M, 1);
bx = zeros(M, 1);
y = zeros(size(x));
e = zeros(size(x));

for n = 1:length(x)
    bx = [x(n); bx(1:end-1)];
    y(n) = h' * bx;
    e(n) = d(n) - y(n);
    h = h + mi * e(n) * bx;
end

ekg_ic_filtered = e;

% === 5. Filtr LMS - Adaptive Linear Prediction (ALP) ===
h_pred = zeros(M, 1);
bx_pred = zeros(M, 1);
y_pred = zeros(size(ekg_noisy_emg));
e_pred = zeros(size(ekg_noisy_emg));

for n = 1:length(ekg_noisy_emg)
    if n > M
        bx_pred = ekg_noisy_emg(n-1:-1:n-M)';
        y_pred(n) = h_pred' * bx_pred;
        e_pred(n) = ekg_noisy_emg(n) - y_pred(n);
        h_pred = h_pred + mi * e_pred(n) * bx_pred;
    else
        e_pred(n) = ekg_noisy_emg(n);
    end
end

ekg_alp_filtered = e_pred;

% === 6. Ocena jakości - MAE ===
mae_ic = mean(abs(ekg - ekg_ic_filtered));
mae_alp = mean(abs(ekg - ekg_alp_filtered));

fprintf('MAE (Interference Cancellation) = %.6f\n', mae_ic);
fprintf('MAE (Linear Prediction)         = %.6f\n', mae_alp);

% === 7. Wykresy porównawcze ===
figure;

subplot(4,1,1);
plot(t, ekg); grid on;
title('Oryginalny sygnał EKG'); xlabel('Czas [s]'); ylabel('Amplituda');

subplot(4,1,2);
plot(t, ekg_noisy_emg); grid on;
title('EKG + zakłócenia: 50 Hz + szum EMG'); xlabel('Czas [s]'); ylabel('Amplituda');

subplot(4,1,3);
plot(t, ekg_ic_filtered); grid on;
title('Po filtrze LMS (Interference Cancellation)'); xlabel('Czas [s]'); ylabel('Amplituda');

subplot(4,1,4);
plot(t, ekg_alp_filtered); grid on;
title('Po filtrze LMS (Linear Prediction)'); xlabel('Czas [s]'); ylabel('Amplituda');
