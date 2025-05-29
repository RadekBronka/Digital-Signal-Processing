clear; close all; clc;

file = 'mowa1.wav';
fs_target = 8000;
N = 160;               % długość ramki (20ms)
p = 10;                % rząd LPC
preemph = [1 -0.95];   % filtr preemfazy
max_len = fs_target * 4;

% === Wczytanie i przygotowanie sygnału ===
[y, fs] = audioread(file);
y = resample(y, fs_target, fs);
y = y(:);  % kolumna

if length(y) > max_len
    y = y(1:max_len);
end

% --- Preemfaza ---
y_pre = filter(preemph, 1, y);

% --- Inicjalizacja ---
numFrames = floor(length(y)/N);
y_syn_all = zeros(size(y));
lpc_coeffs = zeros(p+1, numFrames);
residuals = zeros(N, numFrames);
voiced_flags = false(1, numFrames);
pitch_periods = zeros(1, numFrames);

% === Kodowanie (LPC + resztkowy/szum) ===
for k = 1:numFrames
    idx = (k-1)*N + 1;
    frame = y(idx:idx+N-1);
    frame_pre = filter(preemph, 1, frame);
    frame_win = frame_pre .* hamming(N);

    a = lpc(frame_win, p);
    lpc_coeffs(:, k) = a(:);

    % --- Autokorelacja do detekcji dźwięczności ---
    [r, ~] = xcorr(frame_win, 'coeff');
    mid = ceil(length(r)/2);

    startLag = 5;
    searchRange = 80;
    [rmax, rel_lag] = max(r(mid + startLag : mid + startLag + searchRange - 1));
    pitch_period = rel_lag + (startLag - 1);
    pitch_periods(k) = pitch_period;

    if rmax > 0.3
        voiced = true;
    else
        voiced = false;
    end
    voiced_flags(k) = voiced;

    % --- Sygnał pobudzający ---
    if voiced
        e = filter(a, 1, frame_win);         % pełen sygnał resztkowy
    else
        rng(k);                              % dla deterministyczności
        e = randn(N, 1);
        gain = sqrt(sum(frame_win.^2) / sum(e.^2));  % dokładne dopasowanie energii
        e = gain * e;
    end
    residuals(:, k) = e;

    y_syn = filter(1, a, e);                 % synteza (dla porównania)
    y_syn_all(idx:idx+N-1) = y_syn;
end

% === Dekodowanie (rekonstrukcja z residualu/szumu) ===
y_decoded_raw = zeros(size(y));
y_decoded_deemph = zeros(size(y));

for k = 1:numFrames
    idx = (k-1)*N + 1;
    a = lpc_coeffs(:, k);
    e = residuals(:, k);

    s = filter(1, a, e);
    y_decoded_raw(idx:idx+N-1) = s;
    y_decoded_deemph(idx:idx+N-1) = filter(1, preemph, s);
end

% === Wykresy czasowe ===
figure('Name', 'Porównanie sygnałów - domena czasu', 'Position', [100 100 800 700]);

subplot(5,1,1);
plot(y); title('1. Oryginalny sygnał'); ylabel('Amplituda');

subplot(5,1,2);
plot(y_pre); title('2. Po preemfazie'); ylabel('Amplituda');

subplot(5,1,3);
plot(y_syn_all); title('3. Po kompresji (LPC + residual/szum)'); ylabel('Amplituda');

subplot(5,1,4);
plot(y_decoded_raw); title('4. Po dekodowaniu (bez deemfazy)'); ylabel('Amplituda');

subplot(5,1,5);
plot(y_decoded_deemph); title('5. Po dekodowaniu (z deemfazą)'); xlabel('Próbki'); ylabel('Amplituda');

% === Wykres dźwięczności i pitch period ===
figure('Name', 'Detekcja dźwięczności');
subplot(2,1,1);
stem(voiced_flags, 'filled');
title('Flaga dźwięczności (1 = voiced, 0 = unvoiced)');
ylabel('Voiced');
xlabel('Numer ramki');

subplot(2,1,2);
plot(pitch_periods);
title('Okres tonu podstawowego (tylko informacyjnie)');
ylabel('Pitch period [próbki]');
xlabel('Numer ramki');

% === Odsłuch ===
fprintf('\n=== Odtwarzanie: %s ===\n', file);

fprintf('Oryginał:\n');
sound(y, fs_target);
pause(length(y)/fs_target + 0.5);

fprintf('Po kompresji (LPC + residual/szum):\n');
sound(y_syn_all, fs_target);
pause(length(y_syn_all)/fs_target + 0.5);

fprintf('Po dekodowaniu (bez deemfazy):\n');
sound(y_decoded_raw, fs_target);
pause(length(y_decoded_raw)/fs_target + 0.5);

fprintf('Po dekodowaniu (z deemfazą):\n');
sound(y_decoded_deemph, fs_target);
pause(length(y_decoded_deemph)/fs_target + 0.5);
