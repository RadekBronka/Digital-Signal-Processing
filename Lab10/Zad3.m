clear; close all; clc;

file = 'mowa1.wav';
fs_target = 8000;
N = 160;  % długość ramki 20ms
p = 10;   % rząd LPC
preemph = [1 -0.95];
max_len = fs_target * 4;

[y, fs] = audioread(file);
y = resample(y, fs_target, fs);
y = y(:);  % kolumna

% Przycięcie do 4 sekund
if length(y) > max_len
    y = y(1:max_len);
end

% --- Preemfaza ---
y_pre = filter(preemph, 1, y);

% Inicjalizacja
numFrames = floor(length(y)/N);
y_syn_all = zeros(size(y));
lpc_coeffs = zeros(p+1, numFrames);
residuals = zeros(N, numFrames);

% === Kodowanie (analiza LPC i residual) ===
for k = 1:numFrames
    idx = (k-1)*N + 1;
    frame = y(idx:idx+N-1);
    frame_pre = filter(preemph, 1, frame);
    frame_win = frame_pre .* hamming(N);

    a = lpc(frame_win, p);
    lpc_coeffs(:, k) = a(:);

    e = filter(a, 1, frame_win);
    residuals(:, k) = e;

    y_syn = filter(1, a, e);
    y_syn_all(idx:idx+N-1) = y_syn;
end

% === Dekodowanie ===
y_decoded_raw = zeros(size(y));      % bez deemfazy
y_decoded_deemph = zeros(size(y));   % z deemfazą

for k = 1:numFrames
    idx = (k-1)*N + 1;
    a = lpc_coeffs(:, k);
    e = residuals(:, k);

    s = filter(1, a, e);
    y_decoded_raw(idx:idx+N-1) = s;
    y_decoded_deemph(idx:idx+N-1) = filter(1, preemph, s);  % deemfaza
end

% === Wykresy czasowe ===
figure('Name', 'Porównanie sygnałów - domena czasu', 'Position', [100 100 800 700]);

subplot(5,1,1);
plot(y);
title('1. Oryginalny sygnał');
ylabel('Amplituda');

subplot(5,1,2);
plot(y_pre);
title('2. Po preemfazie');
ylabel('Amplituda');

subplot(5,1,3);
plot(y_syn_all);
title('3. Po kompresji (synteza LPC)');
ylabel('Amplituda');

subplot(5,1,4);
plot(y_decoded_raw);
title('4. Po dekodowaniu (bez deemfazy)');
ylabel('Amplituda');

subplot(5,1,5);
plot(y_decoded_deemph);
title('5. Po dekodowaniu (z deemfazą)');
xlabel('Próbki'); ylabel('Amplituda');

% === Odsłuch ===
fprintf('\n=== Odtwarzanie: %s ===\n', file);

fprintf('Oryginał:\n');
sound(y, fs_target);
pause(length(y)/fs_target + 0.5);

fprintf('Po kompresji (LPC z residualem):\n');
sound(y_syn_all, fs_target);
pause(length(y_syn_all)/fs_target + 0.5);

fprintf('Po dekodowaniu (bez deemfazy):\n');
sound(y_decoded_raw, fs_target);
pause(length(y_decoded_raw)/fs_target + 0.5);

fprintf('Po dekodowaniu (z deemfazą):\n');
sound(y_decoded_deemph, fs_target);
pause(length(y_decoded_deemph)/fs_target + 0.5);
