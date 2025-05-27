clear; close all; clc;

% Parametry
file = 'mowa1.wav';
fs_target = 8000;
N = 160;  % ramka 20ms
p = 10;   % rząd LPC
preemph = [1 -0.95];
max_len = fs_target * 3;  % 3 sekundy

% Wczytanie i przygotowanie sygnału
[y, fs] = audioread(file);
y = resample(y, fs_target, fs);
y = y(:);
y = y(1:min(end, max_len));

numFrames = floor(length(y)/N);
y = y(1:numFrames*N);  % przycinamy do pełnych ramek

% Bufory
lpc_coeffs = zeros(p+1, numFrames);
excitations = zeros(N, numFrames);
VUV = zeros(1, numFrames);
pitch_periods = zeros(1, numFrames);

% === ANALIZA ===
for k = 1:numFrames
    idx = (k-1)*N + 1;
    frame = y(idx:idx+N-1);
    frame_pre = filter(preemph, 1, frame);
    frame_win = frame_pre .* hamming(N);

    a = lpc(frame_win, p);
    lpc_coeffs(:, k) = a;

    y_lp = lowpass(frame_win, 900, fs_target);
    thresh = 0.5 * max(abs(y_lp));
    y_thr = y_lp;
    y_thr(abs(y_lp) < thresh) = 0;

    [r, ~] = xcorr(y_thr, 'coeff');
    mid = ceil(length(r)/2);
    [rmax, rel_lag] = max(r(mid+5:mid+85));
    pitch = rel_lag + 4;

    if rmax > 0.3
        VUV(k) = 1;
        e = zeros(N,1);
        e(1:pitch:N) = 1;
        pitch_periods(k) = pitch;
    else
        VUV(k) = 0;
        e = randn(N,1);
    end

    gain = sqrt(sum(frame_win.^2)/sum(e.^2));
    e = gain * e;
    excitations(:, k) = e;
end

% === DEKODOWANIE – 4 WERSJE ===
y1 = zeros(size(y));  % tylko szum
y2 = zeros(size(y));  % pitch / 2
y3 = zeros(size(y));  % pitch = 80
y4 = zeros(size(y));  % coldvox

[coldvox, fs_cold] = audioread('coldvox.wav');
coldvox = resample(coldvox, fs_target, fs_cold);
coldvox = coldvox(:);
cv_idx = 1;

for k = 1:numFrames
    idx = (k-1)*N + 1;
    a = lpc_coeffs(:, k);
    G = sqrt(sum(excitations(:,k).^2));

    % Eksperyment 1 – ignoruj V/UV, tylko szum
    e1 = randn(N,1) * G;
    y1(idx:idx+N-1) = filter(1, a, e1);

    % Eksperyment 2 – pitch / 2
    if VUV(k)
        T = pitch_periods(k) * 2;
        e2 = zeros(N,1);
        e2(1:T:end) = 1;
    else
        e2 = randn(N,1);
    end
    e2 = e2 * G;
    y2(idx:idx+N-1) = filter(1, a, e2);

    % Eksperyment 3 – stały pitch = 80
    if VUV(k)
        T = 80;
        e3 = zeros(N,1);
        e3(1:T:end) = 1;
    else
        e3 = randn(N,1);
    end
    e3 = e3 * G;
    y3(idx:idx+N-1) = filter(1, a, e3);

    % Eksperyment 4 – coldvox jako szum
    if cv_idx+N-1 <= length(coldvox)
        e4 = coldvox(cv_idx:cv_idx+N-1);
    else
        e4 = zeros(N,1);
    end
    cv_idx = cv_idx + N;
    e4 = e4 * G;
    y4(idx:idx+N-1) = filter(1, a, e4);
end

% === WYKRESY PORÓWNAWCZE ===
t = (0:length(y)-1)/fs_target;
figure;
subplot(5,1,1); plot(t, y); title('Oryginał');
subplot(5,1,2); plot(t, y1); title('Eksperyment 1: tylko szum');
subplot(5,1,3); plot(t, y2); title('Eksperyment 2: pitch/2');
subplot(5,1,4); plot(t, y3); title('Eksperyment 3: pitch = 80');
subplot(5,1,5); plot(t, y4); title('Eksperyment 4: coldvox jako pobudzenie');
xlabel('Czas [s]');

% === ODSŁUCH ===
fprintf('\n=== Odtwarzanie ===\n');
fprintf('Oryginał:\n'); sound(y, fs_target); pause(length(y)/fs_target + 0.5);
fprintf('Eksperyment 1:\n'); sound(y1, fs_target); pause(length(y)/fs_target + 0.5);
fprintf('Eksperyment 2:\n'); sound(y2, fs_target); pause(length(y)/fs_target + 0.5);
fprintf('Eksperyment 3:\n'); sound(y3, fs_target); pause(length(y)/fs_target + 0.5);
fprintf('Eksperyment 4:\n'); sound(y4, fs_target); pause(length(y)/fs_target + 0.5);
