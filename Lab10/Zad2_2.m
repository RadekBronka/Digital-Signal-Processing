clear; close all; clc;

% === PARAMETRY ===
fs_target = 8000;
N = 160; % długość ramki (20 ms)

% === Wczytanie pliku dźwięcznego (głoska dźwięczna) ===
file = 'dzwieczne.wav';
[y, fs] = audioread(file);
y = resample(y, fs_target, fs);
y = y(:);
y = y(1:N);

% Preemfaza i okno Hamminga
preemph = [1 -0.95];
y_pre = filter(preemph, 1, y);
w = hamming(N);
yw = y_pre .* w;

% LPC
p = 10;
a = lpc(yw, p);

% Estymacja tonu podstawowego przez autokorelację
y_lp = lowpass(yw, 900, fs_target);
thresh = 0.5 * max(abs(y_lp));
y_thr = y_lp;
y_thr(abs(y_lp) < thresh) = 0;

[r, lags] = xcorr(y_thr, 'coeff');
mid = ceil(length(r)/2);
startLag = 5;
searchRange = 80;
[rmax, rel_lag] = max(r(mid + startLag : mid + startLag + searchRange - 1));
pitch_period = rel_lag + (startLag - 1);

% Obliczenie sygnału resztkowego
% Transmitancja H(z) = 1 / A(z), więc odwrotna to A(z) / 1
b_inv = 1;       % licznik odwrotnego filtru
a_inv = a;       % mianownik odwrotnego filtru (czyli A(z))

residual = filter(a_inv, b_inv, yw); % przefiltrowanie odwrotnym filtrem

% --- Pobudzenie 1 okres ---
one_period = residual(1:pitch_period);
e1 = zeros(N, 1);
for i = 0:floor(N/pitch_period)-1
    idx = i * pitch_period + 1;
    if idx + pitch_period - 1 <= N
        e1(idx:idx+pitch_period-1) = one_period;
    end
end
gain1 = sqrt(sum(yw.^2)/sum(e1.^2));
e1 = gain1 * e1;
y_syn1 = filter(1, a, e1);

% --- Pobudzenie średni sygnał z kilku okresów ---
num_periods = 5;
if pitch_period * num_periods <= length(residual)
    periods_matrix = reshape(residual(1:(pitch_period*num_periods)), pitch_period, num_periods);
    avg_residual = mean(periods_matrix, 2);
else
    avg_residual = one_period;
end

e2 = zeros(N, 1);
for i = 0:floor(N/pitch_period)-1
    idx = i * pitch_period + 1;
    if idx + pitch_period - 1 <= N
        e2(idx:idx+pitch_period-1) = avg_residual;
    end
end
gain2 = sqrt(sum(yw.^2)/sum(e2.^2));
e2 = gain2 * e2;
y_syn2 = filter(1, a, e2);

% --- Wykresy ---
figure;
subplot(3,1,1);
plot(yw);
title('Oryginalny sygnał (po preemfazie i oknie)');
xlabel('Próbka');
ylabel('Amplituda');

subplot(3,1,2);
plot(y_syn1);
title('Synteza z pobudzeniem: jeden okres sygnału resztkowego');
xlabel('Próbka');
ylabel('Amplituda');

subplot(3,1,3);
plot(y_syn2);
title('Synteza z pobudzeniem: średni sygnał resztkowy z kilku okresów');
xlabel('Próbka');
ylabel('Amplituda');

figure;
plot(one_period, 'b', 'LineWidth', 1.5); hold on;
plot(avg_residual, 'r--', 'LineWidth', 1.5);
legend('Jeden okres', 'Średni sygnał z kilku okresów');
title('Porównanie okresów sygnału resztkowego');
xlabel('Próbka');
ylabel('Amplituda');
grid on;

figure;
Nfft = 512;
f = (0:Nfft-1)*(fs_target/Nfft);

E1 = abs(fft(one_period, Nfft));
E2 = abs(fft(avg_residual, Nfft));

plot(f, 20*log10(E1/max(E1)), 'b', 'LineWidth', 1.5); hold on;
plot(f, 20*log10(E2/max(E2)), 'r--', 'LineWidth', 1.5);
legend('Jeden okres', 'Średni sygnał');
title('Widmo pobudzenia (sygnału resztkowego)');
xlabel('Częstotliwość [Hz]');
ylabel('Wzmocnienie [dB]');
xlim([0 4000]);
grid on;
