clear; close all; clc;

% === PARAMETRY ===
fs_target = 8000;
N = 160; % długość ramki (20 ms)
p = 10;
file = 'dzwiecznadluga2.wav';

% === Wczytanie pliku dźwięcznego ===
[y, fs] = audioread(file);
y = resample(y, fs_target, fs);
y = y(:);
len = length(y);

% === Preemfaza całego sygnału ===
preemph = [1 -0.95];
y_pre = filter(preemph, 1, y);

% === Okno i parametry ===
w = hamming(N);
num_avg_frames = 5;

% === Zbieranie 5 ramek do wyliczenia średniego residuala ===
combined_residual = [];

for i = 1:num_avg_frames
    idx = (i-1)*N + (1:N);
    if idx(end) > length(y_pre)
        break
    end
    y_frame = y_pre(idx);
    y_win = y_frame .* w;
    a = lpc(y_win, p);
    r = filter(a, 1, y_win);
    combined_residual = [combined_residual; r];
end

% === Wyznaczenie jednego okresu z pierwszej ramki ===
y_ref = y_pre(1:N);
y_ref_win = y_ref .* w;
a_ref = lpc(y_ref_win, p);
residual_ref = filter(a_ref, 1, y_ref_win);

y_lp = lowpass(y_ref_win, 900, fs_target);
thresh = 0.5 * max(abs(y_lp));
y_thr = y_lp;
y_thr(abs(y_lp) < thresh) = 0;
[r, ~] = xcorr(y_thr, 'coeff');
mid = ceil(length(r)/2);
startLag = 5;
searchRange = 80;
[~, rel_lag] = max(r(mid+startLag : mid+startLag+searchRange-1));
pitch_period = rel_lag + (startLag - 1);

% === Jeden okres pobudzenia ===
one_period = residual_ref(1:pitch_period);
num_repeats = ceil(N / pitch_period);
e1_template = repmat(one_period, num_repeats, 1);
e1_template = e1_template(1:N);

% === Średni sygnał resztkowy z 5 ramek ===
if pitch_period * num_avg_frames <= length(combined_residual)
    periods_matrix = reshape(combined_residual(1:(pitch_period*num_avg_frames)), pitch_period, num_avg_frames);
    avg_residual = mean(periods_matrix, 2);
else
    avg_residual = one_period;
end

% --- Dla wykresu porównawczego
y_plot = y(1:N);
y_pre_plot = y_pre(1:N);
yw = y_pre_plot .* w;
a = lpc(yw, p);
residual = filter(a, 1, yw);

% === Synteza ===
num_frames = floor(len / N);
y_syn1 = zeros(len,1);
y_syn2 = zeros(len,1);

for k = 1:num_frames
    idx = (k-1)*N + (1:N);
    if idx(end) > len
        break
    end
    y_frame = y(idx);
    y_frame_pre = filter(preemph, 1, y_frame);
    yw = y_frame_pre .* hamming(N);
    a = lpc(yw, p);

    % Estymacja dźwięczności
    y_lp = lowpass(yw, 900, fs_target);
    thresh = 0.5 * max(abs(y_lp));
    y_thr = y_lp;
    y_thr(abs(y_lp) < thresh) = 0;
    [r, ~] = xcorr(y_thr, 'coeff');
    mid = ceil(length(r)/2);
    [rmax, ~] = max(r(mid+5:mid+80));
    is_voiced = rmax > 0.3;

    if is_voiced
        % --- Jeden okres
        e1 = e1_template;
        gain1 = sqrt(sum(yw.^2) / sum(e1.^2));
        e1 = gain1 * e1;
        y_syn1(idx) = filter(1, a, e1);

        % --- Średni okres
        e2 = repmat(avg_residual, ceil(N/pitch_period), 1);
        e2 = e2(1:N);
        gain2 = sqrt(sum(yw.^2) / sum(e2.^2));
        e2 = gain2 * e2;
        y_syn2(idx) = filter(1, a, e2);
    else
        % --- Biały szum
        e = randn(N,1);
        gain = sqrt(sum(yw.^2) / sum(e.^2));
        e = gain * e;
        y_syn1(idx) = filter(1, a, e);
        y_syn2(idx) = filter(1, a, e);
    end
end

% === Wykresy ===
figure;

subplot(4,1,1);
plot(y_plot);
title('Sygnał po imporcie (przed preemfazą i oknem)');
xlabel('Próbka');
ylabel('Amplituda');

subplot(4,1,2);
plot(yw);
title('Oryginalny sygnał (po preemfazie i oknie)');
xlabel('Próbka');
ylabel('Amplituda');

subplot(4,1,3);
plot(y_syn1(1:N));
title('Synteza z pobudzeniem: jeden okres sygnału resztkowego');
xlabel('Próbka');
ylabel('Amplituda');

subplot(4,1,4);
plot(y_syn2(1:N));
title('Synteza z pobudzeniem: średni sygnał resztkowy z kilku okresów');
xlabel('Próbka');
ylabel('Amplituda');
grid on;

% === Widma pobudzeń ===
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

% === Widma porównawcze ===
YW = abs(fft(yw, Nfft));
YS1 = abs(fft(y_syn1(1:N), Nfft));
YS2 = abs(fft(y_syn2(1:N), Nfft));
figure;
plot(f, 20*log10(YW/max(YW)), 'k'); hold on;
plot(f, 20*log10(YS1/max(YS1)), 'b--');
plot(f, 20*log10(YS2/max(YS2)), 'r-.');
legend('Oryginał', 'Synteza: jeden okres', 'Synteza: średni okres');
title('Porównanie widm sygnałów');
xlabel('Częstotliwość [Hz]');
ylabel('Wzmocnienie [dB]');
xlim([0 4000]);
grid on;
