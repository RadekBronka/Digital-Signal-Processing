clear; close all; clc;

% === PARAMETRY ===
fs_target = 8000;
frame_len = 100;        % długość ramki (20 ms)
frame_shift = 80;       % przesunięcie między ramkami (16 ms)
p = 10;                 % rząd LPC

% === Wczytanie pliku ===
file = 'dzwiecznadluga2.wav';
[y, fs] = audioread(file);
y = resample(y, fs_target, fs);
y = y(:);

% === Wybierz fragment do ekstrakcji resztki ===
ref_start = 200;  % dobierz ręcznie na podstawie sygnału
ref_frame = y(ref_start:ref_start+frame_len-1);

% === Preemfaza i okno ===
preemph = [1 -0.95];
ref_pre = filter(preemph, 1, ref_frame);
w = hamming(frame_len);
ref_win = ref_pre .* w;

% === LPC i resztka ===
a_ref = lpc(ref_win, p);
residual_ref = filter(a_ref, 1, ref_win);

% === Estymacja tonu podstawowego ===
y_lp = lowpass(ref_win, 900, fs_target);
thresh = 0.5 * max(abs(y_lp));
y_thr = y_lp;
y_thr(abs(y_lp) < thresh) = 0;
[r, lags] = xcorr(y_thr, 'coeff');
mid = ceil(length(r)/2);
startLag = 5;
searchRange = 80;
[~, rel_lag] = max(r(mid + startLag : mid + startLag + searchRange - 1));
pitch_period = rel_lag + (startLag - 1);

% === Pobudzenie: jeden okres ===
one_period = residual_ref(1:pitch_period);
e1_template = repmat(one_period, ceil(frame_len / pitch_period), 1);
e1_template = e1_template(1:frame_len);
gain1 = sqrt(sum(ref_win.^2) / sum(e1_template.^2));
e1_template = gain1 * e1_template;

% === Pobudzenie: średni okres ===
num_periods = 5;
if pitch_period * num_periods <= length(residual_ref)
    reshaped = reshape(residual_ref(1:(pitch_period*num_periods)), pitch_period, num_periods);
    avg_residual = mean(reshaped, 2);
else
    avg_residual = one_period;
end
e2_template = repmat(avg_residual, ceil(frame_len / pitch_period), 1);
e2_template = e2_template(1:frame_len);
gain2 = sqrt(sum(ref_win.^2) / sum(e2_template.^2));
e2_template = gain2 * e2_template;

% === Dekodowanie całego sygnału ===
L = length(y);
num_frames = floor((L - frame_len)/frame_shift) + 1;
y_syn1 = zeros(L, 1);
y_syn2 = zeros(L, 1);
w = hamming(frame_len);  % okno dla każdej ramki

for i = 1:num_frames
    idx_start = (i-1)*frame_shift + 1;
    idx_end = idx_start + frame_len - 1;
    
    frame = y(idx_start:idx_end);
    frame_pre = filter(preemph, 1, frame);
    frame_win = frame_pre .* w;
    
    a = lpc(frame_win, p);

    % --- Synteza z e1 (jeden okres)
    e1 = e1_template;
    gain1_frame = sqrt(sum(frame_win.^2) / sum(e1.^2));
    e1 = gain1_frame * e1;
    y_s1 = filter(1, a, e1);

    % --- Synteza z e2 (średni okres)
    e2 = e2_template;
    gain2_frame = sqrt(sum(frame_win.^2) / sum(e2.^2));
    e2 = gain2_frame * e2;
    y_s2 = filter(1, a, e2);

    % --- Overlap-add
    y_syn1(idx_start:idx_end) = y_syn1(idx_start:idx_end) + y_s1;
    y_syn2(idx_start:idx_end) = y_syn2(idx_start:idx_end) + y_s2;
end

% === Odsłuch i zapis ===
sound(y, fs_target);
pause(length(y)/fs_target + 0.5);
sound(y_syn1, fs_target);
pause(length(y)/fs_target + 0.5);
sound(y_syn2, fs_target);

audiowrite('syn1_jeden_okres.wav', y_syn1, fs_target);
audiowrite('syn2_sredni_okres.wav', y_syn2, fs_target);

% === Wykresy widmowe (opcjonalne) ===
Nfft = 512;
f = (0:Nfft-1)*(fs_target/Nfft);
Y = abs(fft(y(1:frame_len), Nfft));
Y1 = abs(fft(y_syn1(1:frame_len), Nfft));
Y2 = abs(fft(y_syn2(1:frame_len), Nfft));

figure;
plot(f, 20*log10(Y/max(Y)), 'k'); hold on;
plot(f, 20*log10(Y1/max(Y1)), 'b--');
plot(f, 20*log10(Y2/max(Y2)), 'r-.');
legend('Oryginał', 'Synteza: jeden okres', 'Synteza: średni okres');
xlabel('Częstotliwość [Hz]');
ylabel('Wzmocnienie [dB]');
title('Porównanie widm mowy (1 ramka)');
xlim([0 4000]);
grid on;
