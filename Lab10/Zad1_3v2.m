clear; close all; clc;

files = {'mowa1.wav', 'mowa2.wav', 'mowa3.wav'};
fs_target = 8000;
N = 160;  % długość ramki 20ms
p = 10;   % rząd LPC
preemph = [1 -0.95];
max_len = fs_target * 4;  % 4 sekundy = 32000 próbek

for idx = 1:length(files)
    [y, fs] = audioread(files{idx});
    y = resample(y, fs_target, fs);
    y = y(:);  % upewniamy się, że kolumna

    % Przycięcie do 4 sekund
    if length(y) > max_len
        y = y(1:max_len);
    end

    numFrames = floor(length(y)/N);
    y_syn_all = zeros(size(y));

    % Bufory do "skompresowanych" danych (symulacja)
    lpc_coeffs = zeros(p+1, numFrames);
    excitations = zeros(N, numFrames);

    for k = 1:numFrames
        startIdx = (k-1)*N + 1;
        frame = y(startIdx:startIdx+N-1);
        frame_pre = filter(preemph, 1, frame);
        frame_win = frame_pre .* hamming(N);

        % Estymacja LPC
        a = lpc(frame_win, p);
        lpc_coeffs(:, k) = a(:);

        % Źródło ekscytacji
        y_lp = lowpass(frame_win, 900, fs_target);
        thresh = 0.5 * max(abs(y_lp));
        y_thr = y_lp;
        y_thr(abs(y_lp) < thresh) = 0;

        [r, ~] = xcorr(y_thr, 'coeff');
        mid = ceil(length(r)/2);
        startLag = 5;
        searchRange = 80;

        [rmax, rel_lag] = max(r(mid+startLag:mid+startLag+searchRange-1));
        pitch_period = rel_lag + (startLag - 1);

        if rmax > 0.3
            e = zeros(N,1);
            e(1:pitch_period:N) = 1;
        else
            e = randn(N,1);
        end

        % Dopasowanie energii
        gain = sqrt(sum(frame_win.^2)/sum(e.^2));
        e = gain * e;
        excitations(:, k) = e;

        % Synteza
        y_syn = filter(1, a, e);
        y_syn_all(startIdx:startIdx+N-1) = y_syn;
    end
    
    time_axis = (0:length(y)-1)/fs_target;  % oś czasu w sekundach

    figure;
    subplot(2,1,1);
    plot(time_axis, y);
    title('Sygnał oryginalny (cały)');
    xlabel('Czas [s]');
    ylabel('Amplituda');
    xlim([0 time_axis(end)]);  % pokazujemy cały zakres czasu

    subplot(2,1,2);
    plot(time_axis, y_syn_all);
    title('Sygnał po syntezie LPC (cały)');
    xlabel('Czas [s]');
    ylabel('Amplituda');
    xlim([0 time_axis(end)]);
    % FFT
Nfft = 2^nextpow2(length(y));  % liczba punktów FFT, potęga 2 >= długość sygnału

Y_orig = fft(y, Nfft);
Y_syn = fft(y_syn_all, Nfft);

f = fs_target/2*linspace(0,1,Nfft/2+1);  % oś częstotliwości od 0 do Nyquista

figure;
subplot(2,1,1);
plot(f, 20*log10(abs(Y_orig(1:Nfft/2+1))));
title('Widmo amplitudowe sygnału oryginalnego');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
grid on;

subplot(2,1,2);
plot(f, 20*log10(abs(Y_syn(1:Nfft/2+1))));
title('Widmo amplitudowe sygnału po syntezie LPC');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
grid on;
    % === ODSŁUCH ===
    fprintf('\n=== Odtwarzanie: %s ===\n', files{idx});

    fprintf('Oryginał:\n');
    sound(y, fs_target);
    pause(length(y)/fs_target + 0.5);

    fprintf('Po kompresji (LPC):\n');
    sound(y_syn_all, fs_target);
    pause(length(y_syn_all)/fs_target + 0.5);

   
end
