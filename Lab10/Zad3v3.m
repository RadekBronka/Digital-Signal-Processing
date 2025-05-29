clear; close all; clc;

files = {'mowa1.wav', 'mowa2.wav', 'mowa3.wav'};
fs_target = 8000;
N = 160;               % długość ramki (20ms)
p = 10;                % rząd LPC
preemph = [1 -0.95];
max_len = fs_target * 4;

for fileIdx = 1:length(files)
    file = files{fileIdx};
    fprintf('\nPrzetwarzanie pliku: %s\n', file);

    % === Wczytanie i przygotowanie sygnału ===
    [y, fs] = audioread(file);
    y = resample(y, fs_target, fs);
    y = y(:);
    if length(y) > max_len
        y = y(1:max_len);
    end

    % Preemfaza
    y_pre = filter(preemph, 1, y);

    % Inicjalizacja
    numFrames = floor(length(y)/N);
    y_syn_all = zeros(size(y));
    lpc_coeffs = zeros(p+1, numFrames);
    residuals = zeros(N, numFrames);
    voiced_flags = false(1, numFrames);
    pitch_periods = zeros(1, numFrames);

    % Kodowanie LPC + resztkowy/szum
    for k = 1:numFrames
        idx = (k-1)*N + 1;
        frame = y(idx:idx+N-1);
        frame_pre = filter(preemph, 1, frame);
        frame_win = frame_pre .* hamming(N);

        a = lpc(frame_win, p);
        lpc_coeffs(:, k) = a(:);

        % Autokorelacja do detekcji dźwięczności
        [r, ~] = xcorr(frame_win, 'coeff');
        mid = ceil(length(r)/2);

        startLag = 5;
        searchRange = 80;
        [rmax, rel_lag] = max(r(mid + startLag : mid + startLag + searchRange - 1));
        pitch_period = rel_lag + (startLag - 1);
        pitch_periods(k) = pitch_period;

        voiced = (rmax > 0.3);
        voiced_flags(k) = voiced;

        % Sygnał pobudzający
        if voiced
            e = filter(a, 1, frame_win);
        else
            rng(k); % deterministycznie
            e = randn(N, 1);
            gain = sqrt(sum(frame_win.^2) / sum(e.^2));
            e = gain * e;
        end
        residuals(:, k) = e;

        % Synteza (do porównania)
        y_syn = filter(1, a, e);
        y_syn_all(idx:idx+N-1) = y_syn;
    end

    % Dekodowanie
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

    % Wykresy czasowe
    figure('Name', ['Porównanie sygnałów - ' file], 'Position', [100 100 800 700]);
    subplot(5,1,1); plot(y); title('1. Oryginalny sygnał'); ylabel('Amplituda');
    subplot(5,1,2); plot(y_pre); title('2. Po preemfazie'); ylabel('Amplituda');
    subplot(5,1,3); plot(y_syn_all); title('3. Po kompresji (LPC + residual/szum)'); ylabel('Amplituda');
    subplot(5,1,4); plot(y_decoded_raw); title('4. Po dekodowaniu (bez deemfazy)'); ylabel('Amplituda');
    subplot(5,1,5); plot(y_decoded_deemph); title('5. Po dekodowaniu (z deemfazą)'); xlabel('Próbki'); ylabel('Amplituda');

    % Wykres dźwięczności i pitch period
    figure('Name', ['Detekcja dźwięczności - ' file]);
    subplot(2,1,1);
    stem(voiced_flags, 'filled');
    title('Flaga dźwięczności (1 = voiced, 0 = unvoiced)');
    ylabel('Voiced'); xlabel('Numer ramki');
    subplot(2,1,2);
    plot(pitch_periods);
    title('Okres tonu podstawowego');
    ylabel('Pitch period [próbki]'); xlabel('Numer ramki');

    % Porównanie widma oryginału i syntezy LPC
    segment_len = fs_target * 1; % 1 sekunda
    y_orig_seg = y(1:segment_len);
    y_syn_seg = y_syn_all(1:segment_len);

    NFFT = 2^nextpow2(segment_len);
    Y_orig = fft(y_orig_seg, NFFT);
    Y_syn = fft(y_syn_seg, NFFT);
    f = fs_target/2 * linspace(0,1,NFFT/2+1);

    figure('Name',['Widmo sygnałów - ' file]);
    plot(f, 20*log10(abs(Y_orig(1:NFFT/2+1)) + eps), 'b', 'LineWidth', 1.5); hold on;
    plot(f, 20*log10(abs(Y_syn(1:NFFT/2+1)) + eps), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Częstotliwość [Hz]');
    ylabel('Amplituda [dB]');
    title(['Widmo oryginału i syntezy LPC - ' file]);
    legend('Oryginał', 'Synteza LPC + residual/szum');

    % Odtwarzanie
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

end
