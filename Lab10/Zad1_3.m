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

    % === ODSŁUCH ===
    fprintf('\n=== Odtwarzanie: %s ===\n', files{idx});

    fprintf('Oryginał:\n');
    sound(y, fs_target);
    pause(length(y)/fs_target + 0.5);

    fprintf('Po kompresji (LPC):\n');
    sound(y_syn_all, fs_target);
    pause(length(y_syn_all)/fs_target + 0.5);

    % === DEKODER zgodnie ze schematem źródło–filtr ===
    y_decoded = zeros(size(y));
    for k = 1:numFrames
        startIdx = (k-1)*N + 1;
        a = lpc_coeffs(:, k);
        e = excitations(:, k);

        % === Blok: "Wzmocnienie" ===
        G = sqrt(sum(e.^2));  % Energia ekscytacji (lub gain z kompresji)
        e_scaled = e * G;

        % === Filtr traktu głosowego ===
        s = filter(1, a, e_scaled);
        y_decoded(startIdx:startIdx+N-1) = s;
        
    end

    fprintf('Po dekodowaniu (źródło–filtr):\n');
    sound(y_decoded, fs_target);
    pause(length(y_decoded)/fs_target + 0.5);
end
