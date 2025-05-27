clear; close all; clc;

% === PARAMETRY ===
fs_target = 8000;
N = 160; % długość ramki (20 ms)

% === NAZWY PLIKÓW ===
files = {'dzwieczne.wav', 'bezdzwieczne.wav', 'przejscie.wav'};
labels = {'Dźwięczna', 'Bezdźwięczna', 'Przejściowa'};

% === FILTR PREEMFAZY ===
preemph = [1 -0.95];

for idx = 1:3
    % === a) Wczytanie i preemfaza ===
    [y, fs] = audioread(files{idx});
    y = resample(y, fs_target, fs);   % resampling
    y = y(:);                          % upewnij się, że wektor kolumnowy
    y = y(1:N);                        % wybierz pierwszą ramkę (160 próbek)

    y_pre = filter(preemph, 1, y);     % preemfaza

    % Wykresy (a) sygnał i widmo
    figure('Name', ['(a) ' labels{idx}]);
    subplot(2,2,1); plot(y); title('Sygnał oryginalny'); xlabel('Próbka');
    subplot(2,2,2); pwelch(y, [], [], [], fs_target); title('Widmo oryginalne');

    subplot(2,2,3); plot(y_pre); title('Po preemfazie'); xlabel('Próbka');
    subplot(2,2,4); pwelch(y_pre, [], [], [], fs_target); title('Widmo po preemfazie');

    % === Okno Hamminga ===
    w = hamming(N);
    yw = y_pre .* w;

    % === Filtr LP ===
    y_lp = lowpass(yw, 900, fs_target);

    % === c) Progowanie ===
    thresh = 0.5 * max(abs(y_lp));
    y_thr = y_lp;
    y_thr(abs(y_lp) < thresh) = 0;

    figure('Name', ['(c) ' labels{idx}]);
    subplot(2,1,1); plot(y_lp); title('Po filtrze LP'); xlabel('Próbka');
    subplot(2,1,2); plot(y_thr); title('Po progowaniu'); xlabel('Próbka');

    % === d) Autokorelacja ===
    [r, lags] = xcorr(y_thr, 'coeff');
    mid = ceil(length(r)/2);

    figure('Name', ['(d) Autokorelacja - ' labels{idx}]);
    plot(lags, r); hold on;
    yline(0.3, '--r'); title(['Autokorelacja ' labels{idx}]);
    xlabel('Lagi'); ylabel('Wartość');

   % === e) Decyzja o dźwięczności ===
    startLag = 5; % ignorujemy lagi 1,2,3,4
    searchRange = 80;

    % Szukamy maksimum w zakresie lagów od startLag do startLag+searchRange-1
    [rmax, rel_lag] = max(r(mid + startLag : mid + startLag + searchRange - 1));
    pitch_period = rel_lag + (startLag - 1); % poprawiamy lag o startLag

    f0 = fs_target / pitch_period;

if rmax > 0.3
    fprintf('%s: DŹWIĘCZNA, F0 ≈ %.2f Hz\n', labels{idx}, f0);
    voiced = true;
else
    fprintf('%s: BEZDŹWIĘCZNA\n', labels{idx});
    voiced = false;
end
    % === b) Estymacja współczynników LPC ===
    p = 10;
    a = lpc(yw, p);
    [H, f] = freqz(1, a, 512, fs_target);

    figure('Name', ['(b) |H(f)| - ' labels{idx}]);
    plot(f, abs(H)); title('|H(f)| - Filtr LPC');
    xlabel('Częstotliwość [Hz]'); ylabel('Amplituda');

    % === f) Synteza ramki ===
    if voiced
        e = zeros(N,1);
        e(1:pitch_period:N) = 1;
    else
        e = randn(N,1);
    end

    % Skalowanie energii
    gain = sqrt(sum(yw.^2)/sum(e.^2));
    e = gain * e;

    y_syn = filter(1, a, e);

    figure('Name', ['(f) Synteza - ' labels{idx}]);
    subplot(2,2,1); plot(yw); title('Oryginał (czas)');
    subplot(2,2,2); pwelch(yw, [], [], [], fs_target); title('Oryginał (widmo)');
    subplot(2,2,3); plot(y_syn); title('Synteza (czas)');
    subplot(2,2,4); pwelch(y_syn, [], [], [], fs_target); title('Synteza (widmo)');

    
end

