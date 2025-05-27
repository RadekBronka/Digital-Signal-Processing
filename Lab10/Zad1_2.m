clear; close all; clc;

% === PARAMETRY ===
fs_target = 8000;
N = 160; % długość ramki (20 ms)
preemph = [1 -0.95];

% === PLIKI ===
files = {'dzwieczne.wav', 'bezdzwieczne.wav', 'przejscie.wav'};
labels = {'Dźwięczna', 'Bezdźwięczna', 'Przejściowa'};

% === RÓŻNE ILOŚCI BIEGUNÓW ===
pole_orders = [10, 8, 6, 4, 2];

% === WYNIKI BŁĘDÓW ===
errors = zeros(length(files), length(pole_orders));

for p_idx = 1:length(pole_orders)
    p = pole_orders(p_idx);
    fprintf('\n==== LICZBA BIEGUNÓW: %d ====\n', p);

    for idx = 1:3
        % === a) Wczytanie i preemfaza ===
        [y, fs] = audioread(files{idx});
        y = resample(y, fs_target, fs);
        y = y(:);
        y = y(1:N);
        y_pre = filter(preemph, 1, y);

        % === Okno Hamminga ===
        w = hamming(N);
        yw = y_pre .* w;

        % === Filtr LP ===
        y_lp = lowpass(yw, 900, fs_target);

        % === Progowanie ===
        thresh = 0.5 * max(abs(y_lp));
        y_thr = y_lp;
        y_thr(abs(y_lp) < thresh) = 0;

        % === Autokorelacja ===
        [r, lags] = xcorr(y_thr, 'coeff');
        mid = ceil(length(r)/2);
        startLag = 5;
        searchRange = 80;
        [rmax, rel_lag] = max(r(mid + startLag : mid + startLag + searchRange - 1));
        pitch_period = rel_lag + (startLag - 1);
        f0 = fs_target / pitch_period;

        if rmax > 0.3
            voiced = true;
        else
            voiced = false;
        end

        % === Estymacja LPC ===
        a = lpc(yw, p);

        % === Synteza ===
        if voiced
            e = zeros(N,1);
            e(1:pitch_period:N) = 1;
        else
            e = randn(N,1);
        end

        gain = sqrt(sum(yw.^2)/sum(e.^2));
        e = gain * e;

        y_syn = filter(1, a, e);

        % === Obliczanie błędu średniokwadratowego ===
        mse = mean((yw - y_syn).^2);
        errors(idx, p_idx) = mse;

        fprintf('%s, p = %d: Błąd MSE = %.10f\n', labels{idx}, p, mse);
    end
end

% === TABELA WYNIKÓW ===
fprintf('\n=== PODSUMOWANIE BŁĘDÓW MSE ===\n');
fprintf('             |   p=10         p=8          p=6          p=4          p=2\n');
for idx = 1:3
    fprintf('%-12s|', labels{idx});
    fprintf(' %.10f', errors(idx, :));
    fprintf('\n');
end

% === WYKRESY ===
figure;
for idx = 1:3
    subplot(3,1,idx);
    plot(pole_orders, errors(idx,:), '-o', 'LineWidth', 2);
    title(['Błąd MSE dla: ' labels{idx}]);
    xlabel('Liczba biegunów (p)');
    ylabel('MSE');
    grid on;
end
sgtitle('Porównanie błędów MSE dla różnych liczby biegunów');
