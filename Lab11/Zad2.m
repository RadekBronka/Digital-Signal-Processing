clear; clc;

% Parametry i wczytanie sygnału
[x, fs] = audioread('DontWorryBeHappy.wav');
x = mean(x, 2); % mono
x = x(1:5*fs);  % tylko 5 sekund do testu

N_values = [32, 128];
Q_values = [0.1, 0.5, 1, 2, 4, 8, 16,32,64];  % testowane Q
colors = lines(length(Q_values));

figure;
for ni = 1:length(N_values)
    N = N_values(ni);
    L = N/2;
    h = sin(pi * ((0:N-1) + 0.5) / N)'; % okno

    % Macierz MDCT
    A = zeros(N/2, N);
    for k = 0:(N/2 - 1)
        for n = 0:(N - 1)
            A(k+1, n+1) = sqrt(4/N) * cos(2*pi/N * (k + 0.5) * (n + 0.5 + N/4));
        end
    end
    S = A'; % macierz syntezy

    % Liczba ramek
    numFrames = floor((length(x) - N) / L) + 1;

    snr_all = zeros(1, length(Q_values));
    
    for qi = 1:length(Q_values)
        Q = Q_values(qi);

        yq = zeros(N/2, numFrames);
        for i = 1:numFrames
            idx = (1:N) + (i-1)*L;
            frame = x(idx) .* h;
            X = A * frame;
            yq(:, i) = round(X * Q);  % kwantyzacja
        end

        % Dekodowanie
        xr = zeros(L*(numFrames+1), 1);
        for i = 1:numFrames
            X_hat = yq(:, i) / Q;     % dekwantyzacja
            x_rec = S * X_hat;        % odwrotna MDCT
            x_rec = x_rec .* h;       % okno
            idx = (1:N) + (i-1)*L;
            xr(idx) = xr(idx) + x_rec;
        end

        % Dopasowanie długości
        x_short = x(1:length(xr));
        xr = xr(1:length(x_short));

        % Oblicz SNR
        noise = x_short - xr;
        snr_val = 10*log10(sum(x_short.^2) / sum(noise.^2));
        snr_all(qi) = snr_val;
    end

    % Rysowanie wykresu
    subplot(1, 2, ni);
    plot(Q_values, snr_all, '-o', 'LineWidth', 2);
    xlabel('Q'); ylabel('SNR [dB]');
    title(['N = ', num2str(N)]);
    grid on;
end
sgtitle('Porównanie jakości (SNR) vs. Q dla różnych N');
