function [Y, MSE] = g726_adpcm(x, fs, bitrate)
    a = 0.9545;  % współczynnik predyktora
    Y = zeros(size(x));
    MSE = zeros(1, size(x,2));

    % Liczba bitów na próbkę w zależności od bitrate
    switch bitrate
        case 16000
            nbits = 2;
        case 32000
            nbits = 4;
        otherwise
            error('Obsługiwane tylko 16 kbps i 32 kbps');
    end

    levels = 2^nbits;
    delta = 2 * max(abs(x(:))) / levels;  % krok kwantyzacji

    for ch = 1:size(x,2)
        x_ch = x(:, ch);
        d = x_ch - a * [0; x_ch(1:end-1)];

        % KWANTYZACJA (różnicy)
        q_idx = round(d / delta);
        q_idx = max(min(q_idx, levels/2 - 1), -levels/2);  % ograniczenie zakresu

        % DEKWANTYZACJA
        dq = q_idx * delta;

        % DEKODER
        y = zeros(size(dq));
        for n = 2:length(dq)
            y(n) = dq(n) + a * y(n-1);
        end

        Y(:, ch) = y;
        MSE(ch) = mean((x_ch - y).^2);
    end
end
