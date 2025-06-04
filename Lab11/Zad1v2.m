[x, fs] = audioread('DontWorryBeHappy.wav');  % wczytanie próbki dźwiękowej
x = x(1:min(end, 4*fs), :);                   % przycięcie do 4 sekund
x = double(x);                                % konwersja na double

a = 0.9545;                                    % parametr a kodera
Y = zeros(size(x));                            % inicjalizacja macierzy na zrekonstruowany sygnał
MSE = zeros(1, size(x,2));                     % inicjalizacja MSE dla każdego kanału

% przetwarzanie każdego kanału osobno
for ch = 1:size(x,2)
    x_ch = x(:, ch);
    d = x_ch - a * [0; x_ch(1:end-1)];         % KODER
    dq = lab11_kwant(d);                      % kwantyzator

    % DEKODER
    y = zeros(size(dq));
    for n = 2:length(dq)
        y(n) = dq(n) + a * y(n-1);            % rekonstrukcja (dekoder)
    end
    Y(:, ch) = y;                             % zapisz zrekonstruowany kanał

    % oblicz MSE
    MSE(ch) = mean((x_ch - y).^2);
end

% WYKRES porównawczy obu kanałów
figure;
subplot(2,1,1);
plot(x(:,1), 'b'); hold on; plot(Y(:,1), 'r');
title(['Lewy kanał (MSE = ', num2str(MSE(1), '%.4g'), ')']);
legend('Oryg.', 'Rekonstr.');

subplot(2,1,2);
plot(x(:,2), 'b'); hold on; plot(Y(:,2), 'r');
title(['Prawy kanał (MSE = ', num2str(MSE(2), '%.4g'), ')']);
legend('Oryg.', 'Rekonstr.');

% ODSŁUCH połączonego sygnału stereo po rekonstrukcji
fprintf('Odtwarzam zrekonstruowany sygnał stereo...\n');
sound(Y, fs);

