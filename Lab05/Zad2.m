omega3dB = 2*pi*100;   % pulsacja odcięcia [rad/s]
orders = [2, 4, 6, 8]; % rzędy filtrów

w = linspace(0, 2*omega3dB, 1000); % zakres pulsacji
f = w / (2*pi);                   % częstotliwość [Hz]

%% CHARAKTERYSTYKA AMPLITUDOWA – skala liniowa
figure('Name','Charakterystyka amplitudowa (liniowa)');
hold on;
for i = 1:length(orders)
    N = orders(i);
    [B, A] = butter(N, omega3dB, 's');  % filtr analogowy
    H = tf(B, A);
    [mag, ~] = bode(H, w);
    mag = squeeze(mag);
    plot(f, 20*log10(mag), 'DisplayName', ['N = ' num2str(N)]);
end
title('Charakterystyka amplitudowa (skala liniowa)');
xlabel('Częstotliwość [Hz]');
ylabel('20log_{10}(|H(j\omega)|)');
legend show;
grid on;

%% CHARAKTERYSTYKA AMPLITUDOWA – skala logarytmiczna
figure('Name','Charakterystyka amplitudowa (logarytmiczna)');
hold on;
for i = 1:length(orders)
    N = orders(i);
    [B, A] = butter(N, omega3dB, 's');
    H = tf(B, A);
    [mag, ~] = bode(H, w);
    mag = squeeze(mag);
    semilogx(f, 20*log10(mag), 'DisplayName', ['N = ' num2str(N)]);
end
title('Charakterystyka amplitudowa (skala logarytmiczna)');
xlabel('Częstotliwość [Hz]');
ylabel('20log_{10}(|H(j\omega)|)');
legend show;
grid on;

%% CHARAKTERYSTYKA FAZOWA
figure('Name','Charakterystyka fazowa');
hold on;
for i = 1:length(orders)
    N = orders(i);
    [B, A] = butter(N, omega3dB, 's');
    H = tf(B, A);
    [~, phase] = bode(H, w);
    phase = squeeze(phase);
    plot(f, phase, 'DisplayName', ['N = ' num2str(N)]);
end
title('Charakterystyka fazowa');
xlabel('Częstotliwość [Hz]');
ylabel('Faza [stopnie]');
legend show;
grid on;

%% ODPOWIEDŹ IMPULSOWA – tylko dla N = 4
[B4, A4] = butter(4, omega3dB, 's');
H4 = tf(B4, A4);

figure('Name','Odpowiedź impulsowa – filtr N=4');
impulse(H4);
title('Odpowiedź impulsowa filtru Butterwortha (N = 4)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

%% ODPOWIEDŹ SKOKOWA – tylko dla N = 4
figure('Name','Odpowiedź skokowa – filtr N=4');
step(H4);
title('Odpowiedź skokowa filtru Butterwortha (N = 4)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;
