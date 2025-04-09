clc;
clear;

% --- Parametry podstawowe ---
fc = 96e6; % częstotliwość środkowa stacji FM
fs = 2*(fc + 5e6); % sztucznie przyjęta częstotliwość próbkowania dla normalizacji

% ============================
% 1. FILTR TESTOWY (96 MHz ±1 MHz)
% ============================
passband = [95e6 97e6];
stopband = [94e6 98e6];

Wp = passband / (fs/2);
Ws = stopband / (fs/2);

Rp = 3; Rs = 40; % dopuszczalne zafalowanie i tłumienie

[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[b, a] = ellip(n, Rp, Rs, Wn, 'bandpass');

% Obliczenie odpowiedzi częstotliwościowej
[h, f] = freqz(b, a, 2048, fs);
magnitude = 20*log10(abs(h)); % Amplituda w dB

% Wykres charakterystyki częstotliwościowej
figure;
plot(f, magnitude, 'LineWidth', 1.5);
grid on;
title('Filtr TESTOWY 96 MHz ±1 MHz');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
ylim([-80 5]);

% Zaznaczenie punktów charakterystycznych
yline(-3, '--r', '-3 dB'); % granica pasma przepustowego
yline(-40, '--k', '-40 dB'); % granica pasma zaporowego
xline(passband(1), '--g', 'granica pasma'); % granica pasma przepustowego
xline(passband(2), '--g');
xline(stopband(1), '--m', 'granica zaporowa'); % granica pasma zaporowego
xline(stopband(2), '--m');

% ============================
% 2. FILTR DOCELOWY (96 MHz ±100 kHz)
% ============================
passband = [95.9e6 96.1e6];
stopband = [95.8e6 96.2e6];

Wp = passband / (fs/2);
Ws = stopband / (fs/2);

[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[b, a] = ellip(n, Rp, Rs, Wn, 'bandpass');

% Obliczenie odpowiedzi częstotliwościowej
[h, f] = freqz(b, a, 2048, fs);
magnitude = 20*log10(abs(h)); % Amplituda w dB

% Wykres charakterystyki częstotliwościowej
figure;
plot(f, magnitude, 'LineWidth', 1.5);
grid on;
title('Filtr DOCELOWY 96 MHz ±100 kHz');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
ylim([-80 5]);

% Zaznaczenie punktów charakterystycznych
yline(-3, '--r', '-3 dB'); % granica pasma przepustowego
yline(-40, '--k', '-40 dB'); % granica pasma zaporowego
xline(passband(1), '--g', 'granica pasma'); % granica pasma przepustowego
xline(passband(2), '--g');
xline(stopband(1), '--m', 'granica zaporowa'); % granica pasma zaporowego
xline(stopband(2), '--m');
