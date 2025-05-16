clear; close all; clc;

%% Parametry
fs = 8000;
M = 256;              % długość filtru (tyle co h_real)
mi = 0.01;            % krok LMS
N = 12000;            % długość sygnału

% Rzeczywista odpowiedź impulsowa
h_real = zeros(M, 1);
h_real(256) = 0.8;
h_real(121) = -0.5;
h_real(31)  = 0.1;

%% --- CZĘŚĆ 1: IDENTYFIKACJA Z SYGNAŁU MOWY ---
[x_mowa, fs] = audioread('mowa8000.wav');
x_mowa = x_mowa(1:min(N, end));
x_mowa = x_mowa(:);  % kolumna
d_mowa = filter(h_real, 1, x_mowa) + 0.001*randn(size(x_mowa));

% Inicjalizacja
h = zeros(M,1);
bx = zeros(M,1);
y = zeros(size(x_mowa));
e = zeros(size(x_mowa));

% LMS
for n = 1:length(x_mowa)
    bx = [x_mowa(n); bx(1:M-1)];
    y(n) = h' * bx;
    e(n) = d_mowa(n) - y(n);
    h = h + mi * e(n) * bx;
end

% Wykres – identyfikacja z sygnału mowy
figure;
subplot(2,1,1);
stem(h_real, 'b', 'DisplayName', 'h_{real}');
hold on;
stem(h, 'r--', 'DisplayName', 'h_{estymowana}');
legend;
xlabel('Próbka');
ylabel('Amplituda');
title('Identyfikacja z sygnału mowy');
grid on;

hold off;  % <- DODAJ TO PRZED KOLEJNYM SUBPLOTEM
% MSE
mse_mowa = mean((h - h_real).^2);
disp(['MSE (mowa): ', num2str(mse_mowa)]);


%% --- CZĘŚĆ 2: IDENTYFIKACJA Z BIAŁEGO SZUMU ---
mi = 0.002;
x_szum = randn(N, 1);  % sekwencja treningowa
d_szum = filter(h_real, 1, x_szum) + 0.001*randn(size(x_szum));

% Inicjalizacja
h2 = zeros(M,1);
bx = zeros(M,1);
y2 = zeros(size(x_szum));
e2 = zeros(size(x_szum));

% LMS
for n = 1:length(x_szum)
    bx = [x_szum(n); bx(1:M-1)];
    y2(n) = h2' * bx;
    e2(n) = d_szum(n) - y2(n);
    h2 = h2 + mi * e2(n) * bx;
end

% Wykres – identyfikacja z szumu białego
subplot(2,1,2);
stem(h_real, 'b', 'DisplayName', 'h_{real}','MarkerSize', 5);
hold on;
stem(h2, 'g--', 'DisplayName', 'h_{estymowana\_szum}');
legend;
xlabel('Próbka');
ylabel('Amplituda');
title('Identyfikacja z białego szumu');
grid on;

% MSE
mse_szum = mean((h2 - h_real).^2);
disp(['MSE (szum): ', num2str(mse_szum)]);
