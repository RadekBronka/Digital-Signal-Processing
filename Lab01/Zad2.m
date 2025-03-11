clc; clear; close all;

%% Parametry sygnału
A = 230;      
f = 50;      
fs_analog = 10000; 
fs_sampled = 200; 
T = 0.1;      % Czas trwania [s]


t_analog = 0:1/fs_analog:T;
x_analog = A * sin(2 * pi * f * t_analog);


t_sampled = 0:1/fs_sampled:T;
x_sampled = A * sin(2 * pi * f * t_sampled);

%rekonstrukcja
t_reconstruct = t_analog;  % mialo byc w chw
x_reconstruct = zeros(size(t_reconstruct));

% Rekonstrukcja metodą sinc()
for idx = 1:length(t_reconstruct)
    t = t_reconstruct(idx);
    x_reconstruct(idx) = 0;
    for n = 1:length(t_sampled)  % Sumujemy po próbkach
        T_s = 1/fs_sampled;  % Okres próbkowania
        x_reconstruct(idx) = x_reconstruct(idx) + x_sampled(n) * sinc((t - t_sampled(n)) / T_s);
    end
end

%% Wykresy
figure;
hold on;
plot(t_analog, x_analog, 'b', 'DisplayName', 'Pseudo-analogowy (fs = 10 kHz)');
stem(t_sampled, x_sampled, 'r', 'DisplayName', 'Spróbkowany (fs = 200 Hz)', 'MarkerSize', 5);
plot(t_reconstruct, x_reconstruct, 'g', 'DisplayName', 'Zrekonstruowany (sinc)', 'LineWidth', 1.2);
xlabel('Czas [s]');
ylabel('Amplituda [V]');
title('Rekonstrukcja sygnału za pomocą sinc()');
legend;
grid on;
hold off;

%% Obliczenie błędu rekonstrukcji
error_signal = x_analog - x_reconstruct;

figure;
plot(t_analog, error_signal, 'k');
xlabel('Czas [s]');
ylabel('Błąd rekonstrukcji [V]');
title('Błąd rekonstrukcji');
grid on;
