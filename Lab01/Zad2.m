clear all;
close all;
clearvars;
clc;

% Parametry sygnału z zadania 1A
A = 230;
f = 50;
T = 0.1;

% Częstotliwości próbkowania
fs1 = 10000;  
fs3 = 200;    


t1 = 0:1/fs1:T;  
t3 = 0:1/fs3:T;   %punkty probek dla fs=200Hz
ts = 0:1/fs1:T;   %zmienilem z 0:0.02:T, poniewaz wtedy wychodzily kanciaste wykresy zamiast sinusa, teraz jest sinus

% Sygnał pseudo analog
y_analog = A * sin(2 * pi * f * ts);
% Sygnał próbkowany 200Hz
y3 = A * sin(2 * pi * f * t3);

%Kod z pdf
xhat = zeros(size(ts));

for idx = 1:length(ts)
    t = ts(idx);
    xhat(idx) = 0; 
    for n = -3:3  % Zakres sinc
        sample_idx = round((t * fs3)) + n; 
        if sample_idx >= 1 && sample_idx <= length(y3) 
            %T=1/fs3, nT=t3(sample_idx)
            xhat(idx) = xhat(idx) + y3(sample_idx) * sinc(fs3 * (t - t3(sample_idx)));
        end
    end
end




figure;
plot(ts, y_analog, 'r', 'LineWidth', 1.5); hold on; 
plot(ts, xhat, 'b--', 'LineWidth', 1);  
legend('Sygnał pseudo-analogowy', 'Zrekonstruowany');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

% Błąd rekonstrukcji
error_signal = y_analog - xhat;

figure;
plot(ts, error_signal, 'k', 'LineWidth', 1);
xlabel('Czas [s]');
ylabel('Błąd rekonstrukcji');
title('Błąd rekonstrukcji sygnału (ograniczona suma sinc)');
grid on;
