load('ECG100.mat');      % Załaduj plik
ekg = val(1, :);         % Wybierz pierwszy kanał sygnału
fs = 1000;               % Częstotliwość próbkowania
t = (0:length(ekg)-1)/fs;

% Dodaj zakłócenie 50 Hz
f_noise = 50;
noise = 10* sin(2*pi*f_noise*t);
ekg_noisy = ekg + noise;

% Referencja do filtra LMSa
ref = sin(2*pi*f_noise*t)';

% Parametry filtra LMS
mu = 0.01;
filterOrder = 12;

lms = dsp.LMSFilter('Length', filterOrder, 'StepSize', mu);
[~, err] = lms(ref, ekg_noisy');
ekg_filtered = err';

% Wizualizacja
figure;
subplot(3,1,1);
plot(t, ekg);
title('Oryginalny sygnał EKG');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(3,1,2);
plot(t, ekg_noisy);
title('EKG z zakłóceniem 50 Hz');
xlabel('Czas [s]');
ylabel('Amplituda');

subplot(3,1,3);
plot(t, ekg_filtered);
title('EKG po filtrze LMS');
xlabel('Czas [s]');
ylabel('Amplituda');

mae = mean(abs(ekg - ekg_filtered));
fprintf('MAE = %.6f\n', mae);
