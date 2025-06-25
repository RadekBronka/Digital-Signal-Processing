clc; clear; close all;

%% Parametry
orders = 2:8;            % Rzędy
f0_values = [1, 10, 100]; % Częstotliwości 
Rp = 3;                  % Ripple [dB]


for f0 = f0_values
    w0 = 2*pi*f0;           % Pulsacja graniczna
    figure;
    for N = orders
        % Projekt filtrów
        [zC, pC, kC] = cheb1ap(N, Rp);
        [bC, aC] = zp2tf(zC, pC, kC);
        [zB, pB, kB] = butter(N, 1, 's');
        [bB, aB] = zp2tf(zB, pB, kB);
        % Skala częstotliwości
        [bC, aC] = lp2lp(bC, aC, w0);
        [bB, aB] = lp2lp(bB, aB, w0);
        
        f = linspace(f0/100, f0*100, 1000);
        w = 2*pi*f;                           % Przeliczenie na pulsację
        Hc = freqs(bC, aC, w);
        Hb = freqs(bB, aB, w);
        % Wykres
        subplot(2,4,N);
        semilogx(w/(2*pi), 20*log10(abs(Hc)), 'b', 'LineWidth', 1.5);
        hold on;
        semilogx(w/(2*pi), 20*log10(abs(Hb)), 'r--', 'LineWidth', 1.5);
        grid on;
        title(['N = ' num2str(N) ', f_0 = ' num2str(f0) ' Hz']);
        xlabel('f [Hz]');
        ylabel('|H| [dB]');
        ylim([-80 5]);
        legend('Czebyszew I','Butterworth');
    end
    sgtitle(['Porównanie filtrów analogowych N=2...8, f_0 = ' num2str(f0) ' Hz']);
end
