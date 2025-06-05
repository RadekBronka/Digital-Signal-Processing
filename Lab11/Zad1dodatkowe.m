[x, fs] = audioread('DontWorryBeHappy.wav');
x = x(1:min(end, 4*fs), :);
x = double(x);

[Y_16, MSE_16] = g726_adpcm(x, fs, 16000);
[Y_32, MSE_32] = g726_adpcm(x, fs, 32000);

% Porównanie
figure;
subplot(2,1,1);
plot(x(:,1)); hold on; plot(Y_16(:,1)); plot(Y_32(:,1));
legend('Oryg.', 'Rek. 16 kbps', 'Rek. 32 kbps');
title(['Lewy kanał: MSE (16 kbps) = ', num2str(MSE_16(1), '%.4g'), ...
       ', MSE (32 kbps) = ', num2str(MSE_32(1), '%.4g')]);

subplot(2,1,2);
plot(x(:,2)); hold on; plot(Y_16(:,2)); plot(Y_32(:,2));
legend('Oryg.', 'Rek. 16 kbps', 'Rek. 32 kbps');
title(['Prawy kanał: MSE (16 kbps) = ', num2str(MSE_16(2), '%.4g'), ...
       ', MSE (32 kbps) = ', num2str(MSE_32(2), '%.4g')]);

% Odsłuch
fprintf('Odtwarzam zrekonstruowany sygnał 32 kbps...\n');
sound(Y_32, fs);
