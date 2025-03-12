fs=100;
T=1;
A=230;
t=0:1/fs:T-1/fs;

for f = 0:5:300
    y = A * sin(2 * pi * f * t);  
    plot(t, y, 'b');             
    title(sprintf('Iteracja: %d, Częstotliwość: %d Hz', f/5 + 1, f));
    grid on;
    axis([0 T -A A]);             
    pause(0.2);                   
end

freq_groups = {[5, 105, 205], [95, 195, 295], [95, 105]};
% dla dwoch pierwszych grup wygladaja tak samo, chyba przez aliasing?
figure;
for i = 1:3
    subplot(3,1,i); hold on; grid on;
    for f = freq_groups{i}
        plot(t, A * sin(2 * pi * f * t));
    end
    title(['Porównanie sin ', num2str(freq_groups{i})]);
    legend show;
end

figure;
for i = 1:3
    subplot(3,1,i); hold on; grid on;
    for f = freq_groups{i}
        plot(t, A * cos(2 * pi * f * t));
    end
    title(['Porównanie cos ', num2str(freq_groups{i})]);
    legend show;
end
