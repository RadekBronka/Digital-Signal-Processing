clear all;
close all;
clearvars;
clc;

A=230;
f=50;
T=0.1;

%próbkowanie
fs1 = 10000; 
fs2 = 500;   
fs3 = 200; 

%czas
t1=0:1/fs1:T;
t2=0:1/fs2:T;
t3=0:1/fs3:T;

%sygnały
y1=A*sin(2*pi*f*t1);
y2=A*sin(2*pi*f*t2);
y3=A*sin(2*pi*f*t3);

figure;
hold on;
plot(t1, y1, 'b-', 'LineWidth', 1);  
plot(t2, y2, 'ro', 'MarkerSize', 5, 'LineWidth', 1); 
plot(t3, y3, 'kx', 'MarkerSize', 5, 'LineWidth', 1); 
hold off;

title('Sinusoida 50 Hz dla różnych fs');
legend('fs = 10 kHz', 'fs = 500 Hz', 'fs = 200 Hz');
grid on;
