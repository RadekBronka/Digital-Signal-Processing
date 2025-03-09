clear all;
close all;
clearvars;
clc;

A=230;
f=50;
T=1;

%próbkowanie
fs1 = 10000; 
fs2 = 51;   
fs3 = 50;
fs4 = 49;
 
%fs2 = 26;   
%fs3 = 25;
%fs4 = 24;

%czas
t1=0:1/fs1:T;
t2=0:1/fs2:T;
t3=0:1/fs3:T;
t4=0:1/fs4:T;
%sygnały
y1=A*sin(2*pi*f*t1);
y2=A*sin(2*pi*f*t2);
y3=A*sin(2*pi*f*t3);
y4=A*sin(2*pi*f*t4);


figure;
hold on;
plot(t1, y1, 'b-', 'LineWidth', 1);  
plot(t2, y2, 'g-o', 'MarkerSize', 5, 'LineWidth', 1); 
plot(t3, y3, 'r-o', 'MarkerSize', 5, 'LineWidth', 1); 
plot(t4, y4, 'k-o', 'MarkerSize', 5, 'LineWidth', 1); 
hold off;

title('Sinusoida 50 Hz dla różnych fs');
legend('fs = 10 kHz', 'fs = 51 Hz', 'fs = 50 Hz','fs = 49 Hz');
%legend('fs = 10 kHz', 'fs = 26 Hz', 'fs = 25 Hz','fs = 24 Hz');
grid on;
