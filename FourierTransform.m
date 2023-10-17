% Script - Fourier Transform Lab 2
% functions 
%   Generate square waveform 
%   Fourier expand n cycles 
% 
%         

%Define the parameters- version 1
frequency = 1000; % Frequency of the square wave (in Hz)
duration = 0.002; % Duration of the square wave (in seconds)
sampling_rate = 100000; % Sampling rate (samples per second)
A = 0.5;
Offset = 0.5;
w = 2 * pi * frequency;

% Generate the time vector
t = -0.002:1/sampling_rate:duration;
square_wave_1 = 1/2;
% Generate the even square wave
square_wave = Offset + A * square(2 * pi * frequency * t + 2 * pi * 0.25);
plot(t,square_wave);
hold on;

for n = 1:2:10
    square_wave_1 = square_wave_1 + 2/pi*sin(1/2*n*pi)/n*cos(n*w*t);
    h1 = plot(t,square_wave_1);
    pause(0.5);
    hold on;
end
hold off;
title('square FourierSquareWave')
grid on;

% square_wave_1 = 1/2+2*1/pi*sin(pi/2)*cos(2*pi*frequency*t);
% plot(t,square_wave_1);
% hold on;
% square_wave_2 = 1/2+2*1/pi*sin(pi)/2*cos(2*w*t);
% plot(t,square_wave_2);
% hold on 
% square_wave_3 = square_wave_1+2*1/pi*sin(pi*3/2)/3*cos(3*w*t);
% plot(t,square_wave_3);
% hold off


% %Define the parameters- version 1
% frequency = 1000; % Frequency of the square wave (in Hz)
% duration = 0.002; % Duration of the square wave (in seconds)
% sampling_rate = 100000; % Sampling rate (samples per second)
% A = 0.5;
% Offset = 0.5;
% w = 2 * pi * frequency;
% 
% % Generate the time vector
% t = -0.002:1/sampling_rate:duration;
% square_wave_1 = 1/2;
% % Generate the even square wave
% square_wave = Offset + A * square(2 * pi * frequency * t + 2 * pi * 0.25);
% plot(t,square_wave);
% hold on;
% square_wave_1 = 1/2+2*1/pi*sin(pi/2)*cos(2*pi*frequency*t);
% plot(t,square_wave_1);
% hold on;
% square_wave_2 = 1/2+2*1/pi*sin(pi)/2*cos(2*w*t);
% plot(t,square_wave_2);
% hold on 
% square_wave_3 = square_wave_1+2*1/pi*sin(pi*3/2)/3*cos(3*w*t);
% plot(t,square_wave_3);
% hold off

% frequency = 1000;
% t = linspace(-0.04,0.04,1000);
% x = square(2*pi*frequency*t);
% %f(t) = A*sin(wt+fai) w= 2pif
% plot(t,x)
% plot(t/pi,x,'.-',t/pi,1.15*sin(2*t))
% xlabel('t / \pi')
% grid on

% t1 = linspace(0,0.1,100);
% xt1 = 0;
% A =5;
% f0 = 10;
% for n = 1:2:10
%     xt1 = xt1+4*A/pi*sin(n*2*pi*f0*t1)/n;
%     h1 = plot(t1,xt1);
%     pause(0.5);
%     hold on;
% end
% title('square x(t)')
% grid on;
% XTickLable={'0','T_{0}/2','T_{0}'};
% set(gca,'Xtick',[0:0.05:0.1],'XTickLabel',XTickLable,'TickDir','out');

% x = linspace(-6.18,6.18);
% y1 = sin(x);
% plot(x,y1)
% 
% hold on
% y2 = sin(2*pi*1000*x);
% plot(x,y2)
% hold off



%Plot the square wave
% plot(t, square_wave);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Even Square Wave');
% grid on;
