
% Script - Fourier Transform Lab 2
% functions 
%   Generate square waveform 
%   Fourier expand n cycles 
% 
% case 1:
% square waveform
%case 2:
%saw tooth waveform

Waveform = 2;

switch Waveform 

    case 1;

        %Define the parameters- version 2
        frequency = 1000; % Frequency of the square wave (in Hz)
        duration = 0.002; % Duration of the square wave (in seconds)
        sampling_rate = 100000; % Sampling rate (samples per second)
        A = 1;
        w = 2 * pi * frequency;
        
        % Generate the time vector
        t = -0.002:1/sampling_rate:duration;
        square_wave_1 = 0;
        % Generate the even square wave
        square_wave = A * square(2 * pi * frequency * t + 2 * pi * 0.25);
        plot(t,square_wave);
        hold on;
        
        for n = 1:1:10
            square_wave_1 = square_wave_1 + 2/(n*pi)*sin(n*pi/2)*cos(n*w*t) - (2/(n*pi))*sin(3*pi*n/2)*cos(n*w*t);
            h1 = plot(t,square_wave_1);
            pause(1);
            hold on;
        end
        hold off;
        title('Fourier SquareWave')
        grid on;

    case 2;

        %Define the parameters- version 1
        frequency = 1000; % Frequency of the sawtooth wave (in Hz)
        duration = 0.002; % Duration of the sawtooth wave (in seconds)
        sampling_rate = 100000; % Sampling rate (samples per second)
        A = 1;
        w = 2 * pi * frequency;
        
        % Generate the time vector
        t = -0.002:1/sampling_rate:duration;
        Sawtooth_wave_1 = 0;
        % Generate the even sawtooth wave
        Sawtooth_wave = A * sawtooth(2*pi*frequency*t+pi);
        plot(t,Sawtooth_wave)
        
        y = fft(Sawtooth_wave);
        plot(t,Sawtooth_wave);
        hold on;
        %square_wave_1 = square_wave_1 + 2/pi*sin(1/2*n*pi)/n*cos(n*w*t);
         for n = 1:1:10
             Sawtooth_wave_1 = Sawtooth_wave_1 - 2/(n*pi)*(-1)^n*sin(n*w*t);
             fprintf( 'Frequecy = %d, Amplitude = %d \n',n*frequency,-2/(n*pi)*(-1)^n);
             h1 = plot(t,Sawtooth_wave_1);
             hold on;
        % fprintf('Frequency= %d, Amplitude = %d \n',n*w/(2*pi),2/pi*sin(1/2*n*pi)/n);
            pause(0.5);
            hold on;
        end
        hold off;
        title('Fourier Sawtooth Waveform')
        grid on;
    otherwise
        %Define the parameters- version 2
        frequency = 1000; % Frequency of the square wave (in Hz)
        duration = 0.002; % Duration of the square wave (in seconds)
        sampling_rate = 100000; % Sampling rate (samples per second)
        A = 1;
        Offset = 0.5;
        w = 2 * pi * frequency;

        % Generate the time vector
        t = -0.002:1/sampling_rate:duration;
        square_wave_1 = 0;
        % Generate the even square wave
        square_wave = A * square(2 * pi * frequency * t + 2 * pi * 0.25);
        plot(t,square_wave);
        hold on;
        Y = fft(square_wave);
        f = frequency*(0:(65-1)/2)/65;
        P2 = abs(Y/65);
        P1 = P2(1:(65+1)/2);
        P1(2:end) = 2*P1(2:end);
end

% %Define the parameters- version 2
% frequency = 1000; % Frequency of the square wave (in Hz)
% duration = 0.002; % Duration of the square wave (in seconds)
% sampling_rate = 100000; % Sampling rate (samples per second)
% A = 1;
% Offset = 0.5;
% w = 2 * pi * frequency;
% 
% % Generate the time vector
% t = -0.002:1/sampling_rate:duration;
% square_wave_1 = 0;
% % Generate the even square wave
% square_wave = A * square(2 * pi * frequency * t + 2 * pi * 0.25);
% plot(t,square_wave);
% hold on;
% 
% for n = 1:1:10
%     square_wave_1 = square_wave_1 + 2/(n*pi)*sin(n*pi/2)*cos(n*w*t) - (2/(n*pi))*sin(3*pi*n/2)*cos(n*w*t);
%     h1 = plot(t,square_wave_1);
%     Amplitude = 2/(n*pi)*sin(n*pi/2)-(2/(n*pi))*sin(3*pi*n/2);
%     %fprintf('Frequency= %d, Amplitude = %d \n',n*frequency,Amplitude);
%     pause(1);
%     hold on;
% end
% hold off;
% title('Fourier SquareWave')
% grid on;


% f = 1/(200e-9);%方波信号的频率，5Mhz，200ns
% L = 512;% 每个周期采样的点数
% Fs = L*f;% Sampling frequency ，采样率=信号频率*每个周期的采样数                   
% T = 1/Fs;% Sampling period     两次采样的时间间隔
% t = (0:L-1)*T;% Time vector    通过采样时间间隔和采样的点数生成时间矩阵。这里选取一个周期的方波信号：512个点=1个方波周期，时间从0开始所以最后减去一个点。
% x = square(t*f*2*pi); %square可以方便的生成方波信号，生成方式和正弦信号一样，角速度*时间，官方说明文档有很多实例。
% plot(t,x);%画出方波信号幅值随时间的变化图 

%Calculate the Fourier series coefficients of a periodic square wave with an amplitude of 1 V, offset DC voltage of 0 V and period of 1 ms, assuming the signal is an even function of time.
% t = linspace(0,3);
% x = square(t);
% plot(t/pi,x,'.-',t/pi,sin(t))
% xlabel('t / \pi')
% grid on

% freq=1/1e-3;
% offset=0.5;
% amp=0.5;
% duty=50;
% t=0:0.01:1;%100 seconds
% %t=2*pi*freq
% %t = linspace(0,3);
% sq_wav=offset+amp*square(2*pi*freq*t,duty);
% plot(sq_wav)

% A=0.5;
% offset = 0.5;
% fo=10*pi;
% w=pi/4;
% t=-2:0.01:2;   
% sq=offset+A*square(fo*t+w);
% 
% plot (t, sq)
% axis ([0  1  -2  2])

% % Define the parameters
% frequency = 1000; % Frequency of the square wave (in Hz)
% duration = 1; % Duration of the square wave (in seconds)
% sampling_rate = 1000; % Sampling rate (samples per second)
% 
% % Generate the time vector
% t = 0:1/sampling_rate:duration;
% 
% % Generate the even square wave
% square_wave = square(2 * pi * frequency * t);
% 
% % Plot the square wave
% plot(t, square_wave);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Even Square Wave');
% grid on;

% %%% 正弦波  对应的方波和三角波
% clc
% clear
% close
% 
% %% 参数
% fz = 600000;    %%采样频率
% count = 30;     %%三角波个数及方波个数

% %% 主程序
% xL = [1/fz:(1/fz):1];
% accuracy = fz/count;
% x = 1:1:fz;
% 
% y_sin = sin(2*pi*x/fz);
% y_sin = (0.8*y_sin)/2+0.5; %使正弦波幅值在0.1~0.9之间
% flag = 1;
% y_san(1)=0;
% 
% for i=1:fz
%     y_san(i+1) = y_san(i)+flag * (1/accuracy); 
%     if rem(i,accuracy) == 0
%        flag = -flag;      
%     end
% end
% y_san = y_san(1:fz);
% for i = 1:fz
%    if y_san(i)>y_sin(i)
%       y_fang(i) = 0;
%    else
%       y_fang(i) = 1;
%    end
% end
% figure(1)
% plot(xL,y_san,'k')
% hold on
% plot(xL,y_sin,'k')
% axis([0 1 -0.1 1.1])
% %axis off
% figure(2)
% plot(xL,y_fang,'k')
% hold on
% plot(xL,y_sin,'k')
% axis([0 1 -0.1 1.1])
% %axis off

% % 第一个图
% clear
% clc
% 
% T1 = 1; T2 = 2; T3 = 4/5; T4 = 4/7;
% T0 = 2*pi; % 周期 如果是相加项 取所有周期分子的公约数 再乘pi
% omiga = 2*pi/T0;
% 
% fsin=@(n,t) sin(n*omiga*t).*square(t);
% 
% fcos=@(n,t) cos(n*omiga*t).*square(t);
% 
% N=20;
% 
% Fsin=zeros(1,N+1);
% 
% Fcos=zeros(1,N+1);
% 
% for n=1:N
% 
%     Fsin(n)=quad(@(t)fsin(n,t),-T0/2,T0/2,1e-9)/T0*2;
% 
%     Fcos(n)=quad(@(t)fcos(n,t),-T0/2,T0/2,1e-9)/T0*2;
% 
% end
% 
% n = 0;
% a0 = quad(@(t)fcos(n,t),-T0/2,T0/2,1e-9)/T0;
% 
% % 把极小值弄为0
% Fcos(Fcos<0.01)=0;
% Fsin(Fsin<0.01)=0;
% 
% subplot(211),stem((1:N+1).*omiga,Fcos);title('an');
% grid on
% ylim([-1,1]) % 设置幅值显示范围
% xlabel('n')
% ylabel('幅值')
% 
% subplot(212),stem((1:N+1).*omiga,Fsin);title('bn');
% grid on
% ylim([-1,1]) % 设置幅值显示范围
% xlabel('n')
% ylabel('幅值')
% 
% figure
% subplot(211);
% stem((1:N+1).*omiga,sqrt(Fcos.^2+Fsin.^2));
% grid on
% title('幅值谱','fontsize',13);
% xlabel('\omega','fontsize',12);
% ylabel('A_n','fontsize',12)
% 
% subplot(212);
% phi = atan(Fcos./Fsin);
% phi(isnan(phi))=0; % bn为0时 算出来NaN无效值 这一步去除无效值 
% stem((1:N+1).*omiga,phi);
% grid on
% ylim([-pi pi])
% xlim([0 25])
% title('相位谱','fontsize',13)
% xlabel('\omega','fontsize',12);
% ylabel('\phi_n','fontsize',12)
% 
% 
% 
% % Define the parameters
% amplitude = 0.5;         % Amplitude (1 V)
% offset = 0.5;            % DC Offset (0 V)
% period = 0.001;        % Period (1 ms)
% frequency = 1 / period; % Frequency (1 kHz)
% duration = 0.02;       % Duration of the signal (20 ms)
% sampling_rate = 10e6;  % Sampling rate (10 MSamples per second)
% 
% % Generate the time vector
% t = 0:1/sampling_rate:duration;
% 
% % Generate the square wave signal
% square_wave = amplitude * square(2 * pi * frequency * t) + offset;
% 
% % Plot the square wave
% plot(t, square_wave);
% xlabel('Time (s)');
% ylabel('Amplitude (V)');
% title('Periodic Square Wave');
% grid on;


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
