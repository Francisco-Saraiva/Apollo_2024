clear; clc;
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~ PORT CONFIGURATIONS ~~~~~~~~~~~~~~~~~~~~~~~~

port_name = 'COM12';
baudrate = 115200;
terminator = "LF";

% Define serial port settings
port = serialport(port_name, baudrate);
configureTerminator(port, terminator);

%% ~~~~~~~~~~~~~~~~~~~~~ SIGNAL GENERATION FOR TESTING ~~~~~~~~~~~~~~~~~~

% Define parameters for time domain
Fs = 1000; % Sampling frequency (Hz)
time_end = 1; % end time of the signal (in seconds)
t = 0:1/Fs:(time_end-1/Fs); % Time vector from 0 to 1 second with sampling interval 1/Fs

% PURE SINUSOIDAL SIGNAL
f1 = 60; % Frequency of the sinusoid (Hz)
A1 = 1000; % Amplitude of the sinusoid

sinusoidal_signal = round(A1 * sin(2*pi*f1*t));

% RECTANGULAR SIGNAL
A2 = 1000;  % Amplitude of the signal
LowLimit = 0.2;  % When the impulse starts
HighLimit = 0.6; % When the impulse ends

rectangular_signal = round(A2 * rectangularPulse(LowLimit, HighLimit, t));

% Simply so that the user can know what duty cycle he is using
DutyCycle = HighLimit - LowLimit; % Duty cycle of the rectangular signal

% MORE COMPLEX SIGNAL
f2 = 30; % Frequency of the 2nd sinusoidal (Hz)
A3 = 500; % Amplitude of the 1st sinusoidal
A4 = 1000; % Amplitude of the 2nd sinusoidal
complex_signal = round(A3 * sin(2*pi*f1*t) + A4 * sin(2*pi*f2*t));

% Plotting the signals (visualization)
figure;
subplot(3,1,1);
plot(t, sinusoidal_signal);
title('Pure Sinusoidal Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(t, rectangular_signal);
title('Rectangular Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,3);
plot(t, complex_signal);
title('Complex Signal');
xlabel('Time (s)');
ylabel('Amplitude');


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~ SIGNAL OUTPUT CHOICE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

user_choice = input("Choose the signal: 1-sinusoidal; 2-rectangular; 3-double frequency sinusoidal: ");

if user_choice == 2
    output_signal = rectangular_signal;
elseif user_choice == 3
    output_signal = complex_signal;
else
    output_signal = sinusoidal_signal;
end

output_signal = string(output_signal);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~ SIGNAL GENERATION ~~~~~~~~~~~~~~~~~~~~~~~~

% Generation Time (in seconds)
max_time = 60;

% DON'T CHANGE
timeout = 1/Fs;
counter = 0;


has_read = false;
times = 0;
tic;
while toc < max_time
    for i=1:length(output_signal)
        writeline(port,output_signal(i));
        if (mod(i-12,50) == 0 && times >= 11)
            clear("port")
            figure;
            w = waitforbuttonpress;
            close;
            
            % Once done, re-open the port
            port = serialport(port_name, baudrate);
            configureTerminator(port, terminator);
            %readline(port)
        else
            times = times + 1;
        end
    end
    disp("generated 1 wave");
    counter = counter + 1;
    %readline(port)
end

clear("port")
disp("Generation Complete!")

