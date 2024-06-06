clear; clc;
fs=1000;f0=50; 

b0=1;
b1=-2*cos(2*pi*f0/fs); 
b2=1;
r=0.80;
a0=1; 
a1=-2*r*cos(2*pi*f0/fs); 
a2=(r*r);

b=[b2,b1,b0];
a=[a2,a1,a0];
[h,w]=freqz(b,a);
hf=abs(h);
figure(1);
plot(w*fs/(2*pi),hf);
title('Magnitude Response');
xlabel('Freguency in Hz');
figure(2);
zplane(b,a);legend('zero','plot');

%%

clear; clc;

% Filter specifications
fs = 1000; % Sampling frequency in Hz
f0 = 50; % Notch frequency in Hz
r = 0.80; % Pole radius
b0 = 1; % Coefficient b0

% Compute filter coefficients
b1 = -2 * cos(2 * pi * f0 / fs);
b2 = 1;
a0 = 1;
a1 = -2 * r * cos(2 * pi * f0 / fs);
a2 = r * r;

% Define filter coefficients arrays
b = [b0, b1, b2];
a = [a0, a1, a2];   

% Frequency response
[h, w] = freqz(b, a);

% Plot magnitude response
hf = abs(h);
figure(1);
plot(w * fs / (2 * pi), hf);
title('Magnitude Response');
xlabel('Frequency in Hz');
ylabel('Magnitude');
grid on;

% Plot zero-pole plot
figure(2);
zplane(b, a);
legend('zero', 'pole');
title('Zero-Pole Plot');

% Plot impulse response
figure(3);
impulse_response = impz(b, a);
stem(impulse_response);
title('Impulse Response');
xlabel('Sample');
ylabel('Amplitude');
grid on;