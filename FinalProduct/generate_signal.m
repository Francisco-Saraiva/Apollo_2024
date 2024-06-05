%% This code creates a memory map of the binary file
close all; clear all
%% INPUT PARAMETERS (KNOWN BEFOREHAND)

samp_freq = 30000; %sampling frequency (in Hz)

nchan = 64; %this has to be known in advance (if wrong, data will be loaded incorrectly)

chans_to_open =[64 62 59 57 55 53 51 47 50 52 54 56 58 60 61 63 37 39 41 46 44 33 35 48 36 34 49 43 45 42 40 38 27 25 23 20 22 31 29 18 30 32 15 21 19 24 26 28 2 4 5 7 9 11 13 17 16 14 12 10 8 6 3 1]; %NEW 4chan

%% LOAD FILES

[FileName,PathName,FilterIndex] = uigetfile('.dat'); %select file

contFile=fullfile(PathName,FileName); %create full file path with file name

s=dir(contFile);
file_size=s.bytes; % determine file size in byte

samples=file_size/2/nchan;   % calculate the number of samples

m=memmapfile(contFile,'Format',{'int16' [nchan samples] 'mapped'}); %create memory map of the file

data = m.Data;

% To open a specific electrode site(x) you need to convert the memory map to double data:
% chan(x,:) = double(data.mapped(chans_to_open(x)))
%% MEASURES CALCULATION

% Get desired channel from user
channel_num = input("Please input the number of the channel you want to open:");

channel_data = double(data.mapped(channel_num,:));

total_time = samples / samp_freq; % in seconds
total_time_1 = samples /(60 * samp_freq); % in minutes

% Other stuff we might want...

% uncomment the following lines if you wish to visualize the signal
%{
min_delta = 1/samp_freq;
time_vec = 0:min_delta:total_time-min_delta;

% plot(time_vec, channel_1)  
%}
%% SIGNAL GENERATION: 16BIT

% ~~~~~~~~~~~~~~~~~~~~~~~~ PORT SETTINGS ~~~~~~~~~~~~~~~~~~~~~~~~~
port_name = 'COM12';
baudrate = 115200;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

port = serialport(port_name, baudrate);
channel_str = string(channel_data);

% send the data through port
for i=1:length(channel_str)
    writeline(port,channel_str(i));
    pause(0.001);
end

disp("Generation Complete");
clear("port")