
% Captured for approximately 100uSec
% Phased array probe
% Center frequency 3.5 MHz
% Pitch 1.2 mm
% No of elements 16
% 
% ADC Data
% ADC Sampling rate 80 MHz
% Depth Start 103uSec
%
% 16 channels (left to right)
% Channel number orientation
% Ch0 Ch1 ... Ch14 Ch15
% 
% Frame data (left to right)
% Frame ... Transmit Beam Angle (degrees)
%   0     1     2     3     4     5       6     7
% 26.25 18.75 11.25 3.75 -3.75 -11.25 -18.75 -26.25
%
% Data Packet format
% Int16 RawData[8][16][1024*8]
clc; clear;
fileID = fopen('rawData.bin','r');
data = fread(fileID,'int16');
data = double(data);
fclose(fileID);

data = data./max(abs(data(:)));

dataAcq = permute(reshape(data,[8192 16 8]),[3 2 1]);
% dataAcq = permute(reshape(data,[8192 8 16]),[2 3 1]);
% dataAcq = reshape(data,[8 16 8192]);
