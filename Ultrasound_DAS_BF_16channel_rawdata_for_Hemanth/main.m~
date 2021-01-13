
clc; clear; close all

readData()
clear data

for i = 1:size(dataAcq,1)    
    frames(i).RF_data = squeeze(dataAcq(i,:,:));
end

%%
% Probe details
c_sound = 1540;%1540;
pitch = 2e-3;
num_elements = 16;
f0 = 3.5e6;
fs = 80e6;
image_floor = -60;
image_ceiling = 0;
elem_cntrs=pitch*((1:num_elements)-((num_elements/2)+0.5));

% Reconstruction starts here
num_beams = length(frames);
[num_elements,num_samples] = size(frames(1).RF_data);
sampleFreq = fs;
C = c_sound;
thetas = -26.25:7.5:26.25;
Z_samples = C/sampleFreq/2;
beam_X = zeros(num_beams, num_samples*2);
beam_Z = zeros(num_beams, num_samples*2);

%% Module 1
%%{
t0 = 90e-6;
time_axis = (0:num_samples-1)/fs + t0;

ROI.X_begin = -0.075; %-0.015;
ROI.X_end  = 0.075; %0.015;
ROI.Z_begin = 0.000; %0.005;
ROI.Z_end  = time_axis(end)*c_sound/2; %0.065;

dx = 6e-5;
x_axis = ROI.X_begin:dx:ROI.X_end;

% z_axis = time_axis*c_sound/2; %/2
dz = 3e-5;
z_axis = ROI.Z_begin:dz:ROI.Z_end;

maximum_depth = z_axis(end)*1e3;

[x_matrix,z_matrix] = meshgrid(x_axis,z_axis);
x = x_matrix(:);
z = z_matrix(:);

Nx = length(x_axis);
Nz = length(z_axis);

%%
% Receive apodization
% Dynamically expanding receive aperture
rx_f_number = 0.75; %1.75;
rx_aperture = z/rx_f_number;
rx_aperture_distance = abs(x*ones(1,num_elements)-ones(Nx*Nz,1)*elem_cntrs);
receive_apodization = apodization(rx_aperture_distance,rx_aperture*ones(1,num_elements),'hamming');

%%
beamformed_data = zeros(Nx*Nz,1);
for BEAM = 2:num_beams-1 %1:num_beams
    tx_delay = z*cosd(thetas(BEAM)) + x*sind(thetas(BEAM));
    for E = 1:num_elements-1
        rx_delay = sqrt( (x-elem_cntrs(E)).^2 + z.^2 );
        delay = (tx_delay+rx_delay)/c_sound;
        beamformed_data(:) = beamformed_data(:) + ...
            receive_apodization(:,E).*interp1(time_axis,frames(BEAM).RF_data(E,:),delay,'spline',0);
        
        [BEAM E]
    end
end
beamformed_data(isnan(beamformed_data))=0;

reshaped_beamformed_data = reshape(beamformed_data,[numel(z_axis) numel(x_axis)]);

% compute envelope
envelope_beamformed_data = envelope(reshaped_beamformed_data);

% interpolate the requested grid
x_axisI = x_axis;
dz = 3e-5;
% z_axisI = ROI.Z_begin:dz:ROI.Z_end;
z_axisI = z_axis;
[x_matrixI,z_matrixI] = meshgrid(x_axisI,z_axisI);
% x = x_matrixI(:);
% z = z_matrixI(:);

resampled_envelope_beamformed_data = interp1(z_axis,envelope_beamformed_data,z_axisI,'linear',0);
% resampled_envelope_beamformed_data = interp2(x_matrix,z_matrix,envelope_beamformed_data,x_matrixI,z_matrixI,'linear',0);

mask1 = z>cosd(thetas(end))*sqrt(x.^2+z.^2);
mask2 = sqrt(x.^2+z.^2)>t0*c_sound/2;
mask3 = sqrt(x.^2+z.^2)<time_axis(end)*c_sound/2;
mask = mask1 & mask2 & mask3;

reshaped_mask = reshape(mask,Nz,Nx);

gradientLayer = sqrt(x.^2 + z.^2);
gradientLayer = 1 + exp( gradientLayer./max(gradientLayer(:)) );

reshaped_gradientLayer = reshape(gradientLayer,Nz,Nx);

resampled_envelope_beamformed_data_masked = reshaped_gradientLayer.*reshaped_mask.*resampled_envelope_beamformed_data;

%%
% setting axis limits (mm)
x_lim = [min(x_matrixI(:)) max(x_matrixI(:))]*1e3;
z_lim = [min(z_matrixI(:)) max(z_matrixI(:))]*1e3;

% ploting image reconstruction


% compute dB values
env = resampled_envelope_beamformed_data_masked;
im = 20*log10(env./max(env(:) ));
vrange = [image_floor image_ceiling];

% display image
% imagesc(x_axis*1e3,z_axisI*1e3,im);
imagesc(x_axisI*1e3,z_axisI*1e3,min(max(im,image_floor),image_ceiling));
shading flat; colormap gray; caxis(vrange); colorbar; hold on;
axis equal manual;
xlabel('x [mm]');
ylabel('z [mm]');
set(gca,'YDir','reverse');
axis([x_lim z_lim]);
%}

%% Module 2
%{
S = 1:num_samples*2;
bfdataDAS = [];
for BEAM = 1:num_beams
    t0 = 103e-6;
    time_vector = (0:num_samples-1)/fs + t0;
    a_X = sin(thetas(BEAM)*pi/180);
    a_Z = cos(thetas(BEAM)*pi/180);
    X = Z_samples * S * a_X;
    Z = Z_samples * S * a_Z;
    beam_X(BEAM,:)=X;
    beam_Z(BEAM,:)=Z;
    transmit_delay = sqrt(X.^2 + Z.^2)/C;
    bfdataElem = 0;
    bflines = [];
    for E = 1:num_elements
        receive_delay = sqrt((elem_cntrs(E)-X).^2+Z.^2)/C;
        delay = transmit_delay + receive_delay;
        elemData = interp1(time_vector,frames(BEAM).RF_data(E,:),delay,'spline',0);
        bfdataElem = bfdataElem + elemData;
        bflines = [bflines; elemData];
        [BEAM E]
    end
    bfdataDAS = [bfdataDAS; bfdataElem];
end
BF_dataDAS = reshape(abs(hilbert(bfdataDAS(:,:)'))',size(bfdataDAS)).^2;

%% DAS
max_data   = max(BF_dataDAS(:));
temp_data  = BF_dataDAS/max_data+1e-10;
temp_data  = 10*log10(temp_data);
ROI.X_begin = -0.015;% ROI is extent of scatterers
ROI.X_end  = 0.015;
ROI.Z_begin = 0.005;%0.005;
ROI.Z_end  = 0.150;

% Scan Conversion

X_range = ROI.X_end-ROI.X_begin;
Y_range = maximum_depth;
image_width = 301;
image_height = floor(image_width*Y_range/X_range);

X_data = linspace (ROI.X_begin , ROI.X_end , image_width);
Y_data = linspace (0 , Y_range , image_height);

[XX,YY] = meshgrid(X_data , Y_data);
Recon_image = griddata( beam_X , beam_Z , temp_data , XX, YY);
Recon_image (isnan(Recon_image)) = image_floor;

figure, imagesc(X_data*10^3, Y_data*10^3, min(max(Recon_image,image_floor), image_ceiling));
colormap(gray);
colorbar;
axis image;
xlabel('WIDTH (mm)')
ylabel('DEPTH (mm)')
ylim([0 70])

%}
