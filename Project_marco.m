%-------------------------------------------------------------------------
% Practical session #2a | #dicom #ct
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% 1. Inizialization
%-------------------------------------------------------------------------
clear; clc; close all;
addpath(genpath('../data'));
addpath(genpath('../thirdparty-libraries'));

% Create a log file
diary('log.txt');
diary on


%-------------------------------------------------------------------------
% 2. Working with DICOM images
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% 2.1 Loading and visualizing a single DICOM image
%-------------------------------------------------------------------------
cd '../data/single-patient/CT/dicom';

% Reading the information about the patient and the PET acquisition
info = dicominfo('1-100');

% Printing the name of the patient
info.PatientName

% Visualizing and saving information about image/acquisition
% Information to be considered:
% 1. Voxel dimension (x, y, z) -> useful for computing a lesion volume
% 2. Scale factor between a value visualized and the effective number 
% of counts -> useful for computing the uptake in a specific ROI

% Voxel dimension is saved in the fields PixelSpacing (x,y) and
% SliceThickness (z)
dimx = info.PixelSpacing(1);
dimy = info.PixelSpacing(2);
dimz = info.SliceThickness;

% Calculate the dimension of a single voxel [mm^3]
dim__voxel = dimx*dimy*dimz;

% The scale factor is saved in the field called "RescaleSlope"
slope = info.RescaleSlope;
intercept = info.RescaleIntercept;

% Load a single slice of the 3-dimensional acquisition
img = dicomread('1-100');
img = double(img);

% Multiply by the scaling factor
% IMPORTANT: in this case the intercept is not equal to zero, so we must
% consider it when comuting the img values (img = intercept + img*slope)
img = intercept + img*slope;

% Visualize the image
figure();
imshow(img);

% We can "fix" the visualization based on the output/color scale by
% adjusting the colormap on the required range of values. The following
% command does the task automatically. However, you can modify the
% colormap according to your needs
caxis('auto')

% Visualize the colorbar
colorbar

% IMPORTANT: the "black" pixels in the figure (near the corners of the
% image) of value ~ -3000 seem to fall outside the Hounsfield scale for CT
% scans. These pixels are set to this value by the scanning machine to be
% recognizable, as they fall outside the (circular) field of view of the CT
% scan (as it can be seen in the figure). In order to improve the
% visualization of the CT scan (i.e., in order to correctly fit the
% colormap on the CT-values range), you can easily put these values to
% the lowest value of the Hounsfield scale.
min__ = min(min(img)); disp(min__);
max__ = max(max(img)); disp(max__);
img(img == min__) = -1000;
figure(); imshow(img); caxis('auto'); colorbar


%-------------------------------------------------------------------------
% 2.2 Loading the entire 3-dimensional volume as a DICOM image
%-------------------------------------------------------------------------
files = dir('*.dcm');

for n = 1:size(files,1)
    pathtemp = files(n).name;
    disp(pathtemp);
    info__temp = dicominfo(pathtemp);
    dcm__temp = dicomread(pathtemp);
    dcm__temp = double(info__temp.RescaleIntercept + ...
        dcm__temp.*info__temp.RescaleSlope);
    stackimg(n,:,:) = dcm__temp;
    disp(['The dimension of stackimg after ' num2str(n)...
        ' iterations is ' num2str(size(stackimg))]);
end
nvoxel = size(stackimg,2) * size(stackimg,3);
disp(['Each one of the ' num2str(size(stackimg,1))...
    ' loaded images is made of ' num2str(nvoxel) ' voxels']);

min__ = min(min(min(stackimg)));
max__ = max(max(max(stackimg)));
stackimg(stackimg == min__) = -1000;

% Visualize a single slice of the entire volume
figure();
imshow(squeeze(stackimg(17*3,:,:)),'Colormap',gray);
caxis('auto'); colorbar
% Provo a visualizzare con calma tutte le sezioni
% Quelli interessanti sono i 17*3 e 18*3!!!
for i=1:30
    imshow(rot90(squeeze(stackimg(i*3,:,:)),0),'Colormap',gray);
    caxis('auto'); colorbar
    pause(2)
    set(gcf, 'Name', num2str(i*10), 'NumberTitle', 'Off');
end 

%Mio test
figure();
imshow(rot90(squeeze(stackimg(:,256,:)),0),'Colormap',gray);
caxis('auto'); colorbar
figure();
imshow(rot90(squeeze(stackimg(:,270,:)),0),'Colormap',gray);
caxis('auto'); colorbar
figure();
imshow(rot90(squeeze(stackimg(:,:,256)),0),'Colormap',gray);
caxis('auto'); colorbar


%-------------------------------------------------------------------------
% 3. Closing the script
%-------------------------------------------------------------------------

% Interrupt the input to the log file
close all
clear
clc
diary off