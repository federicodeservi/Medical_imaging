%% HN1280
% Reading the information about the patient and the PET acquisition
info = dicominfo('1-01', UseDictionaryVR=true);

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
% Load a single slice of the 3-dimensional acquisition
img = dicomread('1');
img = double(img);
intercept = info.RescaleIntercept;
% Multiply by the scaling factor
img = img*slope + intercept;



% "Fix" the visualization based on the output/color scale
minimum = min(min(img(img>0))); disp(num2str(minimum));
maximum = max(max(img)); disp(num2str(maximum));



% Normalize the output image between 0 and 1
img = img - minimum;
disp(num2str(max(max(img))));
img = img/maximum;
img = img*1;

% Plot the new image
figure();

imshow(img,'Colormap',hot); colorbar;

%%
% Multiple with for loop

% how many imgs for patient 1?
n_imgs=20
path = 'C:\Users\feder\Desktop\headneck\manifest-1568995398587\HEAD-NECK-RADIOMICS-HN1\HN1280\04-03-2019-PETCT-15018\31340';
files = dir(fullfile(path,'*.dcm'));
addpath(genpath('C:\Users\feder\Google Drive\Universita\Materie\Medical Imaging\code__esempi\thirdparty-libraries'));

for n = 1:n_imgs
    info = dicominfo(fullfile(path,files(n*5).name), UseDictionaryVR=true);
    slope = info.RescaleSlope;
    intercept = info.RescaleIntercept;
    image = dicomread(fullfile(path,files(n*5).name));
    image = double(image);
    image = slope*image + intercept ;
    disp(n);
    ax = subaxis(4, 5, n, 'sh', 0, 'sv', 0.01, 'padding', 0.01, 'margin', 0.01);
    %normalization
    minimum = min(min(img(img>0))); 
    maximum = max(max(img));
    img = img - minimum;
    img = img/maximum;
    img = img*1;
    imshow(image,"Parent", ax,'Colormap',hot);
end

%% 
% from 80 to end

n_imgs=6;

for n = 1:n_imgs
    info = dicominfo(fullfile(path,files(68+(1*n)).name), UseDictionaryVR=true);
    disp(fullfile(path,files(68+(1*n)).name));
    slope = info.RescaleSlope;
    intercept = info.RescaleIntercept;
    image = dicomread(fullfile(path,files(70+(1*n)).name));
    image = double(image);
    image = slope*image+ intercept ;
    disp(n);
    ax = subaxis(3,5, n, 'sh', 0, 'sv', 0.01, 'padding', 0.01, 'margin', 0.01);
    %normalization
    minimum = min(min(img(img>0))); 
    maximum = max(max(img));
    img = img - minimum;
    img = img/maximum;
    img = img*1;
    imshow(image,"Parent", ax,'Colormap',hot);
    caxis('auto');
end
