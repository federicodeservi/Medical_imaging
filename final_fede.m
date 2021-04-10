%%
% Multiple with for loop
clear;
root = 'C:\Users\feder\Desktop\Chest\manifest-1618047023244\Lung-PET-CT-Dx';
% get the folder contents
d = dir(root);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

% which patient to display? 1, 2, 3...?
for k = 1 : length(nameFolds)
    cellContents = nameFolds{k};
    % Truncate and stick back into the cell
    patNames{k} = cellContents(9:end);
    number = num2str(k);
    patNames{k} = strcat('Patient', {' '},patNames{k}, {' '},' number: ',{' '}, number, '\n');
end
patNames = vertcat(patNames{:});
names2Prompt = strjoin(patNames, '\n');

prompt = [sprintf(names2Prompt), sprintf('\n'),sprintf('\n'), 'Please enter the number of the patient to be analyzed'];
name = 'Input patient number';
patNum = str2double(inputdlg(prompt,name,1,{'1'}));




patID = nameFolds{patNum};
pathMain = strcat(root,'\',patID);

d2 = dir(pathMain);
isub2 = [d2(:).isdir]; %# returns logical vector
nameFolds2 = {d2(isub2).name}';
nameFolds2(ismember(nameFolds2,{'.','..'})) = [];

pathMain2 = strcat(pathMain, '\', nameFolds2{1});

d3 = dir(pathMain2);
isub3 = [d3(:).isdir]; %# returns logical vector
nameFolds3 = {d3(isub3).name}';
nameFolds3(ismember(nameFolds3,{'.','..'})) = [];

% final path
path = strcat(pathMain2, '\', nameFolds3{1});
splitstring = split(root, '\Lung-PET-CT-Dx');
pathAnn = strcat(splitstring{1} , '\Annotation\',patID(9:end)) ;

filesxml = dir(fullfile(pathAnn, '*.xml'));
sampleXML = filesxml(1).name;
DOC = xmlread(fullfile(pathAnn, '\',sampleXML));
rootChildNodes = DOC.getChildNodes.item(0);
xmin = str2num(rootChildNodes.getElementsByTagName('xmin').item(0).getChildNodes.item(0).getNodeValue);
ymin = str2num(rootChildNodes.getElementsByTagName('ymin').item(0).getChildNodes.item(0).getNodeValue);
xmax = str2num(rootChildNodes.getElementsByTagName('xmax').item(0).getChildNodes.item(0).getNodeValue);
ymax = str2num(rootChildNodes.getElementsByTagName('ymax').item(0).getChildNodes.item(0).getNodeValue);
w = xmax-xmin;
h=ymax-ymin;
x = xmin;
y=ymin;


% how many imgs for patient?
n_imgs=9;
files = dir(fullfile(path,'*.dcm'));
addpath(genpath('C:\Users\feder\Google Drive\Universita\Materie\Medical Imaging\code__esempi\thirdparty-libraries'));
set(0,'defaultfigurecolor','black');
f = waitbar(0,'Please wait...', 'Position', [200,200,270,70]);
f.Color = 'white';

S.f = figure;

for n = 1:n_imgs
    waitbar((n/n_imgs),f,'Loading your data');
    info = dicominfo(fullfile(path,files(n*5).name), UseDictionaryVR=true);
    slope = info.RescaleSlope;
    intercept = info.RescaleIntercept;
    image = dicomread(fullfile(path,files(n*5).name));
    image = double(image);
    image = slope*image + intercept ;
    disp(n);
    ax = subaxis(3, 3, n, 'sh', 0, 'sv', 0.1, 'padding', 0.001, 'marginleft', 0.01, 'marginright', 0.01);
    %normalization
    min__ = min(min(image)); disp(min__);
    max__ = max(max(image)); disp(max__);
    image(image == min__) = -1000;
    imshow(image,"Parent", ax,'Colormap',gray);caxis('auto');
    if n==1
        title(sprintf('Frame %d',n), 'Color','white');
    else
        title(sprintf('Frame %d',(n-1)*5), 'Color','white');
    end
    hold on;
    rectangle('Position',[x, y, w, h],...
             'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');
    
end
close(f);

set(gcf, 'Position', get(0, 'Screensize'));



pause(2);

dlgTitle    = 'User Question';
dlgQuestion = 'Do you want to visualize all 50 frames?';
choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');

if strcmpi(choice, 'Yes')
    close(gcf);
    close(gcf);
    S2.f = figure;
    f = waitbar(0,'Please wait...', 'Position', [200,200,270,70]);
    f.Color = 'white';
    for n = 1:50
        waitbar((n/50),f,'Loading your data');
        info = dicominfo(fullfile(path,files(n).name), UseDictionaryVR=true);
        slope = info.RescaleSlope;
        intercept = info.RescaleIntercept;
        image = dicomread(fullfile(path,files(n).name));
        image = double(image);
        image = slope*image + intercept ;
        disp(n);
        ax = subaxis(5, 10, n, 'sh', 0, 'sv', 0.02, 'padding', 0.001, 'margintop', 0.025, 'marginbottom', 0.025,'marginleft', 0.01, 'marginright', 0.01);
        %normalization
        min__ = min(min(image)); disp(min__);
        max__ = max(max(image)); disp(max__);
        image(image == min__) = -1000;
        imshow(image,"Parent", ax,'Colormap',gray);caxis('auto');
        title(sprintf('Frame %d',n), 'Color','white');
        
        hold on;
        rectangle('Position',[x, y, w, h],...
                 'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');

    end     
    close(f);
    set(gcf, 'Position', get(0, 'Screensize'));
else
    close all;
end

pause(2);

dlgTitle2    = 'Do you want to visualize a specific frame (1 to 50)? Enter 0 to close and exit.';
dlgQuestion2 = 'Frame Number';
frameNum = str2double(inputdlg(dlgTitle2,dlgQuestion,1,{'1'}));

if frameNum==0
    S2.f = figure;
    close(S2.f);
elseif isempty(frameNum)==1
    S2.f = figure;
    close(S2.f);
else
    S2.f = figure;
    close(S2.f);
    info = dicominfo(fullfile(path,files(frameNum).name), UseDictionaryVR=true);
    slope = info.RescaleSlope;
    intercept = info.RescaleIntercept;
    image = dicomread(fullfile(path,files(frameNum).name));
    image = double(image);
    image = slope*image + intercept ;
    %normalization
    min__ = min(min(image)); disp(min__);
    max__ = max(max(image)); disp(max__);
    image(image == min__) = -1000;
    figure();
    imshow(image,'Colormap',gray);caxis('auto');
    title(sprintf('Frame %d',frameNum), 'Color','white');

    hold on;
    rectangle('Position',[x, y, w, h],...
             'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');
    set(gcf, 'Position', get(0, 'Screensize'));
end







%%
% SEGMENTATION AUTOMATIC 1

dlgTitle3    = 'From which frame do you want to start for segmentation?';
dlgQuestion3 = 'Frame Number';
frameNumSeg = str2double(inputdlg(dlgTitle3,dlgQuestion3,1,{'1'}));

contCond=1;
completedCond=0;
totVolume = 0;
while contCond==1
    if frameNumSeg > 50
       %popup 
        errormsg = sprintf('No such frame available. Select a lower number\n');
        h = msgbox(errormsg);
        contCond=0;
        
    else
        info = dicominfo(fullfile(path,files(frameNumSeg).name));
        slope = info.RescaleSlope;
        intercept = info.RescaleIntercept;
        img = dicomread(fullfile(path,files(frameNumSeg).name));
        img = double(img);
        img = slope*img + intercept;
        dimx = info.PixelSpacing(1);
        dimy = info.PixelSpacing(2);
        dimz = info.SliceThickness;
        dim__voxel = dimx*dimy*dimz;
        selected_img=img;
        figure();
        imshow(img, 'Colormap',gray);
        caxis('auto'); colorbar;
        title(sprintf('Frame %d',frameNumSeg), 'Color','white');
        set(gcf, 'Position', get(0, 'Screensize'));

        dlgTitle4    = 'User Question';
        dlgQuestion4 = 'Do you want to segment thia frame or proceed to the next one?';
        segChoice = questdlg(dlgQuestion4,dlgTitle4,'Segment','Proceed', 'Exit', "Exit");

        if strcmpi(segChoice, 'Segment')
            h = imfreehand;
            % Position of the ROI
            % x- and y-coordinates, respectively, of the n points along the boundary
            % of the freehand region
            pos = getPosition(h); 
            % Creating and visualizing an image that contains only the selected ROI
            mask = createMask(h);
            img__roi = selected_img.*mask;
            max__roi = max(max(img__roi));
            threshold = 0.8 * max__roi;
            img__roithreshold = img__roi;
            img__roithreshold(img__roithreshold > threshold) = 0;
            nvoxel = nnz(img__roithreshold > 0);

            % Calculate the volume, knowing the dimension of a single voxel in mm^3
            slice_volume = nvoxel * dim__voxel;
            totVolume = totVolume + slice_volume;
            completedCond=1;
            dlgTitle5    = 'User Question';
            dlgQuestion5 = 'Do you want to visualize and segment the following frame?';
            choice = questdlg(dlgQuestion5,dlgTitle5,'Yes','No', 'Yes');
            close(gcf);

            if strcmpi(choice, 'Yes')
                contCond=1;
            else
                contCond=0;
            end
            frameNumSeg = frameNumSeg+1;

        elseif strcmpi(segChoice, 'Proceed')
            close(gcf);
            frameNumSeg =frameNumSeg+1;
        else 
            contCond=0;
        end
    end
    
end

%disp(['Estimated total lesion volume: ' num2str(round(totVolume/1000))...
%    ' cc.']);
if completedCond ==1
    popupfinal = sprintf('Estimated total lesion volume: %d cc', (round(totVolume/1000)) );
    result = msgbox(popupfinal);
end
        
%%
% Reading the information about the patient and the PET acquisition
info = dicominfo(fullfile(path,files(17).name));

% Store slope and intercept
slope = info.RescaleSlope;
intercept = info.RescaleIntercept;

% Load a single slice of the 3-dimensional acquisition
img = dicomread(fullfile(path,files(17).name));
img = double(img);

% y = slope*x + intercept
img = slope*img + intercept;

% Visualize the image
figure();
imshow(img); caxis('auto'); colorbar
hold on;
    rectangle('Position',[x, y, w, h],...
             'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');
%%
% SEGMENTATION MANUAL 2

% Voxel dimension is saved in the fields PixelSpacing (x,y) and
% SliceThickness (z)
dimx = info.PixelSpacing(1);
dimy = info.PixelSpacing(2);
dimz = info.SliceThickness;

% Calculate the dimension of a single voxel [mm^3]
dim__voxel = dimx*dimy*dimz;

%%
% SEGMENTATION MANUAL 3

selected_img=img;

figure();
imshow(img, 'Colormap',gray);
caxis('auto'); colorbar;
h = imfreehand;

% Position of the ROI
% x- and y-coordinates, respectively, of the n points along the boundary
% of the freehand region
pos = getPosition(h); 

% Creating and visualizing an image that contains only the selected ROI
mask = createMask(h);
img__roi = selected_img.*mask;

figure();
imshow(img__roi,[],'Colormap',gray);
caxis('auto'); colorbar;

%%% con treshold
% Calculate the threshold
max__roi = max(max(img__roi));
threshold = 0.3 * max__roi;

% Apply the threshold to the ROI
img__roithreshold = img__roi;
img__roithreshold(img__roithreshold > threshold) = 0;
figure();
imshow(img__roithreshold,[],'Colormap',gray);

% Apply the threshold to the entire image
img__threshold = selected_img;
img__threshold(img__threshold > threshold) = 0;
figure();
imshow(img__threshold,[],'Colormap',gray);

% Modify the threshold
threshold = 0.3 * max__roi;
img__threshold = selected_img;
img__threshold(img__threshold > threshold) = 0;
figure();
imshow(img__threshold,[],'Colormap',gray);

% Calculate the volume of the selected ROI through the aplication of an
% intensity filter
% Calculate the number of voxels
nvoxel = nnz(img__roithreshold > 0);
% Alternative: nvoxel = sum(img__threshold(:) > 0);

% Calculate the volume, knowing the dimension of a single voxel in mm^3
volume = nvoxel * dim__voxel;
disp(['Il volume del cancro è pari a ' num2str(round(volume/1000))...
    ' cc.']);



%%

% NON CI INTERESSA PER ORA

for n = 1:50
        info = dicominfo(fullfile(path,files(n).name), UseDictionaryVR=true);
        slope = info.RescaleSlope;
        intercept = info.RescaleIntercept;
        image = dicomread(fullfile(path,files(n).name));
        image = double(image);
        image = slope*image + intercept ;
        stackimg(n,:,:) = image;
        disp(n);
        %normalization
  
        %img(img == min__) = -1000;


end   
    


min__ = min(min(min(stackimg)));
max__ = max(max(max(stackimg)));
stackimg(stackimg == min__) = -1000;


figure();
imshow(squeeze(stackimg(50,:,:)),'Colormap',gray);
caxis('auto'); colorbar

figure();
imshow(rot90(squeeze(stackimg(10,:,:)),0),'Colormap',gray);
caxis('auto'); colorbar
figure();
imshow(rot90(squeeze(stackimg(:,150,:)),0),'Colormap',gray);
caxis('auto'); colorbar
figure();
imshow(rot90(squeeze(stackimg(:,:,256)),0),'Colormap',gray);
caxis('auto'); colorbar
%%
%KMEANS (NON CI INTERESSA PER ORA)
% Use k-means clustering to perform automatic segmentation
xk = x -10;
yk = y-10;
wk = w + 10;
hk = h + 10;
imgcrop = imcrop(img,[xk yk wk hk]);

reshaped__img = reshape(imgcrop,[size(imgcrop,1)*size(imgcrop,2),1]);
idx = kmeans(reshaped__img,2,'replicate',5);
idx__reshaped = reshape(idx,[size(imgcrop,1),size(imgcrop,2)]);
figure; imshow(idx__reshaped); caxis('auto')

% We have now two images to display: the background image and the
% segmented ROI.
% In order to plot these two images as superimposed on the same figure,
% we should create a single image with values on two different color
% scales (in order to make them visually distinguishable).
% This image will then be plotted using two different colormaps stacked one
% over the other. Let us prepare a two-maps-in-one colormap by stacking a
% jet colormap and a hot colormap.
cmap__1 = colormap(gray);
cmap__2 = colormap(jet);
new__cmap = [cmap__1; cmap__2];

% Let us prepare the new image to be plotted. This new image will have the
% original-image values outside our thresholded ROI, and the new values
% inside the thresholded ROI.
% The new values will simply be the original values plus the maximum of the
% original image.
background__img = img./max(max(img));
background__level = mode(idx);
new__img = background__img;
for i = 1:size(img,1)
    for j = 1:size(img,2)
        if idx__reshaped(i,j) ~= background__level
            new__img(i,j) = 1.0001 + background__img(i,j);
        end
    end
end

% Plot the resulting image using the new colormap
imshow(new__img); colormap(new__cmap); caxis('auto'); colorbar

hold on;
    rectangle('Position',[x, y, w, h],...
             'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');















