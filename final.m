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
disp(patID)

popupSel = strcat('Patient selected: ',{' '}, patID(9:end));
selectionPat = msgbox(popupSel);
pause(2);
if ishandle(selectionPat)
    close(selectionPat);
end

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

%VISUALIZATIONS

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
        title(sprintf('Slice %d',n), 'Color','white');
    else
        title(sprintf('Slice %d',(n-1)*5), 'Color','white');
    end
    hold on;
    rectangle('Position',[x, y, w, h],...
             'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');
    
end
close(f);

set(gcf, 'Position', get(0, 'Screensize'));



pause(2);

dlgTitle    = 'User Question';
dlgQuestion = 'Do you want to visualize all 50 slices?';
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
        title(sprintf('Slice %d',n), 'Color','white');
        
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

dlgTitle2    = 'Do you want to visualize a specific slice (1 to 50)? Enter 0 to close and exit.';
dlgQuestion2 = 'Slice Number';
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
    title(sprintf('Slice %d',frameNum), 'Color','white');

    hold on;
    rectangle('Position',[x, y, w, h],...
             'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');
    set(gcf, 'Position', get(0, 'Screensize'));
end







%%
% SEGMENTATION AUTOMATIC 

dlgTitle3    = 'From which slice do you want to start for segmentation?';
dlgQuestion3 = 'Slice Number';
frameNumSeg = str2double(inputdlg(dlgTitle3,dlgQuestion3,1,{'1'}));

contCond=1;
completedCond=0;
totVolume = 0;
while contCond==1
    if frameNumSeg > 50
       %popup 
        errormsg = sprintf('No such slice available. Select a lower number\n');
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
        title(sprintf('Slice %d',frameNumSeg), 'Color','white');
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

if completedCond ==1
    popupfinal = sprintf('Estimated total lesion volume: %d cc', (round(totVolume/1000)) );
    result = msgbox(popupfinal);
end
        
