%% MEDICAL IMAGING PROJECT
% Marco Peracchi & Federico De Servi, 2021

% IMAGE VISUALIZATION TOOL

clear;
%Insert root folder
root = 'C:\Users\feder\Desktop\Chest\manifest-1618047023244\Lung-PET-CT-Dx';
%path for third-party libraries
addpath(genpath('C:\Users\feder\Google Drive\Universita\Materie\Medical Imaging\code__esempi\thirdparty-libraries'));
d = dir(root);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for k = 1 : length(nameFolds)
    cellContents = nameFolds{k};
    patNames{k} = cellContents(9:end);
    number = num2str(k);
    patNames{k} = strcat('Patient', {' '},patNames{k}, {' '},' number: ',{' '}, number, '\n');
end
patNames = vertcat(patNames{:});
names2Prompt = strjoin(patNames, '\n');

%Select which patient to display? 1, 2, 3...?
prompt = [sprintf(names2Prompt), sprintf('\n'),sprintf('\n'), 'Please enter the number of the patient to be analyzed'];
name = 'Input patient number';
patNum = str2double(inputdlg(prompt,name,1,{'1'}));

patID = nameFolds{patNum};
disp(patID)

%Patient selection confirmation popup
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

%View a 9-image overview
n_imgs=9;
files = dir(fullfile(path,'*.dcm'));
set(0,'defaultfigurecolor','black');
f = waitbar(0,'Please wait...', 'Position', [200,200,270,70]);
f.Color = 'white';

S.f = figure;

for n = 1:n_imgs
    waitbar((n/n_imgs),f,'Loading your data');
    info = dicominfo(fullfile(path,files(n*5).name), UseDictionaryVR=true);
    slope = info.RescaleSlope;
    intercept = info.RescaleIntercept;
    dimx = info.PixelSpacing(1);
    dimy = info.PixelSpacing(2);
    dimz = info.SliceThickness;
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
    %rectangle('Position',[x, y, w, h],...
    %         'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');
    
end
close(f);

set(gcf, 'Position', get(0, 'Screensize'));

%Pause 2 sec before asking to view all patient images
pause(2);

%All 50 image view
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
        %rectangle('Position',[x, y, w, h],...
        %         'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');

    end     
    close(f);
    set(gcf, 'Position', get(0, 'Screensize'));
else
    close all;
end

pause(2);

%Specific slice view
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
    %rectangle('Position',[x, y, w, h],...
    %         'LineWidth',0.2,'LineStyle','-', 'EdgeColor', 'r');
    set(gcf, 'Position', get(0, 'Screensize'));
end


% SEGMENTATION AUTOMATIC 

dlgTitle3    = 'From which slice do you want to start for segmentation?';
dlgQuestion3 = 'Slice Number';
frameNumSeg = str2double(inputdlg(dlgTitle3,dlgQuestion3,1,{'1'}));

contCond=1;
completedCond=0;
totVolume = 0;

pathLesionsDCMpart = strcat(root, '\',patID);
pathLesionsDCM = strcat(pathLesionsDCMpart,'\', 'Identified_lesions');

%Asks from which slice to begin segmentation. Then the doctor can segment
%it or proceed to the next. The doctor can continue the segmentation to the
%last slice that presents cancer. At the end the total cancer volume is
%presented.
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
            
            %Saves cancer segmentation as dicom image
            dicomwrite(img__roi, (sprintf('lesion_slice%d.dcm', (frameNumSeg))));
            if ~exist(strcat(pathLesionsDCM, 'dir'))
               mkdir(strcat(pathLesionsDCM))
            end
            filepathlesion = strcat(pwd, '\',(sprintf('lesion_slice%d.dcm', (frameNumSeg))));
            movefile(filepathlesion, strcat(pathLesionsDCM, '\', sprintf('lesion_slice%d.dcm', (frameNumSeg))));
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
        
%%
% FEATURE EXTRACTION

path = pathLesionsDCM;

%patID='Lung_Dx-A0160';
%path= strcat('C:\Users\feder\Desktop\Chest\manifest-1618047023244\Lung-PET-CT-Dx\', patID, '\Identified_lesions');

files = dir(fullfile(path,'*.dcm'));
PixelDimensions=[dimx, dimy, dimz];

sliceNum = [];
volumes=[];
areas=[];
Rs=[];
spherical_disproportions=[];
sphericities=[];
surfacevolume__ratios=[];

fs = waitbar(0,'Please wait. Saving info...', 'Position', [200,200,270,70]);
fs.Color = 'white';

for n = 1:size(files,1)
    waitbar((n/size(files,1)),fs,'Saving your data');
    % Reading the information about the patient and the PET acquisition
    info = dicominfo(fullfile(path,files(n).name));
    % Load a single slice of the 3-dimensional acquisition
    img = dicomread(fullfile(path,files(n).name));
    img = double(img);
    nameSlice = files(n).name;
    selSliceNum=regexp(nameSlice,'\d+','match');
    sliceNum=[sliceNum; str2num(selSliceNum{1})];
    
    % Volume of the lesion
    lesion__voxels = sum(img~=0,'all'); % Returns the number of voxels with non-zero values
    lesion__volume = lesion__voxels * dim__voxel;
    lesion__volume = double(lesion__volume);
    volumes=[volumes; lesion__volume];
    
    % Surface (area) of the lesion
    [lesion__area, surf_mat] = compute__surface(img, double(PixelDimensions));
    lesion__area = double(lesion__area);
    areas = [areas; lesion__area];
    
    % Spherical disproportion
    R_equiv = (lesion__volume*3/(4*pi))^(1/3);
    Rs = [Rs; R_equiv];
    
    spherical__disproportion = lesion__area/(4*pi*(R_equiv)^2);
    spherical_disproportions = [spherical_disproportions; spherical__disproportion];
    
    % Sphericity
    sphericity = ((pi^(1/3))*((6*lesion__volume)^(2/3)))/lesion__area;
    sphericities = [sphericities; sphericity];
    
    % Surface-to-volume ratio
    surfacevolume__ratio = lesion__area/lesion__volume;
    surfacevolume__ratios = [surfacevolume__ratios; surfacevolume__ratio];
    
end
close(fs);

%final table
features = [sliceNum, volumes, areas, Rs, spherical_disproportions, sphericities, surfacevolume__ratios];
T = array2table(features);
T.Properties.VariableNames(1:7) = {'sliceNum','volumes mm3', 'areas', 'Rs', 'spherical_disproportions', 'sphericities', 'surfacevolume__ratios' };
writetable(T,sprintf('%s_features.csv', patID(9:end)));

filepathfeatures = strcat(pwd, '\', sprintf('%s_features.csv', patID(9:end)));

pathFeatures = strsplit(root, 'Lung-PET-CT-Dx');
pathFeatures = pathFeatures{1};
if ~exist(strcat(pathFeatures,'\Features'), 'dir')
      mkdir(strcat(pathFeatures,'\Features'))
end
            
movefile(filepathfeatures, strcat(pathFeatures,'\Features'));


%% Creates a csv containing features for each patient analyzed

path =  strcat(pathFeatures,'\Features');
files = dir(fullfile(path,'*.csv'));

AllIDs = [];
AllVol=[];
AllArea=[];
AllR=[];
AllDispr=[];
AllSpher=[];
AllRatio=[];

for n = 1:size(files,1)
    filneameTemp = files(n).name;
    Temp = readtable(fullfile(path,filneameTemp));
    disp(Temp);
    TempID = files(n).name(1:5);
    
    VolSum = sum(table2array(Temp(:,2)));
    AreaMean = mean(table2array(Temp(:,3)));
    RMean = mean(table2array(Temp(:,4)));
    DisprMean = mean(table2array(Temp(:,5)));
    SpherMean = mean(table2array(Temp(:,6)));
    RatioMean = mean(table2array(Temp(:,7)));
    
    AllIDs = [AllIDs; TempID];
    AllVol = [AllVol; VolSum];
    AllArea = [AllArea; AreaMean];
    AllR = [AllR; RMean];
    AllDispr = [AllDispr; DisprMean];
    AllSpher = [AllSpher; SpherMean];
    AllRatio = [AllRatio; RatioMean];

end


TotFeatures = [AllVol, AllArea, AllR, AllDispr, AllSpher, AllRatio];
TempTot = array2table(TotFeatures);
TempTot.newVar(:,1) = cellstr(AllIDs);
TempTot.Properties.VariableNames(1:7) = {'volumes mm3', 'areas', 'Rs', 'spherical_disproportions', 'sphericities', 'surfacevolume__ratios' ,'PatID'};

writetable(TempTot,'tot_features.csv');

filepathfeaturestot = strcat(pwd, '\', 'tot_features.csv');

pathFeaturesTot = strsplit(root, 'Lung-PET-CT-Dx');
pathFeaturesTot = pathFeaturesTot{1};
if ~exist(strcat(pathFeaturesTot,'\Tot_Features'), 'dir')
      mkdir(strcat(pathFeaturesTot,'\Tot_Features'))
end
            
movefile(filepathfeaturestot, strcat(pathFeaturesTot,'\Tot_Features'));

%% Clustering using kmeans

% Kmeans on tot csv
TotFeaturesTable = readtable(strcat(pathFeaturesTot,'\Tot_Features\','tot_features.csv'));
TotFeatiresMatrix = table2array(TotFeaturesTable(:,1:6));

promptcl = [sprintf('Please enter the desired number of clsuters\n')];
namecl = 'Input cluster number';
n_clusters = str2double(inputdlg(promptcl,namecl,1,{'2'}));

idx=kmeans(TotFeatiresMatrix,n_clusters,'Replicates',1000);
TotFeaturesTableCluster = TotFeaturesTable;
TotFeaturesTableCluster.Cluster(:,1) = idx;

figcl = uifigure('Position',[500 500 840 360]);
uitcl = uitable(figcl,'Data',TotFeaturesTableCluster);
uitcl.Position = [20 20 800 320];
figcl.Color = 'white';

%% Classification

%find median volume
min_vol = min(table2array(TotFeaturesTable(:,1)));
max_vol = max(table2array(TotFeaturesTable(:,1)));
median_treshold = median(table2array(TotFeaturesTable(:,1)));

%Assing label 0 if patient volume < median volume, 1 otherwise
mask_low = TotFeaturesTableCluster.volumesMm3 < median_treshold;
TotFeaturesTableCluster.Label(mask_low) = '0';
mask_high = TotFeaturesTableCluster.volumesMm3 >= median_treshold;
TotFeaturesTableCluster.Label(mask_high) = '1';

%Popup label 
figcl = uifigure('Position',[500 500 940 360]);
uitcl = uitable(figcl,'Data',TotFeaturesTableCluster);
uitcl.Position = [20 20 900 320];
figcl.Color = 'white';


%Model training
labels = TotFeaturesTableCluster.Label;
SVMModel = fitcsvm(TotFeatiresMatrix,labels,'Standardize',true);
CVSVMModel = crossval(SVMModel,'Holdout',0.3);
supervised__accuracy = 1 - kfoldLoss(CVSVMModel);

%% Additional plots




