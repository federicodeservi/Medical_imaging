%% MEDICAL IMAGING PROJECT
% Marco Peracchi & Federico De Servi, 2021

% CLUSTERING AND CLASSIFICATION

%Creates a csv containing features for each patient analyzed
root = 'C:\Users\feder\Desktop\Chest\manifest-1618047023244\Lung-PET-CT-Dx';

%path for third-party libraries
addpath(genpath('..\thirdparty-libraries'));
addpath(genpath('.\thirdparty-libraries'));
pathFeatures = strsplit(root, 'Lung-PET-CT-Dx');
pathFeatures = pathFeatures{1};

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

% Clustering using kmeans

% Kmeans on tot csv
TotFeaturesTable = readtable(strcat(pathFeaturesTot,'\Tot_Features\','tot_features.csv'));
TotFeatiresMatrix = table2array(TotFeaturesTable(:,1:6));

promptcl = [sprintf('Please enter the desired number of clsuters\n')];
namecl = 'Input cluster number';
n_clusters = str2double(inputdlg(promptcl,namecl,1,{'2'}));

idx=kmeans(TotFeatiresMatrix,n_clusters,'Replicates',1000);
%add cluster on non normalized matrix
TotFeaturesTableCluster = TotFeaturesTable;
TotFeaturesTableCluster.Cluster(:,1) = idx;

figcl = uifigure('Position',[500 500 840 360]);
uitcl = uitable(figcl,'Data',TotFeaturesTableCluster);
uitcl.Position = [20 20 800 320];
figcl.Color = 'white';

pause(5);

% Classification
labels = TotFeaturesTableCluster.Cluster;

%find most discriminative features
%First, standardize features
TotFeaturesTableClusterNorm = TotFeaturesTableCluster;
TotFeaturesTableClusterNorm = normalize(TotFeaturesTableClusterNorm,'zscore','DataVariables',{'volumesMm3' 'areas' 'Rs' 'spherical_disproportions' 'sphericities' 'surfacevolume__ratios' });

%LDA to find most discrim features
LDAClassifier = ClassificationDiscriminant.fit(TotFeaturesTableClusterNorm(:,1:6), labels, 'DiscrimType', 'linear');
LDAClassifier.DeltaPredictor

features = {'volumesMm3' 'areas' 'Rs' 'spherical_disproportions' 'sphericities' 'surfacevolume__ratios'};
discrim_power = LDAClassifier.DeltaPredictor;

idx_discrim =find(discrim_power > (max(discrim_power)/2));
top_features = features(idx_discrim);

%popup box
popuplda = strcat(sprintf('The most discriminative features are: \n \n'),sprintf('\n%s', string(top_features)));
resultlda = msgbox(popuplda);
resultlda.Color = 'white';
pause(2);

%start fitting
dlgTitle7    = 'User Question';
dlgQuestion7 = 'Do you want to fit the models?';
choice7 = questdlg(dlgQuestion7,dlgTitle7,'Yes','No', 'Yes');

if strcmpi(choice7, 'Yes')
    close(resultlda);

    %fit models 
    popupfit = sprintf('Model fitting...');
    resultfit = msgbox(popupfit);
    resultfit.Color = 'white';
    if n_clusters ==2
        %SVM model
        SVMModel = fitcsvm(TotFeaturesTableClusterNorm(:,idx_discrim),labels,'Standardize',true);
        CVSVMModel = crossval(SVMModel,'Holdout',0.35);
        svm_acc = 1 - kfoldLoss(CVSVMModel);

        %Classification tree
        tcModel = fitctree(TotFeaturesTableClusterNorm(:,idx_discrim),labels);
        CVtcModel = crossval(tcModel,'Holdout',0.35);
        ctree_acc = 1 - kfoldLoss(CVtcModel);

        close(resultfit);

        %Popup label 
        popupres = sprintf('SVM model accuracy: %d %% \nRandom Forest accuracy: %d %% \n', svm_acc*100, ctree_acc*100);
        resultres = msgbox(popupres);
        resultres.Color = 'white';
    else
        %Classification tree
        tcModel = fitctree(TotFeaturesTableClusterNorm(:,idx_discrim),labels);
        CVtcModel = crossval(tcModel,'Holdout',0.35);
        ctree_acc = 1 - kfoldLoss(CVtcModel);

        close(resultfit);

        %Popup label 
        popupres = sprintf('Random Forest accuracy: %d %% \n', ctree_acc*100);
        resultres = msgbox(popupres);
        resultres.Color = 'white';
    end
end

