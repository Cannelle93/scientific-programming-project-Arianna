
load('Data')

%% Cleaning the Dataset
% since the data represents genes expression, i decided to delete the
% missing values (zeros). In fact a zero value means that the gene hasn't
% been expressed in that specific tissue (brain).

%delete missing data 
x=1;

for i = 1:size(geneEx,2)
    x = x.* geneEx(:,i);
    
end

n = find(x); % finding all the values that are non zero
length(n); 
geneExN = geneEx(n,:); % generate a new dataset without missing values
GeneIDN = GeneID(n);   % list with gene identifiers 
labels = [zeros(1,6) ones(1,6)]; % first 6 samples are offspings of unstressed mothers and the last 6 samples are offsprings of stressed mothers

%% Data Cenetring and Transformation
% Log2 transformation and centering 
%  MeanGene = mean(log(geneExN),2); %compute mean of the log of the data
%  NormData = log(geneExN) - MeanGene; % substract mean to center de dataset 

% Further in this code a volcanoplot will be done, where the negative
% values are neglected. Therfore the dataset shouldn't be centered

% Only log2 trasformation without centering
NormData = log(geneExN);
figure(1)
hist(NormData(:),101 ) %histogram of cleaned and scaled Dataset 


%% PCA
exVar=99;% Explained Variance 
[coeff,score,~,~,explained] = pca(NormData');
idx = find(cumsum(explained)>= exVar,1);
PCAData = score(:,1:idx);


%% PCA plot
figure(2)
plot(score(1:3,1),score(1:3,2),'bo'),hold on,
plot(score(4:6,1),score(4:6,2),'ro'),
plot(score(7:9,1),score(7:9,2),'b+'),
plot(score(10:12,1),score(10:12,3),'r+'),
legend('UM','UF','SM','SF'), hold on
title( 'PCA plot'),
ylabel('component 2    9,04% ')
xlabel('component 1    58,57 %')

hold off

%% Oversampling
%% over sampling by SMOTE technique 
% smote technic allows to generate new synthetic samples between the
% already existing samples by using k-nearest neighbour
%initiation
 UnstressedMale = [];
 UnstressedFemale = [];
 StressedMale = [];
 StressedFemale = [];
 Unstressed = [];
 Stressed = [];
 
 %gender separation
 % Unstressed males
for i=1:10
in_labels = [ones(1,3),zeros(1,9)];
[out_featuresSyn, out_labelsSyn] = ADASYN(NormData', in_labels,1,2,2,true);
UnstressedMale((i-1)*6+1:i*6,:)= out_featuresSyn;

end

%Unstressed females
for i=1:10
in_labels = [zeros(1,3),ones(1,3),zeros(1,6)];

[out_featuresSyn, out_labelsSyn] = ADASYN(NormData', in_labels,1,2,2,true);
UnstressedFemale((i-1)*6+1:i*6,:)= out_featuresSyn;

end

%Stressed males
for i=1:10
in_labels = [zeros(1,6),ones(1,3),zeros(1,3)];

[out_featuresSyn, out_labelsSyn] = ADASYN(NormData', in_labels,1,2,2,true);
StressedMale((i-1)*6+1:i*6,:)= out_featuresSyn;

end

%Stressed females
for i=1:10
in_labels = [zeros(1,9),ones(1,3)];
[out_featuresSyn, out_labelsSyn] = ADASYN(NormData', in_labels,1,2,2,true);
StressedFemale((i-1)*6+1:i*6,:)= out_featuresSyn;

end
%% PCA
% concatenate all the generated data using the smote technic
SmoteData = cat(1,UnstressedMale,UnstressedFemale);
SmoteData = cat(1, SmoteData,StressedMale);
SmoteData = cat( 1,SmoteData,StressedFemale);

exVar=99;% Explained Variance 
[coeff,score,~,~,explained] = pca(SmoteData);
idx = find(cumsum(explained)>= exVar,1);
SmotePCA= score(:,1:idx);

% PCA plot

% figure(3)
% plot3(score(1:60,1),score(1:60,2),score(1:60,4),'bo'),hold on,
% plot3(score(61:120,1),score(61:120,2),score(61:120,3),'ro'),
% plot3(score(121:180,1),score(121:180,2),score(121:180,3),'b+'),
% plot3(score(181:240,1),score(181:240,2),score(181:240,3),'r+'),
% 
% legend('UM','UF','SM','SF')
% title( 'PCA plot'),
% ylabel('component 2')
% xlabel('component 1')
% zlabel('component 3')
% hold off

figure(3)
plot(score(1:60,1),score(1:60,2),'bo'),hold on,
plot(score(61:120,1),score(61:120,2),'ro'),
plot(score(121:180,1),score(121:180,2),'b+'),
plot(score(181:240,1),score(181:240,2),'r+'),
legend('UM','UF','SM','SF'), hold on
title( 'PCA plot Smote with gender separation'),
ylabel('component 2    10,48% ')
xlabel('component 1    54,66 %')

hold off



 %% without gender separation 
 
 % since we have 6 samples from stressed mothers and 6 controls the dataset
 % is balanced. However the ADASYN function applyed before is used to
 % oversample de minority class. Therfore the function has been modified (ADASYN2) so
 % that one could choose wich class is ment to be oversampled eventhough it
 % doesn't represent the minority class. 
 
in_beta=20;
in_kDensity=5;
in_kSMOTE=5;
in_featuresAreNormalized=false;
[out_featuresSyn, out_labelsSyn] = ADASYN2(NormData', labels==1, in_beta, in_kDensity, in_kSMOTE, in_featuresAreNormalized,1,0);
Stressed = out_featuresSyn ;
 [out_featuresSyn, out_labelsSyn] = ADASYN2(NormData', labels==1, in_beta, in_kDensity, in_kSMOTE, in_featuresAreNormalized,0,1);
 Unstressed = out_featuresSyn;

 % concatenated de genertaed datasets
SmoteData2 = cat(1,Unstressed, Stressed);

% PCA
exVar=99;% Explained Variance (100 no data reduction)
[~,score,~,~,explained] = pca(SmoteData2);
idx = find(cumsum(explained)>= exVar,1);
SmotePCA2= score(:,1:idx);
 
% PCA plots
figure(4)
plot(score(1:120,1),score(1:120,2),'go'),hold on,
plot(score(121:240,1),score(121:240,2),'g+'),

legend('U','S')
title( 'PCA plot Smote NO gender separation'),
ylabel('component 2    8,95% ')
xlabel('component 1    55,5 %')

hold off

%% Support vector machine

 % support vector machine males 
X_train = [UnstressedMale;StressedMale];
c_train = [zeros(size(UnstressedMale,1),1);ones(size(StressedMale,1),1)];
NormDatat = NormData';
X_test = NormDatat([1:3,7:9],:);
c_test = [zeros(1,3),ones(1,3)];
SVMModel = fitcsvm(X_train,c_train,'KernelFunction','linear', 'BoxConstraint',20,'KernelScale','auto');

nbSVMM = predict(SVMModel,X_test);
AccMales = sum(nbSVMM' == c_test)/length(c_test)


% support vector machine females
X_train = [UnstressedFemale;StressedFemale];
c_train = [zeros(size(UnstressedFemale,1),1);ones(size(StressedFemale,1),1)];
X_test = NormDatat([4:6,10:12],:);
c_test = [zeros(1,3),ones(1,3)];
SVMModel = fitcsvm(X_train,c_train,'KernelFunction','linear', 'BoxConstraint',20,'KernelScale','auto');

nbSVMM = predict(SVMModel,X_test);
AccFemales = sum(nbSVMM' == c_test)/length(c_test)


% support vector no gender separation
X_train = [Unstressed;Stressed];
c_train = [zeros(size(Unstressed,1),1);ones(size(Stressed,1),1)];
X_test = NormDatat([1:6,7:12],:);
c_test = [zeros(1,6),ones(1,6)];
SVMModel = fitcsvm(X_train,c_train,'KernelFunction','linear', 'BoxConstraint',20,'KernelScale','auto');

nbSVMM = predict(SVMModel,X_test);
AccT = sum(nbSVMM' == c_test)/length(c_test)

% For the crossvalidation the already existing samples have been used as
% test set whereas the newly generated samples have been used as training
% set
% The support vector machine gave an accuracy of 100%, which shows that
% the synthetic samples are strongly depedent on the real samples.
% Even tough one can't use SVM for prediction in this case, it is a nice
% idicator that the oversampling worked well. 

%% oversampling with bootstrapping technic
% bootstrapping is an oversampling technic that is based on random
% permutations


% with gender separation
 
 bootstrapUM = datasample(NormDatat(1:3,:),60);  %unstressed males
 bootstrapUF = datasample(NormDatat(4:6,:),60);  %unstressed females

 bootstrapSM = datasample(NormDatat(7:9,:),60);  %stressed males
 bootstrapSF = datasample(NormDatat(10:12,:),60);%stressed females
 
 % concatenating matrixes
 bootstrapData = cat(1,bootstrapUM,bootstrapUF);
 bootstrapData = cat(1,bootstrapData,bootstrapSM);
 bootstrapData = cat(1,bootstrapData,bootstrapSF);
 
 %PCA
 exVar=99;% Explained Variance 
[coeff,score,~,~,explained] = pca(bootstrapData);
idx = find(cumsum(explained)>= exVar,1);
bootstrapPCA= score(:,1:idx);
 
% PCA plot
figure(6)
plot(score(1:60,1),score(1:60,2),'bo'),hold on,
plot(score(61:120,1),score(61:120,2),'ro'),
plot(score(121:180,1),score(121:180,2),'b+'),
plot(score(181:240,1),score(181:240,2),'r+'),
legend('UM','UF','SM','SF')

title( 'PCA plot Bootstrapping with gender separation'),
ylabel('component 2    9,08% ')
xlabel('component 1    58,33 %')

hold off

 % without gender separation

 bootstrapU = datasample(NormDatat(1:6,:),120); %unstressed
 bootstrapS = datasample(NormDatat(7:12,:),120);%streesed
 
 % concatenating matrixes
 bootstrapData2 = cat(1,bootstrapU, bootstrapS);
 
 %PCA
 exVar=99;% Explained Variance (100 no data reduction)
[coeff,score,~,~,explained] = pca(bootstrapData2);
idx = find(cumsum(explained)>= exVar,1);
bootstrapPCA2= score(:,1:idx);
 
% PCA plot
figure (5)
plot(score(1:120,1),score(1:120,2),'go'),hold on,
plot(score(121:240,1),score(121:240,2),'g+'),
legend('U','S')

title( 'PCA plot Bootstrapping NO gender separation'),
ylabel('component 2    8% ')
xlabel('component 1    59,9 %')

hold off

 
%% SVM for bootstrapping
% SVM for Males
X_train = [bootstrapUM;bootstrapSM];
c_train = [zeros(size(bootstrapUM,1),1);ones(size(bootstrapSM,1),1)];
X_test = NormDatat([1:3,4:6],:);
c_test = [zeros(1,3),ones(1,3)];
SVMModel = fitcsvm(X_train,c_train,'KernelFunction','linear', 'BoxConstraint',20,'KernelScale','auto');

nbSVMM = predict(SVMModel,X_test);
AccBM = sum(nbSVMM' == c_test)/length(c_test)

%SVM for Females
X_train = [bootstrapUF;bootstrapSF];
c_train = [zeros(size(bootstrapUF,1),1);ones(size(bootstrapSF,1),1)];
X_test = NormDatat([4:6,10:12],:);
c_test = [zeros(1,3),ones(1,3)];
SVMModel = fitcsvm(X_train,c_train,'KernelFunction','linear', 'BoxConstraint',20,'KernelScale','auto');

nbSVMM = predict(SVMModel,X_test);
AccBF = sum(nbSVMM' == c_test)/length(c_test)

% without gender separation
X_train = [bootstrapU;bootstrapS];
c_train = [zeros(size(bootstrapU,1),1);ones(size(bootstrapS,1),1)];
X_test = NormDatat([1:6,7:12],:);
c_test = [zeros(1,6),ones(1,6)];
SVMModel = fitcsvm(X_train,c_train,'KernelFunction','linear', 'BoxConstraint',20,'KernelScale','auto');

nbSVMM = predict(SVMModel,X_test);
AccBT = sum(nbSVMM' == c_test)/length(c_test)

% Same as before the SVM only indicated that the oversampling with
% bootstrapping worked well

%% Dataset that has been oversampled considering geneder separation but will be devided only between streesed and unstressed for the analysis 
% Total smote data with gender separation 
TotalSmoteU = cat(1,UnstressedMale,UnstressedFemale);
TotalSmoteS = cat(1,StressedMale,StressedFemale);

% Total bootstrap data with gender separation 
bootstrapUT = cat(1,bootstrapUM, bootstrapUF);
bootstrapST = cat(1,bootstrapSM, bootstrapUF);


%% wilcoxon rank sum test
% Smote Males
for i= 1:17537
p_SM (i) = ranksum(UnstressedMale(:,i),StressedMale (:,i));
end
% Smote Females
for i= 1:17537
p_SF (i) = ranksum(UnstressedFemale(:,i),StressedFemale (:,i));
end
% Total smote with gender separation 
for i= 1:17537
p_TSGS (i) = ranksum(TotalSmoteU(:,i),TotalSmoteS(:,i));
end
% Total Smote no gender separation
for i= 1:17537
p_ST (i) = ranksum(Unstressed(:,i),Stressed (:,i));
end

% bootstrap Males
for i= 1:17537
p_BM (i) = ranksum(bootstrapUM(:,i),bootstrapSM (:,i));
end
% bootstrap Females
for i= 1:17537
p_BF (i) = ranksum(bootstrapUF(:,i),bootstrapSF (:,i));
end
%bootstrap total with gender separation 
for i= 1:17537
p_BTGS (i) = ranksum(bootstrapUT(:,i),bootstrapST (:,i));
end
% bootstrap no gender separation
for i= 1:17537
p_BT (i) = ranksum(bootstrapU(:,i),bootstrapS (:,i));
end
%% adjusting p_value with mafdr
adjustedP_SM = mafdr(p_SM);        %smote males
adjustedP_SF = mafdr(p_SF);        %smote females
adjustedP_TSGS = mafdr(p_TSGS);    %smote Total with gender separation
adjustedP_ST = mafdr(p_ST);        %smote Total 
adjustedP_BM = mafdr(p_BM);        %Bootstrapping males
adjustedP_BF = mafdr(p_BF);        %Bootstrapping femals
adjustedP_BTGS = mafdr(p_BTGS);    %Bootstrapping Total with gender separation
adjustedP_BT = mafdr(p_BT);        %Bootstrapping Total

%% vulcano plot

% Smote males
SMStructure = mavolcanoplot(UnstressedMale', StressedMale',adjustedP_SM','Labels', GeneIDN,'Foldchange',3.5)
% Smote females
SFStructure = mavolcanoplot(UnstressedFemale',StressedFemale',adjustedP_SF','Labels', GeneIDN,'Foldchange',3.5)
% Total smote with gender separation
TSGSStructure = mavolcanoplot(TotalSmoteU',TotalSmoteS',adjustedP_TSGS','Labels', GeneIDN,'Foldchange',3)
% Totla smote no gender separation
STStructure = mavolcanoplot(Unstressed',Stressed',adjustedP_ST','Labels', GeneIDN,'Foldchange',3)

% Bootstrapping males
BMStructure = mavolcanoplot(bootstrapUM', bootstrapSM',adjustedP_BM','Labels', GeneIDN,'Foldchange',3.5)
% Bootstrapping females
BFStructure = mavolcanoplot(bootstrapUF', bootstrapSF',adjustedP_BF','Labels', GeneIDN,'Foldchange',3.5)
% Total bootstrapping with gender separation 
BTGSStructure = mavolcanoplot(bootstrapUT', bootstrapST',adjustedP_BTGS','Labels', GeneIDN,'Foldchange',2)
% Total bootstrapping no gender separation 
BTStructure = mavolcanoplot(bootstrapU', bootstrapS',adjustedP_BT','Labels', GeneIDN,'Foldchange',2.5)

% For the volcanoplots for the Total dataset a smaller Foldchange has been
% selected otherwise we had no genes that were above respectively below the
% set threshold 

%% find over and under expressed genes 
for i = 1:length(SMStructure.PValues)
SM(i)=find(SMStructure.PValues(i) == adjustedP_SM );
end
for i = 1:length(SFStructure.PValues)
SF(i)=find(SFStructure.PValues(i) == adjustedP_SF );
end
for i = 1:length(TSGSStructure.PValues)
TSGS(i)=find(TSGSStructure.PValues(i) == adjustedP_TSGS );
end
 for i = 1:length(STStructure.PValues)
ST(i)=find(STStructure.PValues(i) == adjustedP_ST );

 for i = 1:length(BMStructure.PValues)
BM(i)=find(BMStructure.PValues(i) == adjustedP_BM );
 end

% All under respectivaly over expressed genes identified in SmoteData
TotalS = cat(2,SM,SF);
TotalS = cat(2,TotalS,TSGS);
TotalS = cat(2,TotalS,ST);
TotalS = unique(TotalS) 


 end
 for i = 1:length(BFStructure.PValues)
BF(i)=find(BFStructure.PValues(i) == adjustedP_BF );
 end
 for i = 1:length(BTGSStructure.PValues)
BTGS(i)=find(BTGSStructure.PValues(i) == adjustedP_BTGS );
 end
 for i = 1:length(BTStructure.PValues)
BT(i)=find(BTStructure.PValues(i) == adjustedP_BT );
 end
 
 
% All under respectivaly over expressed genes identified in Data generated
% by bootstrapping
TotalB = cat(2,BM,BF);
TotalB = cat(2, TotalB,BTGS);
TotalB = cat(2,TotalB,BT);
TotalB = unique(TotalB)

% Find same genes over/under expressed genes in both datasets 
Total = find(ismember(TotalS,TotalB));



%% heatmap
HY= GeneIDN(Total);
HX= ['SM';'BM';'SF';'BF';'SA';'BA';'ST';'BT'];
HData = [adjustedP_SM(Total);adjustedP_BM(Total);adjustedP_SF(Total);adjustedP_BF(Total);adjustedP_TSGS(Total);adjustedP_BTGS(Total);adjustedP_ST(Total);adjustedP_BT(Total)]
HData = log(HData')
%recentering HData
meanHDATA = mean(HData,2)
HData= HData- meanHDATA
h = heatmap(HX,HY,HData,'Title','Heatmap','Colormap',parula);

%  %% heatmapsmote
% HY= GeneIDN(TotalS);
% HX= ['SM';'SF';'SA';'ST'];
% HData = [adjustedP_SM(TotalS);adjustedP_SF(TotalS);adjustedP_TSGS(TotalS);adjustedP_ST(TotalS)]
% HData = log(HData')
% %recentering HData
% meanHDATA = mean(HData,2)
% HData= HData- meanHDATA
% h = heatmap(HX,HY,HData,'Title','HeatmapS','Colormap',parula);


%  %% heatmapbootstrapping
% HY= GeneIDN(TotalB);
% HX= ['BM';'BF';'BA';'BT'];
% HData = [adjustedP_BM(TotalB);adjustedP_BF(TotalB);adjustedP_BTGS(TotalB);adjustedP_BT(TotalB)]
% HData = log(HData');
% %recentering HData
% meanHDATA = mean(HData,2)
% HData= HData- meanHDATA
% h = heatmap(HX,HY,HData,'Title','HeatmapB','Colormap',parula);


%% Biplot
exVar=99;% Explained Variance 
[coeff,score,~,~,explained] = pca(SmoteData);
idx = find(cumsum(explained)>= exVar,1);
SmotePCA= score(:,1:idx);

 figure(11)
biplot(coeff(SM,1:2),'Scores',score(1:60,1:2),'Color','b','Marker','o','VarLabels',GeneIDN(SM)), hold on
biplot(coeff(SF,1:2),'Scores',score(61:120,1:2),'Color','r','Marker','o','VarLabels',GeneIDN(SF)), hold on
biplot(coeff(SM,1:2),'Scores',score(121:180,1:2),'Color','b','Marker','+'), hold on
biplot(coeff(SF,1:2),'Scores',score(181:240,1:2),'Color','r','Marker','+'), hold off


% Biplot is difficult to interprete since the under respectively the over
% expressed genes are still too many and due to the smote technic the
% samples of the same class are quasi alined in the biplot
