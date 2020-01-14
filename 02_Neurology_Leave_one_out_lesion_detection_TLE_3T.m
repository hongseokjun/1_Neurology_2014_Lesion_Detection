%% Hong SJ et al, 'Automated detection of cortical dysplasia type II in MRI-negative epilepsy' Neurology, 2014
%  0) The orginal idea was from the 2008 MICCAI paper by Dr. Pierre Besson
%  1) The lesion classifier works on the z-score of features, thus the z-score calculation should be proceeded before carrying out this script
%  2) The algorithm consists of two-step classification: 1> vertex-wise classifier and 2> cluster-wise classifier
%     vertex-wise classifier: First, we randomly select same number of the vertices from non-lesional area as lesions for balanced sample size.                             
%                             Second, we make vertex-wise feature vectors, each of which consists of five surface-based features (RI, PG, CT, MC, SD).
%                             We then feed them into linear discriminant classifier to train first and test a new case (kernel: diagquadratic, which gives highest sensitivity)
%     cluster-wise classifier: We define the cluster by the condition of 6-neighbouring connectivity and non-involvement of mid-line masks
%                              We extract statistical moments (e.g., mean, SD, skewness) of each feature from TP and FP clusters
%                              and finally feed them into the second cluster-wise classifier (kernel: linear, which gives highest specificity).
%  3) This script is for estimating the specificity of the classification in disease controls, TLE with histopathologically confirmed HS, long-term seizure free
%  4) We use the classifier trained using 3T data of FCD patients
%  5) Variable description
%     TemplateSurf: template surface that defines the number of vertices & triangles
%     featVect: feature Vectors, l x m x n matrix, size of l: #vertices, m: #feature
%     types (i.e., thickness, curvature, RI, gradient, etc), n: #FCD patients
%     lesionMasks: lesional masks for all patients, l x n: size of l:#vertices,
%     Please assign FCD lesional vertices as 1 otherwise 0 (i.e., healthy tissues).
%     If you like to exclude any verices from sampling, assign them 2.

% original version of this file is at /local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/lesion_detection/org_script/

clear all;
addpath('/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/lesion_detection/utility');

%% Results used in the paper. For convenience, you can always re-load all the result, so as to reduce the time of all the re-computation 
load('/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/lesion_detection/Result_TLE_3T.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. variable initialization - read files
Group_TLE                = 'TLE';
Prefix_TLE               = 'TLE';
Cases_TLE                = '/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/Demographic_data_TLE.txt';
Group_pat                 = 'FCD';
Prefix_pat                = 'mcd';
Cases_pat                 = '/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/Demographic_data_FCD.txt';
Left_Right                = 'both';
NumMesh                   = 81920;
KERNEL                    = 5;
OUTPATH = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/lesion_detection_frontiers/98_TLE/lesion_detection/';

fid = fopen(Cases_TLE);
demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
case_num_TLE = demo{1};
age_TLE = demo{2};
gender_TLE = demo{3};
fclose(fid);

fid = fopen(Cases_pat);
demo = textscan(fid, '%s%f%s%d%s%d%d', 'Delimiter', ',', 'CollectOutput', 1);
case_num_pat = demo{1};
age_pat = demo{2};
gender_pat = demo{3};
lesion_volume = demo{4};
histo_type = demo{5};
initial = demo{6}(:, 1);
transmantle = demo{6}(:, 2);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. variable initialization - Include data that will be analyzed and 
%%                              Compute age and gender distribution
load([ 'db_FCD_3T_k5.mat' ]);
load('db_lesion_3T.mat');
z_score_database_pat = z_score_database;

load([ 'db_TLE_3T_k5.mat' ]);
z_score_database_TLE = z_score_database;
clear('z_score_database');

mean_age_TLE = mean(age_TLE);
std_age_TLE = std(age_TLE, 0);

mean_age_pat = mean(age_pat);
std_age_pat = std(age_pat, 0);

[h,p,ci,stats] = ttest2(age_TLE,age_pat);

pat_case = { '027', '049', '061', '062', '063', '066', '071', '072', '073', '074', '076', '079', '080_1', '081', '082', '083', '084', '085', '086' };
TLE_case = {  '0361_1', '0363_1', '0369_1', '0375_1', '0380_1', '0390_1', '0392_1', '0394_1', '0395_1', '0396_1', '0404_1', };          

[C, pat_idx, ib]  = intersect(case_num_pat, pat_case);
[C, TLE_idx, ib] = intersect(case_num_TLE, TLE_case);

mean_age_TLE_new = mean(age_TLE(TLE_idx));
std_age_TLE_new = std(age_TLE(TLE_idx), 0);

mean_age_pat_new = mean(age_pat(pat_idx));
std_age_pat_new = std(age_pat(pat_idx), 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. variable initialization - Read AAL map to exclude the FP cluster in the medial brain surface
AAL_data = SurfStatReadData('aal_both_rsl_final.txt');
% figure; SurfStatView(AAL_data, TemplateSurf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. variable initialization - Case exclusion
pat_outlier = [  ];               %% Patient cases that will be excluded in this analysis
pat_case(pat_outlier)
pat_case_temp = pat_case;
pat_case_temp(pat_outlier) =[];

TLE_outlier = [  ];              %% TLE cases that will be excluded in this analysis
TLE_case(TLE_outlier)
TLE_case_temp = TLE_case;
TLE_case_temp(TLE_outlier) =[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. variable initialization - Parameters
smooth_kernel = 5;          %% lesion blurring
binarize_thres = 0.005;     %% lesion blurring threshold
FeatSubSet = [1 2 3 4 5];   %% 1: RI, 2: PG, 3: CT, 4: SD, 5: CV
Siz_thres = 0;              
sampling_scale = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. variable initialization - Database reshape
TemplateSurf = SurfStatReadSurf({'/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/surf_reg_model_left.obj', ...
                                 '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/surf_reg_model_right.obj'});
MaskCut      = SurfStatMaskCut(TemplateSurf);
lesionMasks  = permute(lesion_db, [2 1]);
lesionMask2  = SurfStatSmooth(lesionMasks', TemplateSurf, smooth_kernel);
lesionMasks  = double(lesionMask2>binarize_thres)';

featVect = permute(z_score_database_pat, [2 3 1]);
featVect(:, :, pat_outlier) = [];
lesionMasks(:, pat_outlier) = [];
lesionMasks(MaskCut<1,:)=2;

%% initialization
siz = size(featVect);
numVert = siz(1);
numFeat = siz(2);
numPatient = siz(3);

HTissueMask = lesionMasks;
HTissueMask(:,:) = 0;
sizLesion = sum(sum(lesionMasks == 1));
sizHTissueperSubj = round(sizLesion / numPatient) * sampling_scale;
Candidates = findn(lesionMasks == 0);

%% Random healthy tissue sampling while balancing its sample size with
%% lesional tissue
for k = 1:numPatient
    sizeFullHTissue = sum(Candidates(:,2)==k);
    CandEachPatient =Candidates(Candidates(:,2)==k,1);
    SampleIndex = rand(sizeFullHTissue,1) > 1 - sizHTissueperSubj / sizeFullHTissue;
    sum(SampleIndex~=0)
    HTissueMask(CandEachPatient(SampleIndex),k) = 1;
end

%% TLE database reshape
featVect_TLE = permute(z_score_database_TLE, [2 3 1]);
featVect_TLE(:, :, TLE_outlier) = [];

%% initialization
siz = size(featVect_TLE);
numTLE = siz(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. First classification - Train the LDA classifier using patients' lesional and nonlesional tissue information and
%%                           Test it for TLE data
PostriorAll = zeros(numVert, numTLE);
for j = 1:numTLE
    % exclude the test data
    clear dataLS dataHT
    dataLS(:,:) = featVect(logical(lesionMasks(:,1)==1), :, 1);
    for k = 2:numPatient
        dataLS = cat(1, dataLS, featVect(logical(lesionMasks(:,k)==1), :, k));
    end
    dataHT(:,:) = featVect(logical(HTissueMask(:,1)), :, 1);
    for k = 2:numPatient
        dataHT = cat(1, dataHT, featVect(logical(HTissueMask(:,k)), :, k));
    end
    
    % Classify the test data
    TrainData = cat(1, dataHT, dataLS);
    TestData = featVect_TLE(:, :, j);
    TrainClass = cat(2, ones(1, length(dataHT)), ones(1, length(dataLS))*2)';
    
    [ldaClass, error,POSTERIOR, logp, coeff]  = classify(TestData(:,FeatSubSet), TrainData(:,FeatSubSet), TrainClass, 'diagquadratic');
    PostriorAll(:, j) = POSTERIOR(:,2);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. First classification: compute sensitivity and false positivie rate
interval = 0.005;
range = 0.5 : interval : 1;
specificity_first = zeros(size(range, 2), numTLE);
false_positive_first = zeros(size(range, 2), numTLE);
for post_thres = range
    for j = 1:numTLE
        idx = int32((post_thres-0.5)/interval+1);
        FP = sum(PostriorAll(:, j) >= post_thres & AAL_data' ~= 0);
        TN = sum(PostriorAll(:, j) < post_thres & AAL_data' ~= 0);              
        false_positive_first(idx,j) = FP;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. Second classification initialization
post_thres = 0.93; %% same threshold used in the FCD 3T lesion detection.

idx = find(abs(range - post_thres)<0.000001);
% False positive ratio of the vertex-wise classifier
false_positive_first_temp = false_positive_first(idx, :)/81924

%% visualize
% for j = 1 : numTLE
%     f = figure; SurfStatView((PostriorAll(:, j)'>post_thres).*(AAL_data~=0), TemplateSurf, TLE_case_temp{j});
%     exportfigbo(f,[OUTPATH '/' TLE_case_temp{j} '_FP1.png' ], 'png', 6); close(f);
% end

surf_data = surfGetNeighborsHong(TemplateSurf);
edg=SurfStatEdg(TemplateSurf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10. Second classification to detect FP
%%  -   Feature sets
%%      1: RI, 2: PG, 3: CT, 4: SD, 5: CV
%%             Size     mean            std           skewness          Krutosis            Moment           Asymetry      AAL parcel  x,y,z central coord 
%%               1    2 3 4 5 6     7 8 9 10 11    12 13 14 15 16    17 18 19 20 21     22 23 24 25 26    27 28 29 30 31        32           33 34 35
Subset_temp = { [1]  [2 3 4 5 6]   [7 8 9 10 11]  [12 13 14 15 16]  [17 18 19 20 21]   [22 23 24 25 26]  [27 28 29 30 31]      [32]         [33 34 35] };

%%  -   Feature computation
lesionClass_TLE = PostriorAll > post_thres;
MaxNumClus = 300;
bLesion_Pt_TLE = zeros(MaxNumClus, numTLE);
ClusterFeat_Pt_TLE = zeros(MaxNumClus, size(cell2mat(Subset_temp), 2), numTLE);
featVect_TLE_db = [];
Total_clus_num = zeros(1, numTLE);
Size_FP_clus_num = zeros(1, numTLE);
FP_clus_num = zeros(MaxNumClus, numTLE);
Clusters_TLE = double(lesionClass_TLE);
order = 6;

%% Detect the clusters and remove FPs according to following criteria:
%% 1) FP cluster which doesn't meet the condition of 6-connected cluster 
%%    This is for removing the peculiar shape of clusters, for instance, zigzaged elongated cluster.
%% 2) FP cluster which part is in the mid mask area

%% Collect feature distribution of clusters after prunning process separately in lesion and FP clusters
for j = 1 : numTLE
    l=1;
    
    Clusters_TLE(:,j) = SurfStatCluster(lesionClass_TLE(:,j)', TemplateSurf);
    detected_lesion_size = [ 0 0 ];    
    Total_clus_num(j) = max(Clusters_TLE(:,j));
    
    for k = 1 : Total_clus_num(j)
        Size = sum(Clusters_TLE(:,j)==k);
        flag = 0;
        vert_idx = find(Clusters_TLE(:, j) == k);
        
        % FP cluster which doesn't meet the condition of 6-connected cluster 
        for v = 1 : size(vert_idx, 1)
            nb_v = surf_data.nbr(vert_idx(v), :);
            nb_v = sort(nb_v(nb_v ~= 0));
            if(isequal(nb_v, intersect(nb_v, vert_idx)))
                flag = 1;
                break;
            end
        end
        
        if flag == 1
            % FP cluster which is smaller than size threshold 
            if Size > Siz_thres
                % FP cluster which part is in the mid mask area
                if(sum(Clusters_TLE(:,j)==k & AAL_data' == 0))
                    Clusters_TLE(Clusters_TLE(:,j)==k, j) = 0;
                    Size_FP_clus_num(j) = Size_FP_clus_num(j) + 1;
                    continue;
                else
                    bLesion_Pt_TLE(l,j) = 500 + k;
                    FP_clus_num(l,j) = FP_clus_num(l,j) + Size;
                end

                mahal_dist = sum(featVect_TLE(vert_idx, :, j).^2, 2);
                [maxPost, id] = max(PostriorAll(vert_idx, j));
                id = find(PostriorAll(vert_idx, j) == maxPost);
                
                [maxMahal, id_maxMahal] = max(mahal_dist(id));
                temp = vert_idx(id(id_maxMahal));
                peakCluster = [ temp; find_vertex_multineighbors_by_order(TemplateSurf,temp,order,edg) ]';
                                
                M1 = mean(featVect_TLE(peakCluster,:,j), 1);
                M2 = std(featVect_TLE(peakCluster,:,j), 1);
                M3 = skewness(featVect_TLE(peakCluster,:,j), 1);
                M4 = kurtosis(featVect_TLE(peakCluster,:,j), 1);
                M5 = moment(featVect_TLE(peakCluster,:,j), 5);

                Cx = TemplateSurf.coord(1, peakCluster);
                Cy = TemplateSurf.coord(2, peakCluster);
                Cz = TemplateSurf.coord(3, peakCluster);
                
                TempDist = dist3(mean([Cx; Cy; Cz],2)', [Cx; Cy; Cz]');
                nn = find(TempDist == min(TempDist));
                finalModCoord = [Cx(nn(1)); Cy(nn(1)); Cz(nn(1))];
                vertNN = find(TemplateSurf.coord(1,:)==finalModCoord(1) & TemplateSurf.coord(2,:)==finalModCoord(2)  & TemplateSurf.coord(3,:)==finalModCoord(3));
                AAL_vertNN = AAL_data(vertNN);
                TTT=find(Clusters_TLE(:,j)==k);
                TTT2 = 40962+TTT;
                if sum(TTT2 <= 81924 & TTT2 > 40962)
                    Asym = mean(featVect_TLE(TTT,:,j)) - mean(featVect_TLE(TTT2,:,j));
                else
                    TTT2 = TTT - 40962;
                    Asym = mean(featVect_TLE(TTT,:,j)) - mean(featVect_TLE(TTT2,:,j));
                end
                ClusterFeat_Pt_TLE(l, :, j) = cat(2, Size, M1, M2, M3, M4, M5, Asym, AAL_vertNN, finalModCoord');
                l = l+1;                
            else
                Size_FP_clus_num(j) = Size_FP_clus_num(j) + 1;
            end
        else
            Clusters_TLE(vert_idx, :) = 0;            
            Size_FP_clus_num(j) = Size_FP_clus_num(j) + 1;
        end
    end
end

%% We used the classifier trained by 3T data of FCD patients to test controls
load(['/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/lesion_detection/db_FCD_3T_k5_for_second_classifier.mat']);


%%  -   Training a classifier with different combination of feature sets and
%%      Test it for TLE data ...
size_temp = 0;
for c = 1 : size(Subset_temp, 2)
    cbnk = combnk(1:size(Subset_temp, 2), c);
    size_temp = size_temp + size(cbnk, 1);
end

PostriorAll_second = zeros(MaxNumClus, size_temp, numTLE);
lesionMask_2 = bLesion_Pt;
featVect_2 = ClusterFeat_Pt;
count = 1;
Subset_temp_idx = cell2mat(Subset_temp);
for c = 1 : size(Subset_temp, 2)
    cbnk = combnk(1:size(Subset_temp, 2), c);
    for cb = 1 : size(cbnk, 1)
        if(mod(count,10) == 1)
            disp(['Training done percent: ' num2str(round(count/size_temp*100),2) '%'])
        end
        Subset = cell2mat(Subset_temp(cbnk(cb, :)));
        [c, ia, ib] = intersect(Subset, Subset_temp_idx);
        Subset = ib;
        for j = 1:numTLE
            clear dataLS
            dataLS(:,:) = featVect_2(logical(lesionMask_2(:,1)>0), :, 1);
            for k = 2:numPatient
                dataLS = cat(1, dataLS, featVect_2(logical(lesionMask_2(:,k)>0), :, k));
            end

            % Classify the test data
            TrainData = dataLS;
            TestData =ClusterFeat_Pt_TLE(logical(bLesion_Pt_TLE(:,j)>0), :, j);
            TrainClass = lesionMask_2(logical(lesionMask_2>0));
            TrainClass(TrainClass>500) = 1;
            TrainClass(TrainClass>200 & TrainClass<500) = 2;

            [ldaClass, error, POSTERIOR, logp, coeff ]  = classify(TestData(:, Subset), TrainData(:, Subset), TrainClass, 'linear');
            PostriorAll_second(1:length(POSTERIOR(:,2)), count, j) = POSTERIOR(:,2);            
        end
        count = count + 1;
    end
end

PostriorAll_second_temp = PostriorAll_second;

%%  -   Compute sensitivity and false positivie rate
interval = 0.01;
range2 = 0.0 : interval : 1;

rank_performance = zeros(size_temp, 2, size(range2, 2));
sensitivity_set = zeros(size_temp, numTLE, size(range2, 2));
specificity_set = zeros(size_temp, numTLE, size(range2, 2));
specific_area_set = zeros(size_temp, numTLE, size(range2, 2));
removed_FP_rate_set = zeros(size_temp, numTLE, size(range2, 2));
Subset_set = {};
combnk_set = [];
cluster_set = cell(size_temp, 2, numTLE, size(range2, 2));

for p = 1 : size(range2, 2)
    post_thres_temp = range2(p);

    sensitivity = zeros(1, numTLE);
    specificity = zeros(1, numTLE);
    specific_area = zeros(1, numTLE);
    removed_FP_rate = zeros(1, numTLE);
        
    count = 1;
    for c = 1 : size(Subset_temp, 2)
        cbnk = combnk(1:size(Subset_temp, 2), c);        
        for cb = 1 : size(cbnk, 1)
            Subset = cell2mat(Subset_temp(cbnk(cb, :)));
            for j = 1:numTLE
                FP = sum(bLesion_Pt_TLE(:, j)>500 & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                TN = sum(bLesion_Pt_TLE(:, j)>500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);

                FP_c = find(bLesion_Pt_TLE(:, j)>500 & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                TN_c = find(bLesion_Pt_TLE(:, j)>500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);

                cluster_set(count, :, j, p) = [ {FP_c'}, {TN_c'} ];
                specificity(j) = TN / ( TN + FP );
                specific_area(j) = sum(ClusterFeat_Pt_TLE(bLesion_Pt_TLE(:, j)>500 & PostriorAll_second_temp(:, j) < post_thres_temp, 1, j)) / sum(ClusterFeat_Pt_TLE(bLesion_Pt_TLE(:, j)>500,1,j));
                removed_FP_rate(j) = (TN+Size_FP_clus_num(j))/Total_clus_num(j);
            end
            
            Subset_set{count} =Subset;
            combnk_set(count, :) = [ c, cb ];           
            specificity_set(count, :, p) = specificity;
            specific_area_set(count, :, p) = specific_area;            
            removed_FP_rate_set(count, :, p) = removed_FP_rate;
            count = count + 1
        end        
    end
end

temp_performance = [ squeeze(mean(removed_FP_rate_set, 2)) ];

for p = 1 : size(range2, 2)    
    [b, i] = sort(temp_performance(:, p));
    rank_performance(:, :, p) = [ i, temp_performance(i,  p)];
end

performance = zeros(size(range2, 2), 2);
for p = 1 : size(range2, 2)
    temp = rank_performance(:, 2, p);
    [b, i] = max(temp);
    performance(p, :) = rank_performance(i, :, p);
end

[range2' performance]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters ...
specificity_set(isnan(specificity_set)) = 1;
post_thres_temp = 0.90; %% exactly same threshold used in detection of FCD patients
m = 165;                %% exactly same combination used in detection of FCD patients
p = find(abs(range2-post_thres_temp)<0.00001);
Subset = Subset_set{m}
PostriorAll_second2 = zeros(MaxNumClus, numPatient);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Result 
FP_set = [];
TN_set = [];
FP_quant = [];

FP_set_inv = [];
second_FP_set = [];
second_TN_set = [];
second_FP_quant = [];

Total_clus = [];

for j = 1 : numTLE  
    % The first classification ...
    FP_quant_temp = 0;
    FP_set_temp = [];
    second_FP_set_temp = [];
    second_TN_set_temp = [];

    Total_clus = [ Total_clus max(Clusters_TLE(:, j)) ];
    for k = 1 : max(Clusters_TLE(:, j))
        Size = sum(Clusters_TLE(:, j)==k);
        flag = 0;
        vert_idx = find(Clusters_TLE(:, j) == k);
        for v = 1 : size(vert_idx, 1)
            nb_v = surf_data.nbr(vert_idx(v), :);
            nb_v = sort(nb_v(nb_v ~= 0));
            if(isequal(nb_v, intersect(nb_v, vert_idx)))
                flag = 1;
                break;
            end
        end
        
        cluster_temp = Clusters_TLE(:, j) == k;
        
        if flag == 1
            if(Size > Siz_thres)
                FP_quant_temp = FP_quant_temp + 1;
                FP_set_temp = [ FP_set_temp sum(cluster_temp) ];
                FP_set_inv  = [ FP_set_inv sum(cluster_temp) ];
            else
                if(Size > 0)
                    FP_quant_temp = FP_quant_temp + 1;
                    FP_set_inv   = [ FP_set_inv sum(cluster_temp) ];
                end
            end
        end
    end

    FP_quant = [ FP_quant mean(FP_quant_temp) ];
    FP_set   = [ FP_set mean(FP_set_temp) ];
    TN_set   = [ TN_set sum(Clusters_TLE(:, j) == 0) ];
    
    %% The second classficiation to reduce FP ...
    clusid = bLesion_Pt_TLE(cluster_set{m, 1, j, p}, j) - 500;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 1)
            second_FP_set_temp = [ second_FP_set_temp sum(Clusters_TLE(:, j) == clusid(k)) ];
        end
    end
    
    clusid = bLesion_Pt_TLE(cluster_set{m, 2, j, p}, j) - 500;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 1)
            second_TN_set_temp = [ second_TN_set_temp sum(Clusters_TLE(:, j) == clusid(k)) ];
        end
    end
   
    second_FP_set = [ second_FP_set mean(second_FP_set_temp) ];
    second_TN_set = [ second_TN_set mean(second_TN_set_temp) ];
    second_FP_quant = [ second_FP_quant size(cluster_set{m, 1, j, p}, 2) ];
end

%% False positive mean size and standard deviation
FP_set(isnan(FP_set)) = 0;
mean(FP_set) 
std(FP_set)

sum(FP_set_inv)

%% False positive mean amount and standard deviation
mean(FP_quant) 
std(FP_quant)

%% Second False positive mean size and standard deviation
second_FP_set(isnan(second_FP_set)) = 0;
mean(second_FP_set) 
std(second_FP_set)

%% Second Classification False positive mean amount and standard deviation
second_FP_quant(isnan(second_FP_quant)) = 0;
mean(second_FP_quant) 
std(second_FP_quant)

%% Visualization and report final results of classification
OUTPATH = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/lesion_detection_frontiers/98_TLE/lesion_detection/';
cluster_info = cell(numTLE, 2);
count = 0;
for j = 1 : numTLE
%     f = figure; SurfStatView(Clusters_TLE(:, j).*MaskCut', TemplateSurf, [ TLE_case{j} ' fp:' num2str(false_positive_first(j)*100/81924) '%' ]);
%     exportfigbo(f,[OUTPATH '/' TLE_case{j} '_first_classification.png' ], 'png', 6); close(f);
    
    cluster_temp = zeros(1, 81924);
    
    % False postives
    clusid = bLesion_Pt_TLE(cluster_set{m, 1, j, p}, j) - 500;
    if(~isempty(clusid))
        for t = 1 : size(clusid, 2)
            cluster_temp(Clusters_TLE(:, j) == clusid(t)) = 1;
            cluster_info{j, 2} =  [ cluster_info{j, 2} sum(Clusters_TLE(:, j) == (bLesion_Pt_TLE(cluster_set{m, 1, j, p}, j) - 500)) ];
        end
        count = count + sum(cluster_info{j, 2});
        cluster_info{j, 2}
    end

%     f = figure; SurfStatView(cluster_temp.*MaskCut, TemplateSurf, [ TLE_case{j} ', specificity: ' num2str(specificity_set(m, j)) '%' ]);
%     SurfStatColLim([0 5]);
%     exportfigbo(f,[OUTPATH '/' TLE_case{j} '_second_classification.png' ], 'png', 6); close(f);
end

cluster_info = [ TLE_case_temp' cluster_info ]
