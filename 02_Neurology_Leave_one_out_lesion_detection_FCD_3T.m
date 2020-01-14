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
%  3) Variable description
%     TemplateSurf: template surface that defines the number of vertices & triangles
%     featVect: feature Vectors, l x m x n matrix, size of l: #vertices, m: #feature
%     types (i.e., thickness, curvature, RI, gradient, etc), n: #FCD patients
%     lesionMasks: lesional masks for all patients, l x n: size of l:#vertices,
%     Please assign FCD lesional vertices as 1 otherwise 0 (i.e., healthy tissues).
%     If you like to exclude any verices from sampling, assign them 2.
clear all;
addpath('/local_raid/seokjun/01_project/01_Surface_based_lesion_detection/01_analysis/utility');

%% Results used in the paper. For fast running of the script, you can always re-load all the result to reduce the time of all the re-computation 
load('/local_raid/seokjun/01_project/01_Surface_based_lesion_detection/01_analysis/Result_FCD_3T.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. variable initialization - read files
Group_cont                = 'control';
Prefix_cont               = 'TLE';
Cases_cont                = '/local_raid/seokjun/01_project/01_Surface_based_lesion_detection/01_analysis/Demographic_data_control_3T.txt';
Group_pat                 = 'FCD';
Prefix_pat                = 'mcd';
Cases_pat                 = '/local_raid/seokjun/01_project/01_Surface_based_lesion_detection/01_analysis/Demographic_data_FCD_3T.txt';
Left_Right                = 'both';
NumMesh                   = 81920;
KERNEL                    = 5;
OUTPATH                   = '/local_raid/seokjun/01_project/01_Surface_based_lesion_detection/01_analysis/02_result/';

fid = fopen(Cases_cont);
demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
case_num_cont = demo{1};
age_cont = demo{2};
gender_cont = demo{3};
fclose(fid);

fid = fopen(Cases_pat);
demo = textscan(fid, '%s%f%s%d%s%d%d%s%s', 'Delimiter', ',', 'CollectOutput', 1);
case_num_pat = demo{1};
age_pat = demo{2};
gender_pat = demo{3};
lesion_volume = demo{4};
histo_type = demo{5};
initial = demo{6}(:, 1);
transmantle = demo{6}(:, 2);
lesoin_location = demo{7}(:, 1);
lesoin_laterlization = demo{7}(:, 2);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. variable initialization - Include data that will be analyzed and 
%%                              Compute age and gender distribution
load([ 'db_FCD_3T_k5.mat' ]);
load('db_lesion_3T.mat');

mean_age_cont = mean(age_cont);
std_age_cont = std(age_cont, 0);

mean_age_pat = mean(age_pat);
std_age_pat = std(age_pat, 0);

[h,p,ci,stats] = ttest2(age_cont,age_pat)

pat_case  = { '027', '049', '061', '062', '063', '066', '071', '072', '073', '074', '076', '079', '080_1', '081', '082', '083', '084', '085', '086' };
cont_case = { '301_1', '303_1', '304_1', '305_1', '307_1', '308_1', '309_1', '310_1', '311_1', '314_1', '316_1', '317_1', '321_1', '323_1', '329_1', '332_1', '324_1', '333_1', '334_1', '335_1', '337_1', '338_1', '339_1', '340_1'  };

[C, pat_idx, ib]  = intersect(case_num_pat, pat_case);
[C, cont_idx, ib] = intersect(case_num_cont, cont_case);

mean_age_cont_new = mean(age_cont(cont_idx));
std_age_cont_new = std(age_cont(cont_idx), 0);
sum(strcmp(gender_cont(cont_idx), 'm'))

mean_age_pat_new = mean(age_pat(pat_idx));
std_age_pat_new = std(age_pat(pat_idx), 0);
sum(strcmp(gender_pat(pat_idx), 'm'))
histo_type_temp = histo_type(pat_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. variable initialization - Read AAL map to exclude the FP cluster in the medial brain surface
AAL_data = SurfStatReadData('aal_both_rsl_final.txt');
% figure; SurfStatView(AAL_data, TemplateSurf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. variable initialization - Case exclusion
outlier = [  ];              %% Cases that will be excluded in this analysis
pat_case(outlier)
pat_case_temp = pat_case;
pat_case_temp(outlier) =[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. variable initialization - Parameters
smooth_kernel = 5;          %% lesion blurring
binarize_thres = 0.005;     %% lesion blurring threshold
FeatSubSet = [1 2 3 4 5];   %% 1: RI, 2: PG, 3: CT, 4: SD, 5: CV
Siz_thres = 0;             
sampling_scale = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. variable initialization - Database reshape
TemplateSurf = SurfStatReadSurf({'/local_raid/seokjun/01_project/01_Surface_based_lesion_detection/01_analysis/surf_reg_model_left.obj', ...
                                 '/local_raid/seokjun/01_project/01_Surface_based_lesion_detection/01_analysis/surf_reg_model_right.obj'});
MaskCut      = SurfStatMaskCut(TemplateSurf);
lesionMasks  = permute(lesion_db, [2 1]);
lesional_volume = [ 5400 217 1131 6167 1087 361 8119 223 847 45 585 2545 1366 8672 417 516 264 2679 3305];
lesionMask2  = SurfStatSmooth(lesionMasks', TemplateSurf, smooth_kernel);
lesionMasks  = double(lesionMask2>binarize_thres)';

lesion_size = zeros(1, size(pat_case_temp, 2));
for i = 1 : size(pat_case_temp, 2)
    lesion_size(i) = sum(lesionMasks(:, i) == 1);
end

featVect = permute(z_score_database, [2 3 1]);
featVect(:, :, outlier) = [];
lesionMasks(:, outlier) = [];
lesionMasks(MaskCut<1,:)=2;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. First classification: Leave-one-out
PostriorAll = zeros(numVert, numPatient);
for j = 1:numPatient
    % exclude the test data
    lesionMask_leaveone = lesionMasks;
    lesionMask_leaveone(:, j) = [];
    HTissueMask_leaveone = HTissueMask;
    HTissueMask_leaveone(:, j) = [];
    featVect_leaveone = featVect;
    featVect_leaveone(:, :, j) = [];
    clear dataLS dataHT
    dataLS(:,:) = featVect_leaveone(logical(lesionMask_leaveone(:,1)==1), :, 1);
    for k = 2:(numPatient-1)
        dataLS = cat(1, dataLS, featVect_leaveone(logical(lesionMask_leaveone(:,k)==1), :, k));
    end
    dataHT(:,:) = featVect_leaveone(logical(HTissueMask_leaveone(:,1)), :, 1);
    for k = 2:(numPatient-1)
        dataHT = cat(1, dataHT, featVect_leaveone(logical(HTissueMask_leaveone(:,k)), :, k));
    end
    
    % Classify the test data
    TrainData = cat(1, dataHT, dataLS);
    TestData = featVect(:, :, j);
    TrainClass = cat(2, ones(1, length(dataHT)), ones(1, length(dataLS))*2)';
    
    [ldaClass, error,POSTERIOR, logp, coeff ]  = classify(TestData(:,FeatSubSet), TrainData(:,FeatSubSet), TrainClass, 'diagquadratic');
    PostriorAll(:, j) = POSTERIOR(:,2);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. First classification: compute sensitivity and false positivie rate
interval = 0.005;
range = 0.5 : interval : 1;
sensitivity_first = zeros(size(range, 2), numPatient);
specificity_first = zeros(size(range, 2), numPatient);
Coverage_first = zeros(size(range, 2), numPatient);
Coverage = zeros(size(range, 2), numPatient);
false_positive_first = zeros(size(range, 2), numPatient);
false_positive_first_rate = zeros(size(range, 2), numPatient);
for post_thres = range
    for j = 1:numPatient
        idx = int32((post_thres-0.5)/interval+1);
        TP = sum(lesionMasks(:, j)==1 & PostriorAll(:, j) >= post_thres & AAL_data' ~= 0);
        FP = sum(lesionMasks(:, j)==0 & PostriorAll(:, j) >= post_thres & AAL_data' ~= 0);
        TN = sum(lesionMasks(:, j)==0 & PostriorAll(:, j) < post_thres & AAL_data' ~= 0);
        FN = sum(lesionMasks(:, j)==1 & PostriorAll(:, j) < post_thres & AAL_data' ~= 0);
        
        %% They are not conventional defintion of sensitivity and specificity 
        sensitivity_first(idx,j) = TP > 0;
        specificity_first(idx,j) = TN / sum(lesionMasks(:, j)==0);
        Coverage_first(idx,j) = TP / sum(lesionMasks(:,j)==1);
        Coverage(idx, j) = TP;
        false_positive_first(idx,j) = FP;
        false_positive_first_rate(idx, j) = FP / (FP + TN);
    end
end
[a b] = max(mean(sensitivity_first, 2) - mean(false_positive_first, 2)/81924);

% To mimic the figure 2 of the Besson's paper ISBI 2008
% Visualization of the vertex-wise classifier performance
interval = 0.005;
range = 0.5 : interval : 1;
range_temp = int32(((0.5 - range(1))/interval + 1) : ((1 - range(1))/interval + 1));
figure; 
[ AX, H1, H2 ] = plotyy(range(range_temp), mean(sensitivity_first(range_temp, :), 2)/mean(sensitivity_first(b, :)), ...
                        range(range_temp), mean(false_positive_first(range_temp, :), 2)/mean(false_positive_first(b, :)), 'plot');
set(get(AX(1),'Ylabel'),'String','Ratio Detection rate'); 
set(get(AX(2),'Ylabel'),'String','Ratio FP'); 
xlabel('Threshold for posterior');
title('Tradeoff between detection rate and FP');
set(H1,'LineStyle','--');
set(H2,'LineStyle','--');
bottom = 0; upper = 2; interval_axis = 0.5;
set(AX(1), 'YTick', [ bottom : interval_axis : upper ]);
set(AX(1), 'YLim',  [ bottom upper ]);
set(AX(1), 'YTick', [ bottom : interval_axis : upper ]);
set(AX(2), 'YLim',  [ bottom upper ]);
set(AX(2), 'YTick', [ bottom: interval_axis : upper ]);
grid on; hold on; scatter(range(b), 1.1, 'rv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. Second classification initialization
idx = find(range == range(b));

% Missing cases
pat_case_temp((sensitivity_first(idx, :)==0))

% Sensitivity and false positive ratio of the vertex-wise classifier
sensitivity_first_temp = sensitivity_first(idx, :)
false_positive_first_temp = false_positive_first(idx, :)/89124
post_thres = range(b);

surf_data = surfGetNeighborsHong(TemplateSurf);
edg=SurfStatEdg(TemplateSurf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10. Second classification to remove FP
%%  -   Feature sets
%%      1: RI, 2: PG, 3: CT, 4: SD, 5: CV
%%             Size     mean            std           skewness          Krutosis            Moment           Asymetry      AAL parcel  x,y,z central coord 
%%               1    2 3 4 5 6     7 8 9 10 11    12 13 14 15 16    17 18 19 20 21     22 23 24 25 26    27 28 29 30 31        32           33 34 35
Subset_temp = { [1]  [2 3 4 5 6]   [7 8 9 10 11]  [12 13 14 15 16]  [17 18 19 20 21]   [22 23 24 25 26]  [27 28 29 30 31]      [32]         [33 34 35] };

%%  -   Feature computation
lesionClass = PostriorAll >= post_thres;
MaxNumClus = 300;
bLesion_Pt = zeros(MaxNumClus, numPatient);
ClusterFeat_Pt = zeros(MaxNumClus, size(cell2mat(Subset_temp), 2), numPatient);
featVect_db = [];
Total_clus_num = zeros(1, numPatient);
Size_FP_clus_num = zeros(1, numPatient);
FP_clus_num = zeros(MaxNumClus, numPatient);
Clusters = zeros(81924, numPatient);
order=6;

count_temp = [];
OUTPATH = '/local_raid/seokjun/01_project/IntracorticalAnalysis/03_result/lesion_detection_frontiers/14_result/';

%% Detect the clusters and remove FPs according to following criteria:
%% 1) FP cluster which doesn't meet the condition of 6-connected cluster 
%%    This is for removing the peculiar shape of clusters, for instance, zigzaged elongated cluster.
%% 2) FP cluster which part is in the mid mask area

%% Collect feature distribution of clusters after prunning process separately in lesion and FP clusters
for j = 1 : numPatient
    j
    l=1;
    
    Clusters(:,j) = SurfStatCluster(lesionClass(:,j)', TemplateSurf);        
    Total_clus_num(j) = max(Clusters(:,j));
    
    for k = 1 : Total_clus_num(j)
        Size = sum(Clusters(:,j)==k);
        flag = 0;
        vert_idx = find(Clusters(:, j) == k);
        
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
            if(Size > Siz_thres)
                if sum(lesionMasks(:, j)==1 & Clusters(:,j)==k)
                        bLesion_Pt(l,j) = 200 + k;                                                                                          
                else
                    % FP cluster which part is in the mid mask area
                    if(sum(Clusters(:,j)==k & AAL_data' == 0))
                        Clusters(Clusters(:,j)==k, j) = 0;                        
                        Size_FP_clus_num(j) = Size_FP_clus_num(j) + 1;
                        continue;
                    else
                        bLesion_Pt(l,j) = 500 + k;
                        FP_clus_num(l,j) = FP_clus_num(l,j) + Size;                        
                    end
                end             
                
                mahal_dist = sum(featVect(vert_idx, :, j).^2, 2);
                [maxPost, id] = max(PostriorAll(vert_idx, j));
                id = find(PostriorAll(vert_idx, j) == maxPost);
                
                [maxMahal, id_maxMahal] = max(mahal_dist(id));
                temp = vert_idx(id(id_maxMahal));
                peakCluster = [ temp; find_vertex_multineighbors_by_order(TemplateSurf,temp,order,edg) ]';

                M1 = mean(featVect(peakCluster,:,j), 1);
                M2 = std(featVect(peakCluster,:,j), 1);
                M3 = skewness(featVect(peakCluster,:,j), 1);
                M4 = kurtosis(featVect(peakCluster,:,j), 1);
                M5 = moment(featVect(peakCluster,:,j), 5);

                Cx = TemplateSurf.coord(1, peakCluster);
                Cy = TemplateSurf.coord(2, peakCluster);
                Cz = TemplateSurf.coord(3, peakCluster);

                TempDist = dist3(mean([Cx; Cy; Cz],2)', [Cx; Cy; Cz]');
                nn = find(TempDist == min(TempDist));
                finalModCoord = [Cx(nn(1)); Cy(nn(1)); Cz(nn(1))];
                vertNN = find(TemplateSurf.coord(1,:)==finalModCoord(1) & TemplateSurf.coord(2,:)==finalModCoord(2)  & TemplateSurf.coord(3,:)==finalModCoord(3));
                AAL_vertNN = AAL_data(vertNN);
                TTT=find(Clusters(:,j)==k);
                TTT2 = 40962+TTT;
                if sum(TTT2 <= 81924 & TTT2 > 40962)
                    Asym = mean(featVect(TTT,:,j)) - mean(featVect(TTT2,:,j));
                else
                    TTT2 = TTT - 40962;
                    Asym = mean(featVect(TTT,:,j)) - mean(featVect(TTT2,:,j));
                end            
                ClusterFeat_Pt(l, :, j) = cat(2, Size, M1, M2, M3, M4, M5, Asym, AAL_vertNN, finalModCoord');
                l = l+1;
            else
                Clusters(vert_idx, j) = 0;
                Size_FP_clus_num(j) = Size_FP_clus_num(j) + 1;
            end
        else
            Clusters(vert_idx, j) = 0;            
            Size_FP_clus_num(j) = Size_FP_clus_num(j) + 1;
        end
    end
end

%% Estimate the feature pattern in TP and FP clusters
% CluterFeat_Pt_TP = [];
% CluterFeat_Pt_FP = [];
% for j = 1 : numPatient
%     for k = 1 : size(bLesion_Pt(:, j), 1)
%         if(bLesion_Pt(k, j) < 300 & bLesion_Pt(k, j) > 200)
%             CluterFeat_Pt_TP = [ CluterFeat_Pt_TP; ClusterFeat_Pt(k, :, j) ];
%         elseif(bLesion_Pt(k, j) > 500)
%             CluterFeat_Pt_FP = [ CluterFeat_Pt_FP; ClusterFeat_Pt(k, :, j) ];
%         end
%     end
% end
% 
% feature_evaluation = zeros(1, size(CluterFeat_Pt_TP, 2));
% for i = 1 : size(CluterFeat_Pt_TP, 2)
%     [h,p,ci,stats] = ttest2(CluterFeat_Pt_TP(:, i), CluterFeat_Pt_FP(:, i));
%     feature_evaluation(i) = p;    
% end
% 
% i = 9; figure; scatter(ones(size(CluterFeat_Pt_TP(:, i), 1), 1), CluterFeat_Pt_TP(:, i), 'bo');
% hold on;       scatter(ones(size(CluterFeat_Pt_FP(:, i), 1), 1)*2, CluterFeat_Pt_FP(:, i), 'ro');
% hold on;       scatter(1.2, mean(CluterFeat_Pt_TP(:, i)), 'b*');
% hold on;       scatter(2.2, mean(CluterFeat_Pt_FP(:, i)), 'r*');
% xlim([0 3]);
% 
% [ mean(CluterFeat_Pt_TP(:, i))  std(CluterFeat_Pt_TP(:, i)); mean(CluterFeat_Pt_FP(:, i))  std(CluterFeat_Pt_FP(:, i)) ]

%%  -   Feature sets
%%      1: RI, 2: PG, 3: CT, 4: SD, 5: CV
%%             Size     mean            std           skewness          Krutosis            Moment           Asymetry      AAL parcel  x,y,z central coord 
%%               1    2 3 4 5 6     7 8 9 10 11    12 13 14 15 16    17 18 19 20 21     22 23 24 25 26    27 28 29 30 31        32           33 34 35
%% Subset_temp = { [1]  [2 3 4 5 6]   [7 8 9 10 11]  [12 13 14 15 16]  [17 18 19 20 21]   [22 23 24 25 26]  [27 28 29 30 31]      [32]         [33 34 35] };
%%  -   Training a classifier with different combination of feature sets and
%%      Test it with leave-one-out cross validation
size_temp = 0;
for c = 1 : size(Subset_temp, 2)
    cbnk = combnk(1:size(Subset_temp, 2), c);
    size_temp = size_temp + size(cbnk, 1);
end

PostriorAll_second = zeros(MaxNumClus, size_temp, numPatient);
count = 1;
Subset_temp_idx = cell2mat(Subset_temp);
for c = 1 : size(Subset_temp, 2)
    cbnk = combnk(1:size(Subset_temp, 2), c);
    for cb = 1 : size(cbnk, 1)
        if(mod(count,10) == 1)
            disp(['Training done percent: ' num2str(round(count/size_temp*100),2) '%'])
        end
        Subset = cell2mat(Subset_temp(cbnk(cb, :)));
        [t, ia, ib] = intersect(Subset, Subset_temp_idx);
        Subset = ib;
        for j = 1:numPatient
            % exclude the test data
            lesionMask_leaveone = bLesion_Pt;
            lesionMask_leaveone(:, j) = [];
            featVect_leaveone = ClusterFeat_Pt;
            featVect_leaveone(:, :, j) = [];
            clear dataLS
            dataLS(:,:) = featVect_leaveone(logical(lesionMask_leaveone(:,1)>0), :, 1);
            for k = 2:(numPatient-1)
                dataLS = cat(1, dataLS, featVect_leaveone(logical(lesionMask_leaveone(:,k)>0), :, k));
            end

            % Classify the test data
            TrainData = dataLS;
            TestData =ClusterFeat_Pt(logical(bLesion_Pt(:,j)>0), :, j);
            TrainClass = lesionMask_leaveone(logical(lesionMask_leaveone>0));
            TrainClass(TrainClass>500) = 1;
            TrainClass(TrainClass>200 & TrainClass<500) = 2;

            [ldaClass, error,POSTERIOR,logp,coeff ]  = classify(TestData(:, Subset), TrainData(:, Subset), TrainClass, 'linear');
            PostriorAll_second(1:length(POSTERIOR(:,2)), count, j) = POSTERIOR(:,2);            
        end
        count = count + 1;
    end
end
PostriorAll_second_temp = PostriorAll_second ./ repmat(max(PostriorAll_second), MaxNumClus, 1);

%%  -   Compute sensitivity and false positivie rate
interval = 0.05;
range2 = 0.0 : interval : 1;

rank_performance = zeros(size_temp, 3, size(range2, 2));
sensitivity_set = zeros(size_temp, numPatient, size(range2, 2));
specificity_set = zeros(size_temp, numPatient, size(range2, 2));
specific_area_set = zeros(size_temp, numPatient, size(range2, 2));
removed_FP_rate_set = zeros(size_temp, numPatient, size(range2, 2));
Subset_set = {};
combnk_set = [];
cluster_set = cell(size_temp, 4, numPatient, size(range2, 2));

for p = 1 : size(range2, 2)
    post_thres_temp = range2(p);

    sensitivity = zeros(1, numPatient);
    specificity = zeros(1, numPatient);
    specific_area = zeros(1, numPatient);
    removed_FP_rate = zeros(1, numPatient);
        
    count = 1;
    for c = 1 : size(Subset_temp, 2)
        cbnk = combnk(1:size(Subset_temp, 2), c);        
        for cb = 1 : size(cbnk, 1)
            Subset = cell2mat(Subset_temp(cbnk(cb, :)));
            [t, ia, ib] = intersect(Subset, Subset_temp_idx);
            Subset = ib;
            for j = 1:numPatient
                TP = sum(bLesion_Pt(:, j)>200 & bLesion_Pt(:, j)<500  & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                FP = sum(bLesion_Pt(:, j)>500 & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                TN = sum(bLesion_Pt(:, j)>500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);
                FN = sum(bLesion_Pt(:, j)>200 & bLesion_Pt(:, j)<500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);

                TP_c = find(bLesion_Pt(:, j)>200 & bLesion_Pt(:, j)<500  & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                FP_c = find(bLesion_Pt(:, j)>500 & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                TN_c = find(bLesion_Pt(:, j)>500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);
                FN_c = find(bLesion_Pt(:, j)>200 & bLesion_Pt(:, j)<500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);
                
                cluster_set(count, :, j, p) = [ {TP_c'}, {FP_c'}, {TN_c'}, {FN_c'} ];

                sensitivity(j) = double(TP>0) / double(sum(bLesion_Pt(:, j)>200 & bLesion_Pt(:, j)<500)>0);
                specificity(j) = TN / ( TN + FP );
                specific_area(j) = sum(ClusterFeat_Pt(bLesion_Pt(:, j)>500 & PostriorAll_second_temp(:, j) < post_thres_temp, 1, j)) / sum(ClusterFeat_Pt(bLesion_Pt(:, j)>500,1,j));
                removed_FP_rate(j) = (TN+Size_FP_clus_num(j))/Total_clus_num(j);
            end
            
            Subset_set{count} =Subset;
            combnk_set(count, :) = [ c, cb ];
            sensitivity_set(count, :, p) = sensitivity;
            specificity_set(count, :, p) = specificity;
            specific_area_set(count, :, p) = specific_area;            
            removed_FP_rate_set(count, :, p) = removed_FP_rate;
            count = count + 1
        end        
    end
end

sensitivity_set_temp = sensitivity_set;
sensitivity_set_temp(:, isnan(mean(mean(sensitivity_set, 3), 1)), :) = 0;
specificity_set_temp = specificity_set;
temp_performance = [ squeeze(mean(sensitivity_set_temp, 2)) squeeze(mean(removed_FP_rate_set, 2)) ];

for p = 1 : size(range2, 2)    
    [b, i] = sort(temp_performance(:, p));
    rank_performance(:, :, p) = [ i, temp_performance(i, [p p + size(range2, 2) ])];
end

performance = zeros(size(range2, 2), 3);
for p = 1 : size(range2, 2)
    temp = rank_performance(:, 2, p)+ rank_performance(:, 3, p);
    [b, i] = max(temp);
    performance(p, :) = rank_performance(i, :, p);
end

%% Final Result
[range2' performance]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters ...
post_thres_temp = 0.90;
m = 165;
p = find(abs(range2-post_thres_temp)<0.00001);
Subset = Subset_set{m}
PostriorAll_second2 = zeros(MaxNumClus, numPatient);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table report ...
[a b] = max(mean(sensitivity_first, 2) - mean(false_positive_first, 2)/81924);
idx = find(range == range(b));
pat_case_temp((sensitivity_first(idx, :)==0))
post_thres = range(b);
sensitivity_first_temp = sensitivity_first(idx, :)
false_positive_first_temp = false_positive_first(idx, :)/89124

TP_set = [];
FP_set = [];
TN_set = [];
FN_set = [];
FP_quant = [];
coverage_set = [];

TP_set_inv = [];
FP_set_inv = [];
FP_set_inv2 = [];

second_TP_set = [];
second_FP_set = [];
second_TN_set = [];
second_FN_set = [];
second_FP_quant = [];
second_coverage_set = [];
second_FP_set_inv = [];

Total_clus = [];
Siz_thres = 0;
for j = 1:numPatient  
    % The first classification ...
    FP_quant_temp = 0;
    TP_set_temp = [];
    FP_set_temp = [];
    coverage_set_temp = [];
    second_TP_set_temp = [];
    second_FP_set_temp = [];
    second_TN_set_temp = [];
    second_FN_set_temp = [];
    second_coverage_set_temp = [];
    
    Total_clus = [ Total_clus max(Clusters(:, j)) ];
    for k = 1 : max(Clusters(:, j))
        Size = sum(Clusters(:, j)==k);
        flag = 0;
        vert_idx = find(Clusters(:, j) == k);
        for v = 1 : size(vert_idx, 1)
            nb_v = surf_data.nbr(vert_idx(v), :);
            nb_v = sort(nb_v(nb_v ~= 0));
            if(isequal(nb_v, intersect(nb_v, vert_idx)))
                flag = 1;
                break;
            end
        end
        cluster_temp = Clusters(:, j) == k;
        if flag == 1
            if(Size > Siz_thres)                            
                if(sum(lesionMasks(:, j) == 1 & cluster_temp == 1))
                    TP_set_temp = [ TP_set_temp sum(cluster_temp) ];
                    coverage_set_temp  = [ coverage_set_temp sum(lesionMasks(:, j) == 1 & cluster_temp == 1) ];
                else
                    FP_quant_temp = FP_quant_temp + 1;
                    FP_set_temp = [ FP_set_temp sum(cluster_temp) ];
                    FP_set_inv   = [ FP_set_inv sum(cluster_temp) ];
                end  
            end
        end
    end

    FP_quant = [ FP_quant mean(FP_quant_temp) ];
    TP_set   = [ TP_set mean(TP_set_temp) ];
    coverage_set = [ coverage_set mean(coverage_set_temp) ];
    TP_set_inv = [ TP_set_inv TP_set_temp ];
    FP_set   = [ FP_set mean(FP_set_temp) ];    
    TN_set   = [ TN_set sum(lesionMasks(:, j) == 0 & Clusters(:, j) == 0) ];
    FN_set   = [ FN_set sum(lesionMasks(:, j) == 1 & Clusters(:, j) == 0) ];
    
    %% The second classficiation to reduce FP ...
    clusid = bLesion_Pt(cluster_set{m, 1, j, p}, j) - 200;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 2)
            second_TP_set_temp = [ second_TP_set_temp sum(Clusters(:, j) == clusid(k)) ];
            second_coverage_set_temp = [ second_coverage_set_temp sum(Clusters(:, j) == clusid(k) & lesionMasks(:, j) == 1) ];
        end
    end
    
    clusid = bLesion_Pt(cluster_set{m, 2, j, p}, j) - 500;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 2)
            second_FP_set_temp = [ second_FP_set_temp sum(Clusters(:, j) == clusid(k)) ];            
        end
    end
    
    clusid = bLesion_Pt(cluster_set{m, 3, j, p}, j) - 500;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 2)
            second_TN_set_temp = [ second_TN_set_temp sum(Clusters(:, j) == clusid(k)) ];
        end
    end
    
    clusid = bLesion_Pt(cluster_set{m, 4, j, p}, j) - 500;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 2)
            second_FN_set_temp = [ second_FN_set_temp sum(Clusters(:, j) == clusid(k)) ];
        end
    end
    
    second_TP_set = [ second_TP_set mean(second_TP_set_temp) ];
    second_FP_set = [ second_FP_set mean(second_FP_set_temp) ];
    second_FP_set_inv = [ second_FP_set_inv  second_FP_set_temp ];
    second_TN_set = [ second_TN_set mean(second_TN_set_temp) ];
    second_FN_set = [ second_FN_set mean(second_FN_set_temp) ];    
    second_FP_quant = [ second_FP_quant size(cluster_set{m, 2, j, p}, 2) ];
    second_coverage_set = [ second_coverage_set mean(second_coverage_set_temp) ]; 
end

%% First classification sensitivity and specificity
mean(sensitivity_first(b, :))
mean(specificity_first(b, :))

%% Lesional cluster mean size and standard deviation
TP_set(isnan(TP_set)) = 0;
mean(TP_set)
std(TP_set)
sum(TP_set)

%% False positive mean size and standard deviation
FP_set(isnan(FP_set)) = 0;
mean(FP_set) 
std(FP_set)

sum(FP_set_inv)
sum(FP_set_inv)/1545163

%% False positive mean amount and standard deviation
mean(FP_quant) 
std(FP_quant)

%% Second classification sensitivity and specificity
sensitivity_second = mean(sensitivity_set_temp(m, :, p))
specificity_second = mean(specificity_set_temp(m, :, p))

%% Second Lesional cluster mean size and standard deviation
second_TP_set(isnan(second_TP_set)) = 0;
mean(second_TP_set)
std(second_TP_set)

second_coverage_set(isnan(second_coverage_set)) = 0;
mean(second_coverage_set)
std(second_coverage_set)

%% Second False positive mean size and standard deviation
second_FP_set(7) = 0;
second_FP_set(isnan(second_FP_set)) = 0;
mean(second_FP_set) 
std(second_FP_set)
sum(second_FP_set)/1545163

%% Second Classification False positive mean amount and standard deviation
second_FP_quant(7) = 0;
second_FP_quant(isnan(second_FP_quant)) = 0;
mean(second_FP_quant(logical(sensitivity_set_temp(m, :, p))))
std(second_FP_quant(logical(sensitivity_set_temp(m, :, p))))

%% Visualization and report final results of classification
OUTPATH = '/local_disk/seokjun/01_project/IntracorticalAnalysis/03_result/lesion_detection_Neurology/15_result/';
temp_zscore_data  = load('db_FCD_3T_k5.mat');
temp_zscore_data  = sum(temp_zscore_data.z_score_database, 3);
temp_zscore_data2 = load('db_FCD_3T_k5.mat');
temp_zscore_data2 = temp_zscore_data2.z_score_database;

cluster_info = cell(numPatient, 5);
contribution_features_TP = zeros(numPatient, 5);
contribution_features_FP = zeros(numPatient, 5, 5);
count = 0;
detected_case = [];
lesionMask = lesionMasks';
for j = 1 : numPatient
%     f = figure; SurfStatView(Clusters(:, j).*MaskCut', TemplateSurf, [ pat_case{j} ', cov:' num2str(Coverage_first(j)) '%, fp:' num2str(false_positive_first(j)*100/81924) '%' ]);
%     exportfigbo(f,[OUTPATH '/' pat_case{j} '_first_classification.png' ], 'png', 6); close(f);
  
    lesion_label = SurfStatReadData({ ['/data/noel/noel6/CTFCD-1.2.0_64/Lesion/lesion_surf/mcd_' pat_case_temp{j} '_label_union_left_rsl.txt'], ...
                                      ['/data/noel/noel6/CTFCD-1.2.0_64/Lesion/lesion_surf/mcd_' pat_case_temp{j} '_label_union_right_rsl.txt'] });
    cluster_temp = zeros(1, 81924);
    cluster_info{j, 5} = 0;
    %% True postives
    for k = 1 : size(cluster_set{m, 1, j, p}, 2)
        if(k == 1)
            detected_case = [ detected_case j ];
        end
        
        cluster_temp(Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200)) = p+1;        
        cluster_info{j, 1} =  [ cluster_info{j, 1} sum(Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200)) ];
        cluster_info{j, 3} =  [ cluster_info{j, 3} mean(power((temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 1)).^2 + ...
                                                              (temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 2)).^2 + ...
                                                              (temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 3)).^2 + ...
                                                              (temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 4)).^2 + ...
                                                              (temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 5)).^2, 0.5)) ];
        contribution_features_TP(j, :) = [ mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 1))) ...
                                           mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 2))) ... 
                                           mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 3))) ...
                                           mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 4))) ...
                                           mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200), 5))) ];
                                       
                                        
        cluster_info{j, 5} = cluster_info{j, 5} + (sum((Clusters(:, j) == (bLesion_Pt(cluster_set{m, 1, j, p}(k), j) - 200))' & (lesion_label==1))/sum((lesion_label==1)));
    end
    
    %% False postives
    for k = 1 : size(cluster_set{m, 2, j, p}, 2)
        cluster_temp(Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500)) = k;        
        cluster_info{j, 2} =  [ cluster_info{j, 2} sum(Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500)) ];
        cluster_info{j, 4} =  [ cluster_info{j, 4} mean(power((temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 1)).^2 + ...
                                                              (temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 2)).^2 + ...
                                                              (temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 3)).^2 + ...
                                                              (temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 4)).^2 + ...
                                                              (temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 5)).^2, 0.5)) ];
        contribution_features_FP(j, :, k) = [ mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 1))) ...
                                              mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 2))) ... 
                                              mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 3))) ...
                                              mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 4))) ...
                                              mean(abs(temp_zscore_data2(j, Clusters(:, j) == (bLesion_Pt(cluster_set{m, 2, j, p}(k), j) - 500), 5))) ];
    end
    
    count = count + sum(cluster_info{j, 2});
    
%     f = figure; SurfStatView(sum(featVect(:, :, j), 2).*MaskCut', TemplateSurf, [ pat_case{j} ]); SurfStatColLim([-5 5]);
%     exportfigbo(f,[OUTPATH '/' pat_case{j} '_composite_scoremap.png' ], 'png', 6); close(f);
%     f = figure; SurfStatView(featVect(:, 4, j).*MaskCut', TemplateSurf, [ pat_case{j} ]); SurfStatColLim([-5 5]);
%     exportfigbo(f,[OUTPATH '/' pat_case{j} '_composite_scoremap.png' ], 'png', 6); close(f);
%     f = figure; SurfStatView(cluster_temp.*MaskCut, TemplateSurf, [ pat_case{j} ', specificity: ' num2str(specificity_set(m, j)) '%' ]);
%     SurfStatColLim([0 5]);
%     f = figure; SurfStatView1(temp_zscore_data2(j, :, 2).*MaskCut,TemplateSurf);    SurfStatColLim([-5 5]); cameramenu;
%     f = figure; SurfStatView1(cluster_temp.*MaskCut, TemplateSurf);  SurfStatColLim([0 5]); cameramenu;
%     exportfigbo(f,[OUTPATH '/' pat_case{j} '_second_classification.png' ], 'png', 6); close(f);
end

cluster_info{7, 2} = [];
cluster_info{7, 4} = [];
cluster_info = [ transpose(num2cell(1:numPatient)) pat_case_temp' histo_type_temp cluster_info ]

coverage_set = cell2mat(cluster_info(:, 8)); coverage_set_temp = coverage_set(detected_case);
histo_type_temp2 = histo_type_temp(detected_case);
[ mean(coverage_set_temp(strcmp(histo_type_temp2, 'FCDIIb'))) std(coverage_set_temp(strcmp(histo_type_temp2, 'FCDIIb'))) ]
[ mean(coverage_set_temp(strcmp(histo_type_temp2, 'FCDIIa'))) std(coverage_set_temp(strcmp(histo_type_temp2, 'FCDIIa'))) ]
mean(coverage_set_temp)
