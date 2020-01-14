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
%  3) This script is for testing generalizability of our 3T-based classifier on different field strength and scanners (i.e., independent second dataset of 1.5T MRI)
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
load('/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/lesion_detection/Result_FCD_1.5T.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. variable initialization - read files
Group_pat15T              = 'FCD';
Prefix_pat15T             = 'mcd';
Cases_pat15T              = '/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/Demographic_data_FCD_1.5T.txt';
Group_pat                 = 'FCD';
Prefix_pat                = 'mcd';
Cases_pat                 = '/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/Demographic_data_FCD.txt';
Left_Right                = 'both';
NumMesh                   = 81920;
KERNEL                    = 5;
OUTPATH                   = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/lesion_detection_frontiers/06_result/';

fid = fopen(Cases_pat15T);
demo = textscan(fid, '%s%f%f', 'Delimiter', ',', 'CollectOutput', 1);
case_num_pat15T = demo{1};
gender_pat15T = demo{2}(:, 1); gender_pat15T_temp = cell(size(gender_pat15T, 1), 1); gender_pat15T_temp(gender_pat15T==0) = {'m'}; gender_pat15T_temp(gender_pat15T==1) = {'f'}; gender_pat15T = gender_pat15T_temp;
age_pat15T = demo{2}(:, 2);
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
z_score_database_FCD = z_score_database;
lesion_db_FCD        = lesion_db;

load([ 'db_FCD_1.5T_k5.mat' ]);
load('db_lesion_1.5T.mat');
z_score_database_FCD15 = z_score_database;
lesion_db_FCD15        = lesion_db_fcd15T;

mean_age_pat = mean(age_pat);
std_age_pat = std(age_pat, 0);
mean_age_pat15T = mean(age_pat15T);
std_age_pat15T = std(age_pat15T, 0);

[h,p,ci,stats] = ttest2(age_pat15T,age_pat);

pat_case    = { '027', '049', '061', '062', '063', '066', '071', '072', '073', '074', '076', '079', '080_1', '081', '082', '083', '084', '085', '086' };
pat15T_case = { '001', '002', '004', '006', '01215T', '014', '018', '019', '023', '024', '026', '027', '028', '029', '030', '031', '032', ...
                '035', '036', '039', '041', '046', '047', '048', '04915T', '050', '051', '05515T', '05715T', '05815T', '05915T', '06015T', '06115T' };
lesion_vol = [       10702       25098        1870        4785       24055        1272        3325        6424        2941        4583        1949 ...
                     5432        1539        6796       11727        8084        2033         728        4786        2933       12076        1446  ...                   
                     741         898         217        1902        8376        4619         309        1415        2469         285        1024  ];

[C, pat_idx, ib]  = intersect(case_num_pat, pat_case);
[C, pat15T_idx, ib] = intersect(case_num_pat15T, pat15T_case);

mean_age_pat15T_new = mean(age_pat15T(pat15T_idx));
std_age_pat15T_new = std(age_pat15T(pat15T_idx), 0);
sum(strcmp(gender_pat15T(pat15T_idx), 'm'))

mean_age_pat_new = mean(age_pat(pat_idx));
std_age_pat_new = std(age_pat(pat_idx), 0);
sum(strcmp(gender_pat(pat_idx), 'm'))

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

outlier = [  ];              %% Cases that will be excluded in this analysis
pat15T_case(outlier)
pat15T_case_temp = pat15T_case;
pat15T_case_temp(outlier) =[];

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

lesionMasksFCD     = permute(lesion_db_FCD, [2 1]);
lesionMasksFCD15T  = permute(lesion_db_FCD15, [2 1]);

lesion_sizeFCD15T_org = zeros(1, size(pat15T_case_temp, 2));
for i = 1 : size(pat15T_case_temp, 2)
    lesion_sizeFCD15T_org(i) = sum(lesionMasksFCD15T(:, i) > 0);
end

lesionMaskTemp     = SurfStatSmooth(lesionMasksFCD', TemplateSurf, smooth_kernel);
lesionMasksFCD     = double(lesionMaskTemp>binarize_thres)';
lesionMaskTemp     = SurfStatSmooth(lesionMasksFCD15T', TemplateSurf, smooth_kernel);
lesionMasksFCD15T  = double(lesionMaskTemp>binarize_thres)';

lesion_sizeFCD = zeros(1, size(pat_case_temp, 2));
for i = 1 : size(pat_case_temp, 2)
    lesion_sizeFCD(i) = sum(lesionMasksFCD(:, i) == 1);
end
lesion_sizeFCD15T = zeros(1, size(pat15T_case_temp, 2));
for i = 1 : size(pat15T_case_temp, 2)
    lesion_sizeFCD15T(i) = sum(lesionMasksFCD15T(:, i) == 1);
end

featVectFCD = permute(z_score_database_FCD, [2 3 1]);
featVectFCD(:, :, outlier) = [];
lesionMasksFCD(:, outlier) = [];
lesionMasksFCD(MaskCut<1,:)=2;

featVectFCD15T = permute(z_score_database_FCD15, [2 3 1]);
featVectFCD15T(:, :, outlier) = [];
lesionMasksFCD15T(:, outlier) = [];
lesionMasksFCD15T(MaskCut<1,:)=2;

siz = size(featVectFCD);
numVertFCD    = siz(1);
numFeatFCD    = siz(2);
numPatientFCD = siz(3);

siz = size(featVectFCD15T);
numVertFCD15T    = siz(1);
numFeatFCD15T    = siz(2);
numPatientFCD15T = siz(3);

HTissueMaskFCD = lesionMasksFCD;
HTissueMaskFCD(:,:) = 0;
sizLesionFCD = sum(sum(lesionMasksFCD == 1));
sizHTissueperSubjFCD = round(sizLesionFCD / numPatientFCD) * sampling_scale;
CandidatesFCD = findn(lesionMasksFCD == 0);

%% Random healthy tissue sampling while balancing its sample size with
%% lesional tissue
for k = 1:numPatientFCD
    sizeFullHTissueFCD = sum(CandidatesFCD(:,2)==k);
    CandEachPatientFCD =CandidatesFCD(CandidatesFCD(:,2)==k,1);
    SampleIndex = rand(sizeFullHTissueFCD,1) > 1 - sizHTissueperSubjFCD / sizeFullHTissueFCD;
    sum(SampleIndex~=0)
    HTissueMaskFCD(CandEachPatientFCD(SampleIndex),k) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. First classification: Leave-one-out
PostriorAll = zeros(numVertFCD15T, numPatientFCD15T);
for j = 1:numPatientFCD15T
    clear dataLS dataHT
    dataLS(:,:) = featVectFCD(logical(lesionMasksFCD(:,1)==1), :, 1);
    for k = 2:(numPatientFCD-1)
        dataLS = cat(1, dataLS, featVectFCD(logical(lesionMasksFCD(:,k)==1), :, k));
    end
    dataHT(:,:) = featVectFCD(logical(HTissueMaskFCD(:,1)), :, 1);
    for k = 2:(numPatientFCD-1)
        dataHT = cat(1, dataHT, featVectFCD(logical(HTissueMaskFCD(:,k)), :, k));
    end
    
    % Classify the test data
    TrainData = cat(1, dataHT, dataLS);
    TestData = featVectFCD15T(:, :, j);
    TrainClass = cat(2, ones(1, length(dataHT)), ones(1, length(dataLS))*2)';
    
    [ldaClass, error,POSTERIOR, logp, coeff ]  = classify(TestData(:,FeatSubSet), TrainData(:,FeatSubSet), TrainClass, 'diagquadratic');
    PostriorAll(:, j) = POSTERIOR(:,2);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. First classification: compute sensitivity and false positivie rate
interval = 0.001;
range = 0.5 : interval : 1;
sensitivity_first = zeros(size(range, 2), numPatientFCD15T);
specificity_first = zeros(size(range, 2), numPatientFCD15T);
Coverage_first = zeros(size(range, 2), numPatientFCD15T);
Coverage = zeros(size(range, 2), numPatientFCD15T);
false_positive_first = zeros(size(range, 2), numPatientFCD15T);
false_positive_first_rate = zeros(size(range, 2), numPatientFCD15T);
for post_thres = range
    for j = 1:numPatientFCD15T
        idx = int32((post_thres-0.5)/interval+1);
        TP = sum(lesionMasksFCD15T(:, j)==1 & PostriorAll(:, j) >= post_thres & AAL_data' ~= 0);
        FP = sum(lesionMasksFCD15T(:, j)==0 & PostriorAll(:, j) >= post_thres & AAL_data' ~= 0);
        TN = sum(lesionMasksFCD15T(:, j)==0 & PostriorAll(:, j) < post_thres & AAL_data' ~= 0);
        FN = sum(lesionMasksFCD15T(:, j)==1 & PostriorAll(:, j) < post_thres & AAL_data' ~= 0);

        sensitivity_first(idx,j)= TP>0;
        specificity_first(idx,j) = TN / sum(lesionMasksFCD15T(:, j)==0);
        Coverage_first(idx,j) = TP / sum(lesionMasksFCD15T(:,j)==1);
        Coverage (idx, j) = TP;
        false_positive_first(idx,j) = FP;
        false_positive_first_rate(idx, j) = FP / (FP + TN);
    end
end
[a b] = max(mean(sensitivity_first, 2) - mean(false_positive_first, 2)/81924);

% Visualization of the vertex-wise classifier performance
% interval = 0.001;
% range = 0.5 : interval : 1;
% range_temp = int32(((0.5 - range(1))/interval + 1) : ((1 - range(1))/interval + 1));
% figure; 
% [ AX, H1, H2 ] = plotyy(range(range_temp), mean(sensitivity_first(range_temp, :), 2)/mean(sensitivity_first(b, :)), ...
%                         range(range_temp), mean(false_positive_first(range_temp, :), 2)/mean(false_positive_first(b, :)), 'plot');
% set(get(AX(1),'Ylabel'),'String','Ratio Detection rate'); 
% set(get(AX(2),'Ylabel'),'String','Ratio FP'); 
% xlabel('Threshold for posterior');
% title('Tradeoff between detection rate and FP');
% set(H1,'LineStyle','--');
% set(H2,'LineStyle','--');
% bottom = 0; upper = 2; interval_axis = 0.5;
% set(AX(1), 'YTick', [ bottom : interval_axis : upper ]);
% set(AX(1), 'YLim',  [ bottom upper ]);
% set(AX(1), 'YTick', [ bottom : interval_axis : upper ]);
% set(AX(2), 'YLim',  [ bottom upper ]);
% set(AX(2), 'YTick', [ bottom: interval_axis : upper ]);
% grid on; hold on; scatter(range(b), 1.1, 'rv');
% 
% figure; plot(mean(false_positive_first(range_temp, :), 2)/max(mean(false_positive_first(range_temp, :), 2)), ...
%              mean(sensitivity_first(range_temp, :), 2)/max(mean(sensitivity_first(range_temp, :), 2)));
% figure; plot(mean(false_positive_first_rate(range_temp, :), 2), mean(sensitivity_first(range_temp, :), 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. Second classification initialization
idx = find(range == range(b));

% Missing cases
pat15T_case_temp((sensitivity_first(idx, :)==0))

% Sensitivity and false positive ratio of the vertex-wise classifier
sensitivity_first_temp = sensitivity_first(idx, :)
false_positive_first_temp = false_positive_first(idx, :)/89124
post_thres = range(b);

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
lesionClassFCD15T = PostriorAll >= post_thres;
MaxNumClus = 500;
bLesion_PtFCD15T = zeros(MaxNumClus, numPatientFCD15T);
ClusterFeat_PtFCD15T = zeros(MaxNumClus, size(cell2mat(Subset_temp), 2), numPatientFCD15T);
featVect_db = [];
Total_clus_numFCD15T = zeros(1, numPatientFCD15T);
Size_FP_clus_numFCD15T = zeros(1, numPatientFCD15T);
FP_clus_numFCD15T = zeros(MaxNumClus, numPatientFCD15T);
ClustersFCD15T = zeros(81924, numPatientFCD15T);
order=6;

count_temp = [];
OUTPATH = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/lesion_detection_frontiers/14_result/';

%% Detect the clusters and remove FPs according to following criteria:
%% 1) FP cluster which doesn't meet the condition of 6-connected cluster 
%%    This is for removing the peculiar shape of clusters, for instance, zigzaged elongated cluster.
%% 2) FP cluster which part is in the mid mask area

%% Collect feature distribution of clusters after prunning process separately in lesion and FP clusters
for j = 1 : numPatientFCD15T
    j
    l=1;
    
    ClustersFCD15T(:,j) = SurfStatCluster(lesionClassFCD15T(:,j)', TemplateSurf);        
    Total_clus_numFCD15T(j) = max(ClustersFCD15T(:,j));
    
    for k = 1 : Total_clus_numFCD15T(j)
        Size = sum(ClustersFCD15T(:,j)==k);
        flag = 0;
        vert_idx = find(ClustersFCD15T(:, j) == k);
        
        %% Define cluster ... 
        for v = 1 : size(vert_idx, 1)
            nb_v = surf_data.nbr(vert_idx(v), :);
            nb_v = sort(nb_v(nb_v ~= 0));
            if(isequal(nb_v, intersect(nb_v, vert_idx)))
                flag = 1;
                break;
            end
        end
                
        if flag == 1
            if(Size > Siz_thres)
                if sum(lesionMasksFCD15T(:, j)==1 & ClustersFCD15T(:,j)==k)                                    
                        bLesion_PtFCD15T(l,j) = 200 + k;                                                                                          
                else
                    %% Remove clusters in medial wall ...
                    if(sum(ClustersFCD15T(:,j)==k & AAL_data' == 0))
                        ClustersFCD15T(ClustersFCD15T(:,j)==k, j) = 0;                        
                        Size_FP_clus_numFCD15T(j) = Size_FP_clus_numFCD15T(j) + 1;
                        continue;
                    else
                        bLesion_PtFCD15T(l,j) = 500 + k;
                        FP_clus_numFCD15T(l,j) = FP_clus_numFCD15T(l,j) + Size;                        
                    end
                end             
                
                mahal_dist = sum(featVectFCD15T(vert_idx, :, j).^2, 2);
                [maxPost, id] = max(PostriorAll(vert_idx, j));
                id = find(PostriorAll(vert_idx, j) == maxPost);
                
                [maxMahal, id_maxMahal] = max(mahal_dist(id));
                temp = vert_idx(id(id_maxMahal));                
                peakCluster = [ temp; find_vertex_multineighbors_by_order(TemplateSurf,temp,order,edg) ]';

                M1 = mean(featVectFCD15T(peakCluster,:,j), 1);
                M2 = std(featVectFCD15T(peakCluster,:,j), 1);
                M3 = skewness(featVectFCD15T(peakCluster,:,j), 1);
                M4 = kurtosis(featVectFCD15T(peakCluster,:,j), 1);
                M5 = moment(featVectFCD15T(peakCluster,:,j), 5);

                Cx = TemplateSurf.coord(1, peakCluster);
                Cy = TemplateSurf.coord(2, peakCluster);
                Cz = TemplateSurf.coord(3, peakCluster);

                TempDist = dist3(mean([Cx; Cy; Cz],2)', [Cx; Cy; Cz]');
                nn = find(TempDist == min(TempDist));
                finalModCoord = [Cx(nn(1)); Cy(nn(1)); Cz(nn(1))];
                vertNN = find(TemplateSurf.coord(1,:)==finalModCoord(1) & TemplateSurf.coord(2,:)==finalModCoord(2)  & TemplateSurf.coord(3,:)==finalModCoord(3));
                AAL_vertNN = AAL_data(vertNN);
                TTT=find(ClustersFCD15T(:,j)==k);
                TTT2 = 40962+TTT;
                if sum(TTT2 <= 81924 & TTT2 > 40962)
                    Asym = mean(featVectFCD15T(TTT,:,j)) - mean(featVectFCD15T(TTT2,:,j));
                else
                    TTT2 = TTT - 40962;
                    Asym = mean(featVectFCD15T(TTT,:,j)) - mean(featVectFCD15T(TTT2,:,j));
                end            
                ClusterFeat_PtFCD15T(l, :, j) = cat(2, Size, M1, M2, M3, M4, M5, Asym, AAL_vertNN, finalModCoord');
                l = l+1;
            else       
                ClustersFCD15T(vert_idx, j) = 0;            
                Size_FP_clus_numFCD15T(j) = Size_FP_clus_numFCD15T(j) + 1;
            end
        else
            ClustersFCD15T(vert_idx, j) = 0;            
            Size_FP_clus_numFCD15T(j) = Size_FP_clus_numFCD15T(j) + 1;
        end
    end
end

%% Estimate the feature pattern in TP and FP clusters
% CluterFeat_Pt_TP = [];
% CluterFeat_Pt_FP = [];
% for j = 1 : numPatientFCD15T
%     for k = 1 : size(bLesion_PtFCD15T(:, j), 1)
%         if(bLesion_PtFCD15T(k, j) < 300 & bLesion_PtFCD15T(k, j) > 200)
%             CluterFeat_Pt_TP = [ CluterFeat_Pt_TP; ClusterFeat_PtFCD15T(k, :, j) ];
%         elseif(bLesion_PtFCD15T(k, j) > 500)
%             CluterFeat_Pt_FP = [ CluterFeat_Pt_FP; ClusterFeat_PtFCD15T(k, :, j) ];
%         end
%     end
% end
% mean(CluterFeat_Pt_TP, 1)
% mean(CluterFeat_Pt_FP, 1)

%%  -   Feature sets
%%      1: RI, 2: PG, 3: CT, 4: SD, 5: CV
%%             Size     mean            std           skewness          Krutosis            Moment           Asymetry      AAL parcel  x,y,z central coord 
%%               1    2 3 4 5 6     7 8 9 10 11    12 13 14 15 16    17 18 19 20 21     22 23 24 25 26    27 28 29 30 31        32           33 34 35
%% Subset_temp = { [1]  [2 3 4 5 6]   [7 8 9 10 11]  [12 13 14 15 16]  [17 18 19 20 21]   [22 23 24 25 26]  [27 28 29 30 31]      [32]         [33 34 35] };

%% We used the classifier trained by 3T data of FCD patients to test 1.5T FCD
load(['/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/lesion_detection/db_FCD_3T_k5_for_second_classifier.mat']);

%%  -   Training a classifier with different combination of feature sets and
%%      Test it with leave-one-out cross validation

bLesion_PtFCD = bLesion_Pt;
ClusterFeat_PtFCD = ClusterFeat_Pt;

size_temp = 0;
for c = 1 : size(Subset_temp, 2)
    cbnk = combnk(1:size(Subset_temp, 2), c);
    size_temp = size_temp + size(cbnk, 1);
end

PostriorAll_second = zeros(MaxNumClus, size_temp, numPatientFCD15T);
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
        for j = 1:numPatientFCD15T            
            clear dataLS
            dataLS(:,:) = ClusterFeat_PtFCD(logical(bLesion_PtFCD(:,1)>0), :, 1);
            for k = 2:numPatientFCD
                dataLS = cat(1, dataLS, ClusterFeat_PtFCD(logical(bLesion_PtFCD(:,k)>0), :, k));
            end

            % Classify the test data
            TrainData = dataLS;
            TestData =ClusterFeat_PtFCD15T(logical(bLesion_PtFCD15T(:,j)>0), :, j);
            TrainClass = bLesion_PtFCD(logical(bLesion_PtFCD>0));
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
interval = 0.01;
range2 = 0.0 : interval : 1;

rank_performance = zeros(size_temp, 3, size(range2, 2));
sensitivity_set = zeros(size_temp, numPatientFCD15T, size(range2, 2));
specificity_set = zeros(size_temp, numPatientFCD15T, size(range2, 2));
specific_area_set = zeros(size_temp, numPatientFCD15T, size(range2, 2));
removed_FP_rate_set = zeros(size_temp, numPatientFCD15T, size(range2, 2));
Subset_set = {};
combnk_set = [];
cluster_set = cell(size_temp, 4, numPatientFCD15T, size(range2, 2));

for p = 1 : size(range2, 2)
    post_thres_temp = range2(p);

    sensitivity = zeros(1, numPatientFCD15T);
    specificity = zeros(1, numPatientFCD15T);
    specific_area = zeros(1, numPatientFCD15T);
    removed_FP_rate = zeros(1, numPatientFCD15T);
        
    count = 1;
    for c = 1 : size(Subset_temp, 2)
        cbnk = combnk(1:size(Subset_temp, 2), c);        
        for cb = 1 : size(cbnk, 1)
            Subset = cell2mat(Subset_temp(cbnk(cb, :)));
            [t, ia, ib] = intersect(Subset, Subset_temp_idx);
            Subset = ib;
            for j = 1:numPatientFCD15T
                TP = sum(bLesion_PtFCD15T(:, j)>200 & bLesion_PtFCD15T(:, j)<500  & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                FP = sum(bLesion_PtFCD15T(:, j)>500 & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                TN = sum(bLesion_PtFCD15T(:, j)>500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);
                FN = sum(bLesion_PtFCD15T(:, j)>200 & bLesion_PtFCD15T(:, j)<500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);

                TP_c = find(bLesion_PtFCD15T(:, j)>200 & bLesion_PtFCD15T(:, j)<500  & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                FP_c = find(bLesion_PtFCD15T(:, j)>500 & PostriorAll_second_temp(:, count, j) > post_thres_temp);
                TN_c = find(bLesion_PtFCD15T(:, j)>500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);
                FN_c = find(bLesion_PtFCD15T(:, j)>200 & bLesion_PtFCD15T(:, j)<500 & PostriorAll_second_temp(:, count, j) < post_thres_temp);
                
                cluster_set(count, :, j, p) = [ {TP_c'}, {FP_c'}, {TN_c'}, {FN_c'} ];

                sensitivity(j) = double(TP>0) / double(sum(bLesion_PtFCD15T(:, j)>200 & bLesion_PtFCD15T(:, j)<500)>0);
                specificity(j) = TN / ( TN + FP );
                specific_area(j) = sum(ClusterFeat_PtFCD15T(bLesion_PtFCD15T(:, j)>500 & PostriorAll_second_temp(:, j) < post_thres_temp, 1, j)) / sum(ClusterFeat_PtFCD15T(bLesion_PtFCD15T(:, j)>500,1,j));
                removed_FP_rate(j) = (TN+Size_FP_clus_numFCD15T(j))/Total_clus_numFCD15T(j);
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

%% Results ...
[range2' performance]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters ...
post_thres_temp = 0.94;
m = 368;
p = find(abs(range2-post_thres_temp)<0.00001);
Subset = Subset_set{m}
PostriorAll_second2 = zeros(MaxNumClus, numPatientFCD15T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table report ...
[a b] = max(mean(sensitivity_first, 2) - mean(false_positive_first, 2)/81924);
idx = find(range == range(b));
pat15T_case_temp((sensitivity_first(idx, :)==0))
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
for j = 1:numPatientFCD15T 
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
    
    Total_clus = [ Total_clus max(ClustersFCD15T(:, j)) ];
    for k = 1 : max(ClustersFCD15T(:, j))
        Size = sum(ClustersFCD15T(:, j)==k);
        flag = 0;
        vert_idx = find(ClustersFCD15T(:, j) == k);
        for v = 1 : size(vert_idx, 1)
            nb_v = surf_data.nbr(vert_idx(v), :);
            nb_v = sort(nb_v(nb_v ~= 0));
            if(isequal(nb_v, intersect(nb_v, vert_idx)))
                flag = 1;
                break;
            end
        end
        cluster_temp = ClustersFCD15T(:, j) == k;
        if flag == 1
            if(Size > Siz_thres)                            
                if(sum(lesionMasksFCD15T(:, j) == 1 & cluster_temp == 1))
                    TP_set_temp = [ TP_set_temp sum(cluster_temp) ];
                    coverage_set_temp  = [ coverage_set_temp sum(lesionMasksFCD15T(:, j) == 1 & cluster_temp == 1) ];
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
    TN_set   = [ TN_set sum(lesionMasksFCD15T(:, j) == 0 & ClustersFCD15T(:, j) == 0) ];
    FN_set   = [ FN_set sum(lesionMasksFCD15T(:, j) == 1 & ClustersFCD15T(:, j) == 0) ];
    
    %% The second classficiation to reduce FP ...
    clusid = bLesion_PtFCD15T(cluster_set{m, 1, j, p}, j) - 200;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 2)
            second_TP_set_temp = [ second_TP_set_temp sum(ClustersFCD15T(:, j) == clusid(k)) ];
            second_coverage_set_temp = [ second_coverage_set_temp sum(ClustersFCD15T(:, j) == clusid(k) & lesionMasksFCD15T(:, j) == 1) ];
        end
    end
    
    clusid = bLesion_PtFCD15T(cluster_set{m, 2, j, p}, j) - 500;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 2)
            second_FP_set_temp = [ second_FP_set_temp sum(ClustersFCD15T(:, j) == clusid(k)) ];            
        end
    end
    
    clusid = bLesion_PtFCD15T(cluster_set{m, 3, j, p}, j) - 500;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 2)
            second_TN_set_temp = [ second_TN_set_temp sum(ClustersFCD15T(:, j) == clusid(k)) ];
        end
    end
    
    clusid = bLesion_PtFCD15T(cluster_set{m, 4, j, p}, j) - 500;
    if(~isempty(clusid))
        for k = 1 : size(clusid, 2)
            second_FN_set_temp = [ second_FN_set_temp sum(ClustersFCD15T(:, j) == clusid(k)) ];
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
second_FP_set(isnan(second_FP_set)) = 0;
mean(second_FP_set) 
std(second_FP_set)
sum(second_FP_set)/1545163

%% Second Classification False positive mean amount and standard deviation
second_FP_quant(isnan(second_FP_quant)) = 0;
mean(second_FP_quant(logical(sensitivity_set_temp(m, :, p))))
std(second_FP_quant(logical(sensitivity_set_temp(m, :, p))))

[a c] = sort(lesion_vol);
[a c] = sort(lesion_sizeFCD15T_org);
pat15T_case_temp(c)
[ lesion_vol(c); lesion_sizeFCD15T_org(c); sensitivity_set_temp(m, c, p); second_FP_quant(c) ]

idx_range = [2 4 : 5  8 : 18];
pat15T_case_temp(c(idx_range))
[ lesion_vol(c(idx_range)); lesion_sizeFCD15T_org(c(idx_range)); sensitivity_set_temp(m, c(idx_range), p); second_FP_quant(c(idx_range)) ]
sum(sensitivity_set_temp(m, c(idx_range), p))/length(idx_range)

mean(sensitivity_first(b, c(idx_range)))
mean(second_FP_quant(c(idx_range)))
std(second_FP_quant(c(idx_range)))
mean(FP_quant(c(idx_range)))
std(FP_quant(c(idx_range)))

%% Visualization and report final results of classification
OUTPATH = '/local_raid/seokjun/01_project/IntracorticalAnalysis/03_result/lesion_detection_frontiers/17_result_1.5T_FCD/';
temp_zscore_data = load('db_FCD_1.5T_k5.mat');
temp_zscore_data = sum(temp_zscore_data.z_score_database, 3);
temp_zscore_data2 = load('db_FCD_1.5T_k5.mat');
temp_zscore_data2 = temp_zscore_data2.z_score_database;

cluster_info = cell(numPatientFCD15T, 4);
contribution_features_TP = zeros(numPatientFCD15T, 5);
contribution_features_FP = zeros(numPatientFCD15T, 5, 5);
count = 0;
for j = 1 : numPatientFCD15T
    cluster_temp = zeros(1, 81924);
    %% True postives
    for k = 1 : size(cluster_set{m, 1, j, p}, 2)
        cluster_temp(ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200)) = p+1;        
        cluster_info{j, 1} =  [ cluster_info{j, 1} sum(ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200)) ];
        cluster_info{j, 3} =  [ cluster_info{j, 3} mean(power((temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 1)).^2 + ...
                                                              (temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 2)).^2 + ...
                                                              (temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 3)).^2 + ...
                                                              (temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 4)).^2 + ...
                                                              (temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 5)).^2, 0.5)) ];
        contribution_features_TP(j, :) = [ mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 1))) ...
                                           mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 2))) ... 
                                           mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 3))) ...
                                           mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 4))) ...
                                           mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 1, j, p}(k), j) - 200), 5))) ];
            
    end
    
    %% False postives
    for k = 1 : size(cluster_set{m, 2, j, p}, 2)
        cluster_temp(ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500)) = k;        
        cluster_info{j, 2} =  [ cluster_info{j, 2} sum(ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500)) ];
        cluster_info{j, 4} =  [ cluster_info{j, 4} mean(power((temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 1)).^2 + ...
                                                              (temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 2)).^2 + ...
                                                              (temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 3)).^2 + ...
                                                              (temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 4)).^2 + ...
                                                              (temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 5)).^2, 0.5)) ];
        contribution_features_FP(j, :, k) = [ mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 1))) ...
                                              mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 2))) ... 
                                              mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 3))) ...
                                              mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 4))) ...
                                              mean(abs(temp_zscore_data2(j, ClustersFCD15T(:, j) == (bLesion_PtFCD15T(cluster_set{m, 2, j, p}(k), j) - 500), 5))) ];
    end
    
    count = count + sum(cluster_info{j, 2});
    
%     f = figure; SurfStatView(Clusters(:, j).*MaskCut', TemplateSurf, [ pat_case{j} ', cov:' num2str(Coverage_first(j)) '%, fp:' num2str(false_positive_first(j)*100/81924) '%' ]);
%     exportfigbo(f,[OUTPATH '/' pat_case{j} '_first_classification.png' ], 'png', 6); close(f);    
%    f = figure; SurfStatView(cluster_temp.*MaskCut, TemplateSurf);  SurfStatColLim([0 1]); 
%    exportfigbo(f,[OUTPATH '/' pat15T_case_temp{j} '_second_classification.png' ], 'png', 6); close(f);
end

cluster_info_temp1 = [ transpose(num2cell(1:numPatientFCD15T)) num2cell(lesion_vol') pat15T_case_temp' cluster_info ]
cluster_info_temp2 = [ transpose(num2cell(c(idx_range))) num2cell(lesion_vol(c(idx_range)))' pat15T_case_temp(c(idx_range))' cluster_info(c(idx_range), :) ]
cluster_info_temp2 = [ transpose(num2cell(sort(c(idx_range)))) num2cell(lesion_vol(sort(c(idx_range))))' pat15T_case_temp(sort(c(idx_range)))' cluster_info(sort(c(idx_range)), :) ]
