clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OutDir                    = '/local_raid/seokjun/01_project/IntracorticalAnalysis/03_result/'
Group_cont                = 'control'
Prefix_cont               = 'TLE'
Cases_cont                = '/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/Demographic_data_control_1.5T.txt'
Group_pat                 = 'FCD'
Prefix_pat                = 'mcd'
Cases_pat                 = '/local_raid/seokjun/01_project/IntracorticalAnalysis/01_analysis/Demographic_data_FCD_1.5T.txt'

Left_Right                = 'both'
NumIntSurf                = 3
NumSubSurf                = 3
NumMesh                   = 81920
SamplingSpace             = 'native'
img_contrast              = {'t1'}
Kernel                    = 5
Parametric                = 'quadratic'
average_surface_dir       = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/average_surfaces/'

visualization = 1;
load('/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');

fid = fopen(Cases_cont);
demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
case_num_cont = demo{1};
fclose(fid);

fid = fopen(Cases_pat);
demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
case_num_pat = demo{1};
fclose(fid);

excluded_control_cases = { '201', '204', '205', '207', '209', '211', '214', '215', '218', '220', '221', '222', '228', '234', '223', '235', '236', '237', '241', '243', '250', '258', '261' }';
[C, ia, ib] = intersect(case_num_cont, excluded_control_cases);
case_num_cont_new = case_num_cont; case_num_cont_new(ia) = [];
case_num_cont = case_num_cont_new;

excluded_patient_cases = { '016', '020', '5315T' }';
[C, ia, ib] = intersect(case_num_pat, excluded_patient_cases);
case_num_pat_new = case_num_pat; case_num_pat_new(ia) = [];
case_num_pat = case_num_pat_new;

%% Parameters related to the analysis
if(NumMesh == 81920)
    vertexnum = 81924;
elseif(NumMesh == 327680)
    vertexnum = 327684;
end

if(strcmp(Left_Right, 'left'))    
    vertexnum = vertexnum/2;
    Left_Right = { 'left'; 1:vertexnum };
elseif(strcmp(Left_Right, 'right')) 
    vertexnum = vertexnum/2;
    Left_Right = { 'right'; vertexnum+1:(vertexnum*2) };    
elseif(strcmp(Left_Right, 'both'))
    Left_Right = {'left', 'right'; 1:vertexnum/2 vertexnum/2+1:vertexnum };
end

TemplateSurf = cell(NumSubSurf+NumIntSurf+2, 1);
TemplateSurf{1} = SurfStatReadSurf({[ average_surface_dir 'average_gray_surface_left_',  num2str(NumMesh) '_t1.obj' ],  ...
                                    [ average_surface_dir 'average_gray_surface_right_', num2str(NumMesh) '_t1.obj' ] });
TemplateSurf{2} = SurfStatReadSurf({[ average_surface_dir 'average_intracortical_surface_1_left_',  num2str(NumMesh) '_t1.obj' ],  ...
                                    [ average_surface_dir 'average_intracortical_surface_1_right_', num2str(NumMesh) '_t1.obj' ] });
TemplateSurf{3} = SurfStatReadSurf({[ average_surface_dir 'average_intracortical_surface_2_left_',  num2str(NumMesh) '_t1.obj' ],  ...
                                    [ average_surface_dir 'average_intracortical_surface_2_right_', num2str(NumMesh) '_t1.obj' ] });
TemplateSurf{4} = SurfStatReadSurf({[ average_surface_dir 'average_intracortical_surface_3_left_',  num2str(NumMesh) '_t1.obj' ],  ...
                                    [ average_surface_dir 'average_intracortical_surface_3_right_', num2str(NumMesh) '_t1.obj' ] });
TemplateSurf{5} = SurfStatReadSurf({[ average_surface_dir 'average_white_surface_left_',  num2str(NumMesh) '_t1.obj' ],  ...
                                    [ average_surface_dir 'average_white_surface_right_', num2str(NumMesh) '_t1.obj' ] });
TemplateSurf{6} = SurfStatReadSurf({[ average_surface_dir 'average_white_surface_1_left_',  num2str(NumMesh) '_t1.obj' ],  ...
                                    [ average_surface_dir 'average_white_surface_1_right_', num2str(NumMesh) '_t1.obj' ] });
TemplateSurf{7} = SurfStatReadSurf({[ average_surface_dir 'average_white_surface_2_left_',  num2str(NumMesh) '_t1.obj' ],  ...
                                    [ average_surface_dir 'average_white_surface_2_right_', num2str(NumMesh) '_t1.obj' ] });
TemplateSurf{8} = SurfStatReadSurf({[ average_surface_dir 'average_white_surface_3_left_',  num2str(NumMesh) '_t1.obj' ],  ...
                                    [ average_surface_dir 'average_white_surface_3_right_', num2str(NumMesh) '_t1.obj' ] });
TemplateSurf_standard = SurfStatReadSurf({'/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/surf_reg_model_left.obj', ...
                                          '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/surf_reg_model_right.obj'});
TemplateMask = SurfStatMaskCut(TemplateSurf_standard);
AAL_data = SurfStatReadData1('/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/aal_both_rsl_final.txt');
    
%% Read control database ...
%% Modality: T1, FLAIR, IR

%% Feature: RI_IntraCortical (vertexnum x 5 x numControl), 
%%          RI_SubCortical   (vertexnum x 3 x numControl), 
%%          PG_IntraCortical (vertexnum x 5 x numControl),
%%          PG_GM_WM         (vertexnum x numControl),
%%          PG_SubCortical   (vertexnum x 3 x numControl),
%%          TG_IntraCortical (vertexnum x 5 x numControl),
%%          TG_SubCortical   (vertexnum x 3 x numControl),
%%          CT, MC, SD       (vertexnum x numControl, only T1)

%% Abb: RI relative intensity
%%      PG perpendicular gradient
%%      TG tangential gradient
%%      CT cortical thickness
%%      MC mean curvature
%%      SD sulcal depth

%% T1
%% Intensity-based features
T1_RI_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
T1_PG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
T1_PG_gw_IntCortical_temp   = zeros(size(case_num_cont, 1), vertexnum);
T1_CT_midCortical_temp      = zeros(size(case_num_cont, 1), vertexnum);
T1_SD_wmCortical_temp       = zeros(size(case_num_cont, 1), vertexnum);
T1_MC_midCortical_temp      = zeros(size(case_num_cont, 1), vertexnum);
T1_RI_wCortical_temp        = zeros(size(case_num_cont, 1), vertexnum);

data_dir = '/local_raid/seokjun/01_project/IntracorticalAnalysis/03_result/';
file_postfix = [ '_quadratic_sm_5_rsl.txt' ];
for j = 1 : size(case_num_cont, 1)
    fprintf([ case_num_cont{j} ' : ' ]);
    prefix_path = [ data_dir '/' case_num_cont{j} '/measurement/' Prefix_cont '_' case_num_cont{j} ]; 
    postfix_surf = '_native_t1_';
    
    %% Intensity-based features: corrected RI, pg, tg
    for i = 1 : NumIntSurf+2    
        if(i == 1)
            basename{1}  = [ prefix_path '_gray_surface_left_'  num2str(NumMesh) postfix_surf ];
            basename{2}  = [ prefix_path '_gray_surface_right_'  num2str(NumMesh) postfix_surf ];
        elseif(i == 5)
            basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
            basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
        else
            basename{1}  = [ prefix_path '_intracortical_surface_' num2str(i-1) '_left_'  num2str(NumMesh) postfix_surf ];
            basename{2}  = [ prefix_path '_intracortical_surface_' num2str(i-1) '_right_'  num2str(NumMesh) postfix_surf ];
        end
        
        T1_RI_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI_corrected' file_postfix ], [ basename{2} 'RI_corrected' file_postfix ] } ); 
        T1_PG_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
        %T1_TG_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );        
    end
    T1_RI_wCortical_temp(j, 1:vertexnum) = mean(T1_RI_IntCortical_temp([2 3 4 5], 1:vertexnum, j), 1);
    fprintf([ 'T1: RI_corrected, PG, TG in Inc, ' ]);    
    
    %% Morphological / Intensity-based features: MC, SD, PG_GW
    basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
    basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
    T1_PG_gw_IntCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'pg_GM_WM' file_postfix ], [ basename{2} 'pg_GM_WM' file_postfix ] } );
    
    basename{1}  = [ prefix_path '_intracortical_surface_2_left_'  num2str(NumMesh) postfix_surf ];
    basename{2}  = [ prefix_path '_intracortical_surface_2_right_'  num2str(NumMesh) postfix_surf ];
    T1_CT_midCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'ct' file_postfix ], [ basename{2} 'ct' file_postfix ] } );
    
    basename{1}  = [ prefix_path '_intracortical_surface_2_left_'  num2str(NumMesh) postfix_surf ];
    basename{2}  = [ prefix_path '_intracortical_surface_2_right_'  num2str(NumMesh) postfix_surf ];
    T1_MC_midCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'mc' file_postfix ], [ basename{2} 'mc' file_postfix ] } );
    
    basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
    basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
    T1_SD_wmCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'sd2' file_postfix ], [ basename{2} 'sd2' file_postfix ] } );    
    fprintf([ 'CT, MC, SD, PG_gw in Inc\n' ]);
end

%% Mean of T1 features
T1_RI_IntCortical_mean      = mean(T1_RI_IntCortical_temp(:, 1:vertexnum, :), 3);
T1_PG_IntCortical_mean      = mean(T1_PG_IntCortical_temp(:, 1:vertexnum, :), 3);
T1_PG_gw_IntCortical_mean   = mean(T1_PG_gw_IntCortical_temp(:, 1:vertexnum), 1);
T1_CT_midCortical_mean      = mean(T1_CT_midCortical_temp(:, 1:vertexnum), 1);
T1_MC_midCortical_mean      = mean(T1_MC_midCortical_temp(:, 1:vertexnum), 1);
T1_SD_wmCortical_mean       = mean(T1_SD_wmCortical_temp(:, 1:vertexnum), 1);
T1_RI_wCortical_mean        = mean(T1_RI_wCortical_temp(:, 1:vertexnum), 1);

T1_RI_IntCortical_std       = std(T1_RI_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
T1_PG_IntCortical_std       = std(T1_PG_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
T1_PG_gw_IntCortical_std    = std(T1_PG_gw_IntCortical_temp(:, 1:vertexnum), 0, 1);
T1_CT_midCortical_std       = std(T1_CT_midCortical_temp(:, 1:vertexnum), 0, 1);
T1_MC_midCortical_std       = std(T1_MC_midCortical_temp(:, 1:vertexnum), 0, 1);
T1_SD_wmCortical_std        = std(T1_SD_wmCortical_temp(:, 1:vertexnum), 0, 1);
T1_RI_wCortical_std         = std(T1_RI_wCortical_temp(:, 1:vertexnum), 0, 1);

%% 1st mean and SD examination ...
IntMask = repmat(TemplateMask, [5, 1]);
SubMask = repmat(TemplateMask, [3, 1]);

% figure; CSFSurfStatView(T1_RI_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'RI mean'); CSFSurfStatViewColLim([-200 200]);
% figure; CSFSurfStatView(T1_RI_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'RI std');  CSFSurfStatViewColLim([0 50]);
% figure; CSFSurfStatView(T1_PG_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'PG mean'); CSFSurfStatViewColLim([0 300]);
% figure; CSFSurfStatView(T1_PG_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'PG std');  CSFSurfStatViewColLim([0 150]);
% figure; CSFSurfStatView(T1_TG_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'TG mean'); CSFSurfStatViewColLim([0 50]);
% figure; CSFSurfStatView(T1_TG_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'TG std');  CSFSurfStatViewColLim([0 20]);
% figure; CSFSurfStatView(T1_PG_gw_IntCortical_mean.*TemplateMask, TemplateSurf{3}, 'MF', 'PG mean'); CSFSurfStatViewColLim([0 300]);
% figure; CSFSurfStatView(T1_PG_gw_IntCortical_std.*TemplateMask,  TemplateSurf{3}, 'MF', 'PG std');  CSFSurfStatViewColLim([0 150]);
% figure; CSFSurfStatView(T1_CT_midCortical_mean.*TemplateMask, TemplateSurf{3}, 'MF', 'CT mean'); CSFSurfStatViewColLim([0 5]);
% figure; CSFSurfStatView(T1_CT_midCortical_std.*TemplateMask,  TemplateSurf{3}, 'MF', 'CT std');  CSFSurfStatViewColLim([0 1]);
% figure; CSFSurfStatView(T1_MC_midCortical_mean.*TemplateMask, TemplateSurf{3}, 'MF', 'MC mean'); CSFSurfStatViewColLim([0 0.3]);
% figure; CSFSurfStatView(T1_MC_midCortical_std.*TemplateMask,  TemplateSurf{3}, 'MF', 'MC std');  CSFSurfStatViewColLim([0 0.15]);
% figure; CSFSurfStatView(T1_SD_wmCortical_mean.*TemplateMask, TemplateSurf{3},  'MF', 'SD mean'); CSFSurfStatViewColLim([0 25]);
% figure; CSFSurfStatView(T1_SD_wmCortical_std.*TemplateMask,  TemplateSurf{3},  'MF', 'SD std');  CSFSurfStatViewColLim([0 10]);
% 
% figure; CSFSurfStatView(T1_RI_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'RI mean'); CSFSurfStatViewColLim([-100 100]);
% figure; CSFSurfStatView(T1_RI_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'RI std');  CSFSurfStatViewColLim([0 100]);
% figure; CSFSurfStatView(T1_PG_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'PG mean'); CSFSurfStatViewColLim([0 100]);
% figure; CSFSurfStatView(T1_PG_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'PG std');  CSFSurfStatViewColLim([0 50]);
% figure; CSFSurfStatView(T1_TG_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'TG mean'); CSFSurfStatViewColLim([0 40]);
% figure; CSFSurfStatView(T1_TG_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'TG std');  CSFSurfStatViewColLim([0 15]);

%% Outlier detection
sid = 1;
temp_mean = T1_RI_IntCortical_mean(sid, :);
temp_std  = T1_RI_IntCortical_std(sid, :);
temp_mean_new = temp_mean;
temp_std_new = temp_std;

vid_max = 100;
sub_idx = zeros(size(case_num_cont, 1), 19);

%% RI outliers
for j = 1 : 5
    sid = j;
    for i = 1 : vid_max
        disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
        vid = i;
        [b,idx,outliers] = deleteoutliers(T1_RI_IntCortical_temp(sid, vid, :), 0.1);
        sub_idx(idx, j) = sub_idx(idx, j) + 1;
    end    
end

%% PG outliers
for j = 6 : 10
    sid = j-5;
    for i = 1 : vid_max
        disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
        vid = i;
        [b,idx,outliers] = deleteoutliers(T1_PG_IntCortical_temp(sid, vid, :), 0.1);
        sub_idx(idx, j) = sub_idx(idx, j) + 1;
    end    
end


%% CT,MC,SD,PG_gw outliers
j = 16;
for i = 1 : vid_max
    disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
    vid = i;
    [b,idx,outliers] = deleteoutliers(T1_CT_midCortical_temp(:, vid), 0.1);
    sub_idx(idx, j) = sub_idx(idx, j) + 1;
end
j = 17;
for i = 1 : vid_max
    disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
    vid = i;
    [b,idx,outliers] = deleteoutliers(T1_PG_gw_IntCortical_temp(:, vid), 0.1);
    sub_idx(idx, j) = sub_idx(idx, j) + 1;
end
j = 18;
for i = 1 : vid_max
    disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
    vid = i;
    [b,idx,outliers] = deleteoutliers(T1_MC_midCortical_temp(:, vid), 0.1);
    sub_idx(idx, j) = sub_idx(idx, j) + 1;
end
j = 19;
for i = 1 : vid_max
    disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
    vid = i;
    [b,idx,outliers] = deleteoutliers(T1_SD_wmCortical_temp(:, vid), 0.1);
    sub_idx(idx, j) = sub_idx(idx, j) + 1;
end

[outlier_score, idx_outlier ] = sort(sum(sub_idx, 2)/1900, 'descend');
strvcat(case_num_cont{idx_outlier})

%% z-score
T1_RI_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
T1_PG_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
T1_CT_midCortical_z         = zeros(size(case_num_cont, 1), vertexnum);
T1_MC_midCortical_z         = zeros(size(case_num_cont, 1), vertexnum);
T1_SD_wmCortical_z          = zeros(size(case_num_cont, 1), vertexnum);
T1_PG_gw_IntCortical_z      = zeros(size(case_num_cont, 1), vertexnum);
T1_RI_wCortical_z           = zeros(size(case_num_cont, 1), vertexnum);

%% T1: calculate z-score of patients w.r.t control's distribution
for j = 1 : size(case_num_cont, 1)
    fprintf([ case_num_cont{j} ' : ' ]);
    
    %% Intensity-based features: corrected RI, pg, tg      
    T1_RI_IntCortical_z(:, 1:vertexnum, j) = (T1_RI_IntCortical_temp(:, 1:vertexnum, j) - T1_RI_IntCortical_mean) ./ T1_RI_IntCortical_std;
    T1_PG_IntCortical_z(:, 1:vertexnum, j) = (T1_PG_IntCortical_temp(:, 1:vertexnum, j) - T1_PG_IntCortical_mean) ./ T1_PG_IntCortical_std;
    T1_RI_IntCortical_z(isnan(T1_RI_IntCortical_z)) = 0; T1_RI_IntCortical_z(isinf(T1_RI_IntCortical_z)) = 0;
    T1_PG_IntCortical_z(isnan(T1_PG_IntCortical_z)) = 0; T1_PG_IntCortical_z(isinf(T1_PG_IntCortical_z)) = 0;
    fprintf([ 'z-score of T1: RI_corrected, PG in Inc, ' ]);
    
    %% Morphological / Intensity-based features: MC, SD, PG_GW
    T1_PG_gw_IntCortical_z(j, 1:vertexnum) = (T1_PG_gw_IntCortical_temp(j, 1:vertexnum) - T1_PG_gw_IntCortical_mean) ./ T1_PG_gw_IntCortical_std;
    T1_CT_midCortical_z(j, 1:vertexnum)    = (T1_CT_midCortical_temp(j, 1:vertexnum) - T1_CT_midCortical_mean) ./ T1_CT_midCortical_std;
    T1_MC_midCortical_z(j, 1:vertexnum)    = (T1_MC_midCortical_temp(j, 1:vertexnum) - T1_MC_midCortical_mean) ./ T1_MC_midCortical_std;
    T1_SD_wmCortical_z(j, 1:vertexnum)     = (T1_SD_wmCortical_temp(j, 1:vertexnum) - T1_SD_wmCortical_mean) ./ T1_SD_wmCortical_std;
    T1_RI_wCortical_z(j, 1:vertexnum)     = (T1_RI_wCortical_temp(j, 1:vertexnum) - T1_RI_wCortical_mean) ./ T1_RI_wCortical_std;
    
    T1_PG_gw_IntCortical_z(isnan(T1_PG_gw_IntCortical_z)) = 0; T1_PG_gw_IntCortical_z(isinf(T1_PG_gw_IntCortical_z)) = 0;
    T1_CT_midCortical_z(isnan(T1_CT_midCortical_z)) = 0; T1_CT_midCortical_z(isinf(T1_CT_midCortical_z)) = 0;
    T1_MC_midCortical_z(isnan(T1_MC_midCortical_z)) = 0; T1_MC_midCortical_z(isinf(T1_MC_midCortical_z)) = 0;   
    T1_SD_wmCortical_z(isnan(T1_SD_wmCortical_z)) = 0; T1_SD_wmCortical_z(isinf(T1_SD_wmCortical_z)) = 0;     
    T1_RI_wCortical_z(isnan(T1_RI_wCortical_z)) = 0; T1_RI_wCortical_z(isinf(T1_RI_wCortical_z)) = 0;     
    fprintf([ 'CT, MC, SD, PG_gw in Inc\n' ]);    
end

% 1st mean and SD examination ...
IntMask = repmat(TemplateMask, [5, 1]);
SubMask = repmat(TemplateMask, [3, 1]);
OUTPATH_temp = '/local_raid/seokjun/01_project/IntracorticalAnalysis/03_result/zscoremap_frontiers/t1_1.5T';
for j = 1 : size(case_num_cont, 1)
    idx = find(strcmp(case_num_cont, case_num_cont{j}));
    f = figure; CSFSurfStatView(T1_RI_wCortical_z(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3}, 'MF', 'RI z-score');         CSFSurfStatViewColLim([-5 5]);
    exportfigbo(f,[OUTPATH_temp '/' Prefix_cont '_' case_num_cont{j} '_z-score_T1_Intra_RIw.png' ], 'png', 6); close(f);    
%     f = figure; CSFSurfStatView(T1_RI_IntCortical_z(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'RI z-score');         CSFSurfStatViewColLim([-5 5]);
%     exportfigbo(f,[OUTPATH_temp '/' Prefix_cont '_' case_num_cont{j} '_z-score_T1_Intra_RI.png' ], 'png', 6); close(f);
%     f = figure; CSFSurfStatView(T1_PG_IntCortical_z(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'PG z-score');         CSFSurfStatViewColLim([-5 5]);
%     exportfigbo(f,[OUTPATH_temp '/' Prefix_cont '_' case_num_cont{j} '_z-score_T1_Intra_PG.png' ], 'png', 6); close(f);
%     f = figure; CSFSurfStatView(T1_PG_gw_IntCortical_z(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3}, 'MF', 'PG z-score');     CSFSurfStatViewColLim([-5 5]);
%     exportfigbo(f,[OUTPATH_temp '/' Prefix_cont '_' case_num_cont{j} '_z-score_T1_GM-WM_PG.png' ], 'png', 6); close(f);
%     f = figure; CSFSurfStatView(T1_CT_midCortical_z(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3}, 'MF', 'CT z-score');          CSFSurfStatViewColLim([-5 5]);
%     exportfigbo(f,[OUTPATH_temp '/' Prefix_cont '_' case_num_cont{j} '_z-score_T1_Mid_CT.png' ], 'png', 6); close(f);
%     f = figure; CSFSurfStatView(T1_MC_midCortical_z(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3}, 'MF', 'MC z-score');          CSFSurfStatViewColLim([-5 5]);
%     exportfigbo(f,[OUTPATH_temp '/' Prefix_cont '_' case_num_cont{j} '_z-score_T1_Mid_MC.png' ], 'png', 6); close(f);
%     f = figure; CSFSurfStatView(T1_SD_wmCortical_z(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3},  'MF', 'SD z-score');          CSFSurfStatViewColLim([-5 5]);
%     exportfigbo(f,[OUTPATH_temp '/' Prefix_cont '_' case_num_cont{j} '_z-score_T1_WM_SD.png' ], 'png', 6); close(f);

end

OUTPATH = '/local_raid/seokjun/01_project/IntracorticalAnalysis/03_result/zscoremap_frontiers/t1_1.5T/';
post_fix = [ '_' Parametric '_sm_' num2str(Kernel) '_cont.png' ];
z_score_database = zeros(size(case_num_cont, 2), vertexnum, 5);
f = figure; 
for i = 1 : size(case_num_cont, 1)
    
    z_score_database(i, :, 1) = T1_RI_wCortical_z(i, :);
    z_score_database(i, :, 2) = -T1_PG_gw_IntCortical_z(i, :);
    z_score_database(i, :, 3) = T1_CT_midCortical_z(i, :);
    z_score_database(i, :, 4) = T1_SD_wmCortical_z(i, :);
    z_score_database(i, :, 5) = T1_MC_midCortical_z(i, :);   
    
    SurfStatView(T1_RI_wCortical_z(i, :)-T1_PG_gw_IntCortical_z(i, :)+T1_CT_midCortical_z(i, :)+T1_SD_wmCortical_z(i, :)+T1_MC_midCortical_z(i, :), TemplateSurf{3}, case_num_cont{i});
    SurfStatColLim([-10 10]);
    exportfigbo(f,[OUTPATH '/TLE_' case_num_cont{i} '_z-score_total.png'], 'png', 6); 
end
close(f);