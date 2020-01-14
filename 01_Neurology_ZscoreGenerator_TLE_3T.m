clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OutDir                    = '/local_raid/seokjun/01_project/IntracorticalAnalysis/03_result/'
Group_pat                 = 'control'
Prefix_pat                = 'TLE'
Cases_cont                = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/Demographic_data_control.txt'
Group_cont                = 'tle'
Prefix_cont               = 'TLE'
Cases_pat                 = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/01_analysis/Demographic_data_TLE.txt'

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
age_cont = demo{2};
gender_cont = demo{3};
fclose(fid);

fid = fopen(Cases_pat);
demo = textscan(fid, '%s%f%s%d%s%d%d', 'Delimiter', ',', 'CollectOutput', 1);
case_num_pat = demo{1};
age_pat = demo{2};
gender_pat = demo{3};
fclose(fid);


mean_age_cont = mean(age_cont);
std_age_cont = std(age_cont, 0);

mean_age_pat = mean(age_pat);
std_age_pat = std(age_pat, 0);

[h,p,ci,stats] = ttest2(age_cont,age_pat);

pat_case  = {  '0361_1', '0363_1', '0364_1', '0369_1', '0370_1', '0372_1', '0373_1', '0375_1', '0379_1', '0380_1', '0390_1', '0392_1', ...
               '0394_1', '0395_1', '0396_1', '0397_1', '0401_1', '0404_1', '0415_1', '0420_1', '0422_1', '0423_1', '0427_2', '0428_1' };

cont_case = { '301_1', '303_1', '304_1', '305_1', '306_1', '307_1', '308_1', '309_1', '310_1', '311_1', '314_1', '316_1', '317_1', '321_1', '323_1', '329_1', '332_1', ...
              '324_1', '333_1', '334_1', '335_1', '337_1', '338_1', '339_1', '340_1'  };

[C, pat_idx, ib]  = intersect(case_num_pat, pat_case);
[C, cont_idx, ib] = intersect(case_num_cont, cont_case);

mean_age_cont_new = mean(age_cont(cont_idx));
std_age_cont_new = std(age_cont(cont_idx), 0);

mean_age_pat_new = mean(age_pat(pat_idx));
std_age_pat_new = std(age_pat(pat_idx), 0);

%% Read control database ...
%% RI_IntraCortical, RI_SubCortical, PG_IntraCortical, PG_SubCortical, TG_IntraCortical, TG_SubCortical, CT, MC, SD
%% RI_idx_exp_outlier, PG_idx_exp_outlier, TG_idx_exp_outlier, CT_idx_exp_outlier, MC_idx_exp_outlier, SD_idx_exp_outlier
load(['/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/average_featuremaps/Control_' ...
       img_contrast{1}, '_feature_set_' num2str(NumMesh), '_' SamplingSpace '_sm_' num2str(Kernel) '_' Parametric '.mat']);
   
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

RI_IntraCortical_temp = [];
RI_SubCortical_temp = [];
PG_IntraCortical_temp = [];
PG_SubCortical_temp = [];
TG_IntraCortical_temp = [];
TG_SubCortical_temp = [];

RI_IntraCortical_z = [];
RI_SubCortical_z = [];
PG_IntraCortical_z = [];
PG_SubCortical_z = [];
TG_IntraCortical_z = [];
TG_SubCortical_z = [];

RI_IntraCortical_mean = [];
RI_SubCortical_mean = [];
PG_IntraCortical_mean = [];
PG_SubCortical_mean = [];
TG_IntraCortical_mean = [];
TG_SubCortical_mean = [];
RI_IntraCortical_std = [];
RI_SubCortical_std = [];
PG_IntraCortical_std = [];
PG_SubCortical_std = [];
TG_IntraCortical_std = [];
TG_SubCortical_std = [];

RI_wCortical_z = zeros(size(pat_idx', 1), vertexnum);
RI_wSubCortical_z = zeros(size(pat_idx', 1), vertexnum);
RI_wCortical_mean = zeros(size(cont_idx', 1), vertexnum);
RI_wCortical_mean_temp = zeros(size(cont_idx', 1), vertexnum);
RI_wSubCortical_mean = zeros(size(cont_idx', 1), vertexnum);

s_temp = [ 2 3 4 5 ];
for i = 1 : NumIntSurf+2
    RI_IntraCortical_temp{i} = zeros(size(pat_idx', 1), vertexnum);
    PG_IntraCortical_temp{i} = zeros(size(pat_idx', 1), vertexnum);
    TG_IntraCortical_temp{i} = zeros(size(pat_idx', 1), vertexnum);    
    RI_IntraCortical_z{i} = zeros(size(pat_idx', 1), vertexnum);
    PG_IntraCortical_z{i} = zeros(size(pat_idx', 1), vertexnum);
    TG_IntraCortical_z{i} = zeros(size(pat_idx', 1), vertexnum);
    
    RI_IntraCortical_mean{i} = mean(RI_IntraCortical{i}(cont_idx, :), 1);
    PG_IntraCortical_mean{i} = mean(PG_IntraCortical{i}(cont_idx, :), 1);
    TG_IntraCortical_mean{i} = mean(TG_IntraCortical{i}(cont_idx, :), 1);    
    RI_IntraCortical_std{i} = std(RI_IntraCortical{i}(cont_idx, :), 0, 1);
    PG_IntraCortical_std{i} = std(PG_IntraCortical{i}(cont_idx, :), 0, 1);
    TG_IntraCortical_std{i} = std(TG_IntraCortical{i}(cont_idx, :), 0, 1); 

    if( ismember(i, s_temp) )
        RI_wCortical_mean_temp = RI_wCortical_mean_temp + RI_IntraCortical{i}(cont_idx, :);
    end    
end

for i = 1 : NumSubSurf
    RI_SubCortical_temp{i} = zeros(size(pat_idx', 1), vertexnum);
    PG_SubCortical_temp{i} = zeros(size(pat_idx', 1), vertexnum);
    TG_SubCortical_temp{i} = zeros(size(pat_idx', 1), vertexnum);
    RI_SubCortical_z{i} = zeros(size(pat_idx', 1), vertexnum);
    PG_SubCortical_z{i} = zeros(size(pat_idx', 1), vertexnum);
    TG_SubCortical_z{i} = zeros(size(pat_idx', 1), vertexnum);
    
    RI_SubCortical_mean{i} = mean(RI_SubCortical{i}(cont_idx, :), 1);
    PG_SubCortical_mean{i} = mean(PG_SubCortical{i}(cont_idx, :), 1);
    TG_SubCortical_mean{i} = mean(TG_SubCortical{i}(cont_idx, :), 1);        
    RI_SubCortical_std{i} = std(RI_SubCortical{i}(cont_idx, :), 0, 1);
    PG_SubCortical_std{i} = std(PG_SubCortical{i}(cont_idx, :), 0, 1);
    TG_SubCortical_std{i} = std(TG_SubCortical{i}(cont_idx, :), 0, 1);     
    
    RI_wSubCortical_mean = RI_wSubCortical_mean + RI_SubCortical{i}(cont_idx, :);
end

PG2_temp = zeros(size(pat_idx', 1), vertexnum);
PG2_z = zeros(size(pat_idx', 1), vertexnum);
PG2_mean = mean(PG2(cont_idx, :));
PG2_std  = std(PG2(cont_idx, :), 0);

RI_wCortical_std = std(RI_wCortical_mean_temp, 0);
RI_wCortical_mean = mean(RI_wCortical_mean_temp / length(s_temp), 1);

RI_wSubCortical_std = std(RI_wSubCortical_mean, 0);
RI_wSubCortical_mean = mean(RI_wSubCortical_mean / (NumSubSurf), 1);

if(strcmp(img_contrast{1}, 't1'))    
    CT_temp = zeros(size(pat_idx', 1), vertexnum);
    MC_temp = zeros(size(pat_idx', 1), vertexnum);
    SD_temp = zeros(size(pat_idx', 1), vertexnum);
    SD2_temp = zeros(size(pat_idx', 1), vertexnum);    
    
    CT_z = zeros(size(pat_idx', 1), vertexnum);
    MC_z = zeros(size(pat_idx', 1), vertexnum);
    SD_z = zeros(size(pat_idx', 1), vertexnum);
    SD2_z = zeros(size(pat_idx', 1), vertexnum);
    
    CT_mean = mean(CT(cont_idx,:), 1);
    MC_mean = mean(abs(MC(cont_idx,:)), 1);
    SD_mean = mean(SD(cont_idx,:), 1);
    SD2_mean = mean(SD2(cont_idx,:), 1);
    CT_std = std(CT(cont_idx,:), 0, 1);
    MC_std = std(abs(MC(cont_idx,:)), 0, 1);
    SD_std = std(SD(cont_idx,:), 0, 1);    
    SD2_std = std(SD2(cont_idx,:), 0);
end

clear('RI_IntraCortical', 'PG_IntraCortical', 'TG_IntraCortical', 'RI_SubCortical', 'PG_SubCortical', 'TG_SubCortical');

for i = 1 : size(pat_idx', 1)
    for l_r = 1 : size(Left_Right, 2)
        %% RI, PG, TG
        Postfix1 = [ Left_Right{1, l_r} '_' num2str(NumMesh) '_' SamplingSpace '_' img_contrast{1} ];
        Postfix2 = [ Parametric '_sm_' num2str(Kernel) '_rsl.txt' ];
        for s = 1 : NumIntSurf + 2
            switch s
                case 1                        
                    postfix_temp = [ OutDir '/' case_num_pat{pat_idx(i)} '/measurement/' Prefix_pat '_' case_num_pat{pat_idx(i)} '_gray_surface_' Postfix1 ]
                    RI_IntraCortical_temp{s}(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_RI_corrected_' Postfix2 ] });
                    PG_IntraCortical_temp{s}(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_pg_' Postfix2 ] });
                    
                case {2, 3, 4}
                    postfix_temp = [ OutDir '/' case_num_pat{pat_idx(i)} '/measurement/' Prefix_pat '_' case_num_pat{pat_idx(i)} '_intracortical_surface_' num2str(s-1) '_' Postfix1 ]
                    RI_IntraCortical_temp{s}(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_RI_corrected_' Postfix2 ] });
                    PG_IntraCortical_temp{s}(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_pg_' Postfix2 ] });
                case 5
                    postfix_temp = [ OutDir '/' case_num_pat{pat_idx(i)} '/measurement/' Prefix_pat '_' case_num_pat{pat_idx(i)} '_white_surface_' Postfix1 ]
                    RI_IntraCortical_temp{s}(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_RI_corrected_' Postfix2 ] });
                    PG_IntraCortical_temp{s}(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_pg_' Postfix2 ] });
            end
        end

        for s = 1 : NumSubSurf
            postfix_temp = [ OutDir '/' case_num_pat{pat_idx(i)} '/measurement/' Prefix_pat '_' case_num_pat{pat_idx(i)} '_intracortical_surface_' num2str(s) '_' Postfix1 ]
            RI_SubCortical_temp{s}(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_RI_' Postfix2 ] });
            PG_SubCortical_temp{s}(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_pg_' Postfix2 ] });
        end      
        
        postfix_temp = [ OutDir '/' case_num_pat{pat_idx(i)} '/measurement/' Prefix_pat '_' case_num_pat{pat_idx(i)} '_white_surface_' Postfix1 ]
        PG2_temp(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_pg_GM_WM_' Postfix2 ] });

        if(strcmp(img_contrast{1}, 't1'))
            postfix_temp = [ OutDir '/' case_num_pat{pat_idx(i)} '/measurement/' Prefix_pat '_' case_num_pat{pat_idx(i)} '_intracortical_surface_2_' Postfix1 ]
            CT_temp(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_ct_' Postfix2 ] });                
            MC_temp(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_mc_' Postfix2 ] });
            postfix_temp = [ OutDir '/' case_num_pat{pat_idx(i)} '/measurement/' Prefix_pat '_' case_num_pat{pat_idx(i)} '_white_surface_' Postfix1 ];
            SD2_temp(i, Left_Right{2, l_r}) = SurfStatReadData({ [ postfix_temp '_sd2_' Postfix2 ] });
        end
    end
end  

for i = 1 : size(pat_idx', 1)
    RI_wCortical_temp = zeros(1, vertexnum);
    RI_wSubCortical_temp = zeros(1, vertexnum);    
    temp = 0;
    for l_r = 1 : size(Left_Right, 2)
        %% RI, PG, TG
        for s = 1 : NumIntSurf + 2
            s
            RI_IntraCortical_z{s}(i, Left_Right{2, l_r}) = (RI_IntraCortical_temp{s}(i, Left_Right{2, l_r}) - RI_IntraCortical_mean{s}(Left_Right{2, l_r}))./RI_IntraCortical_std{s}(Left_Right{2, l_r});
            PG_IntraCortical_z{s}(i, Left_Right{2, l_r}) = (PG_IntraCortical_temp{s}(i, Left_Right{2, l_r}) - PG_IntraCortical_mean{s}(Left_Right{2, l_r}))./PG_IntraCortical_std{s}(Left_Right{2, l_r});
            
            if(ismember(s, s_temp))
                RI_wCortical_temp(Left_Right{2, l_r}) = RI_wCortical_temp(Left_Right{2, l_r}) + RI_IntraCortical_temp{s}(i, Left_Right{2, l_r});                
                temp = temp + 1;
            end
        end

        for s = 1 : NumSubSurf
            s
            RI_SubCortical_z{s}(i, Left_Right{2, l_r}) = (RI_SubCortical_temp{s}(i, Left_Right{2, l_r}) - RI_SubCortical_mean{s}(Left_Right{2, l_r}))./RI_SubCortical_std{s}(Left_Right{2, l_r});
            PG_SubCortical_z{s}(i, Left_Right{2, l_r}) = (PG_SubCortical_temp{s}(i, Left_Right{2, l_r}) - PG_SubCortical_mean{s}(Left_Right{2, l_r}))./PG_SubCortical_std{s}(Left_Right{2, l_r});
            
            RI_wSubCortical_temp(Left_Right{2, l_r}) = RI_wSubCortical_temp(Left_Right{2, l_r}) + RI_SubCortical_temp{s}(i, Left_Right{2, l_r});
        end

        if(strcmp(img_contrast{1}, 't1'))
            CT_z(i, Left_Right{2, l_r}) = (CT_temp(i, Left_Right{2, l_r}) - CT_mean(Left_Right{2, l_r}))./CT_std(Left_Right{2, l_r});
            MC_z(i, Left_Right{2, l_r}) = (abs(MC_temp(i, Left_Right{2, l_r})) - MC_mean(Left_Right{2, l_r}))./MC_std(Left_Right{2, l_r});
            SD2_z(i, Left_Right{2, l_r}) = (SD2_temp(i, Left_Right{2, l_r}) - SD2_mean(Left_Right{2, l_r}))./SD2_std(Left_Right{2, l_r});
        end
    end
    
    RI_wCortical_temp = RI_wCortical_temp / (temp / size(Left_Right, 2));
    RI_wSubCortical_temp = RI_wSubCortical_temp / NumSubSurf;
    
    RI_wCortical_z(i, :) = (RI_wCortical_temp - RI_wCortical_mean)./RI_wCortical_std;
    RI_wSubCortical_z(i, :) = (RI_wSubCortical_temp - RI_wSubCortical_mean)./RI_wSubCortical_std;    
    
    PG2_z(i, :) = (PG2_temp(i, :) - PG2_mean)./PG2_std;
end

if(visualization)
    vis_basedir = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/average_surfaces/';    
    for s = 1 : NumIntSurf + 2
        switch s
            case 1                        
                Postfix1 = [ 'average_gray_surface_' ];
                Postfix2 = [ '_' num2str(NumMesh) '_' img_contrast{1} ];
                SurfIntra{s} = SurfStatReadSurf({[ vis_basedir '/' Postfix1 'left' Postfix2 '.obj' ], [ vis_basedir '/' Postfix1 'right' Postfix2 '.obj'  ]});
            case {2, 3, 4}
                Postfix1 = [ 'average_intracortical_surface_' num2str(s-1) '_' ];
                Postfix2 = [ '_' num2str(NumMesh) '_' img_contrast{1} ];
                SurfIntra{s} = SurfStatReadSurf({[ vis_basedir '/' Postfix1 'left' Postfix2 '.obj' ], [ vis_basedir '/' Postfix1 'right' Postfix2 '.obj'  ]});
            case 5
                Postfix1 = [ 'average_white_surface_' ];
                Postfix2 = [ '_' num2str(NumMesh) '_' img_contrast{1} ];
                SurfIntra{s} = SurfStatReadSurf({[ vis_basedir '/' Postfix1 'left' Postfix2 '.obj' ], [ vis_basedir '/' Postfix1 'right' Postfix2 '.obj'  ]});
        end
    end

    for s = 1 : NumSubSurf
        Postfix1 = [ 'average_white_surface_' num2str(s) '_'];
        Postfix2 = [ '_' num2str(NumMesh) '_' img_contrast{1} ];
        SurfSub{s} = SurfStatReadSurf({[ vis_basedir '/' Postfix1 'left' Postfix2 '.obj' ], [ vis_basedir '/' Postfix1 'right' Postfix2 '.obj'  ]});
    end
end

%% RI_wCortical_mean: add RI of 2~5th surfaces at each individiual and make the average ...
%% RI_IntraCortical_mean: make the average of RI at each surface across subjects ...
figure; SurfStatView(RI_wCortical_mean, SurfIntra{3}); SurfStatColLim([-100 50]);
figure; SurfStatView(RI_wSubCortical_mean, SurfSub{2}); SurfStatColLim([-100 50]);
figure; SurfStatView(PG2_mean, SurfIntra{3}); SurfStatColLim([0 150]);
figure; SurfStatView(RI_wCortical_std, SurfIntra{3}); SurfStatColLim([0 150]);
figure; SurfStatView(RI_wSubCortical_std, SurfSub{2}); SurfStatColLim([0 150]);
figure; SurfStatView(PG2_std, SurfIntra{3}); SurfStatColLim([0 150]);

OUTPATH = '/data/noel/noel6/seokjun/01_project/IntracorticalAnalysis/03_result/zscoremap_frontiers/';
post_fix = [ '_' Parametric '_sm_' num2str(Kernel) '_pat.png' ];
z_score_database = zeros(size(pat_idx, 2), vertexnum, 5);
for i = 1 : size(pat_idx', 1)
    RI_intra_z = [];
    PG_intra_z = [];
    TG_intra_z = [];
    for s = 1 : NumIntSurf + 2
        RI_intra_z_temp = [];
        PG_intra_z_temp = [];
        TG_intra_z_temp = [];
        for l_r = 1 : size(Left_Right, 2)
            RI_intra_z_temp = [ RI_intra_z_temp RI_IntraCortical_z{s}(i, Left_Right{2, l_r}) ];
            PG_intra_z_temp = [ PG_intra_z_temp PG_IntraCortical_z{s}(i, Left_Right{2, l_r}) ];
            TG_intra_z_temp = [ TG_intra_z_temp TG_IntraCortical_z{s}(i, Left_Right{2, l_r}) ];
        end
        RI_intra_z = [ RI_intra_z; RI_intra_z_temp ];
        PG_intra_z = [ PG_intra_z; PG_intra_z_temp ];
        TG_intra_z = [ TG_intra_z; TG_intra_z_temp ];
    end

    RI_sub_z = [];
    PG_sub_z = [];
    TG_sub_z = [];
    for s = 1 : NumSubSurf
        RI_sub_z_temp = [];
        PG_sub_z_temp = [];
        TG_sub_z_temp = [];
        for l_r = 1 : size(Left_Right, 2)
            RI_sub_z_temp = [ RI_sub_z_temp RI_SubCortical_z{s}(i, Left_Right{2, l_r}) ];
            PG_sub_z_temp = [ PG_sub_z_temp PG_SubCortical_z{s}(i, Left_Right{2, l_r}) ];
            TG_sub_z_temp = [ TG_sub_z_temp TG_SubCortical_z{s}(i, Left_Right{2, l_r}) ];
        end
        RI_sub_z = [ RI_sub_z; RI_sub_z_temp ];
        PG_sub_z = [ PG_sub_z; PG_sub_z_temp ];
        TG_sub_z = [ TG_sub_z; TG_sub_z_temp ];
    end

    if(strcmp(img_contrast{1}, 't1'))
        CT_z_temp = [];
        MC_z_temp = [];
        SD_z_temp = [];
        for l_r = 1 : size(Left_Right, 2)
            CT_z_temp = [ CT_z_temp CT_z(i, Left_Right{2, l_r}) ];
            MC_z_temp = [ MC_z_temp MC_z(i, Left_Right{2, l_r}) ];
            SD_z_temp = [ SD_z_temp SD_z(i, Left_Right{2, l_r}) ];
        end
    end
    SD2_z(i, isnan(SD2_z(i, :))) = 0;
    SD2_z(i, isinf(SD2_z(i, :))) = 0;
    
    z_score_database(i, :, 1) = RI_wCortical_z(i, :);
    z_score_database(i, :, 2) = -PG2_z(i, :);
    z_score_database(i, :, 3) = CT_z_temp;
    z_score_database(i, :, 4) = SD2_z(i, :);
    z_score_database(i, :, 5) = MC_z_temp;   
    
%     f = figure; SurfStatView(RI_wCortical_z(i, :)-PG2_z(i, :)+CT_z_temp+SD2_z(i, :)+MC_z_temp, SurfIntra{3}, case_num_pat{pat_idx(i)}); cameramenu; axis off; material dull;
%     SurfStatColLim([-10 10]);
%     exportfigbo(f,[OUTPATH '/' img_contrast{1} '/' case_num_pat{pat_idx(i)} '_zscore_Total1.png'], 'png', 6); close(f);
% 
%     f = figure; SurfStatView(mean(RI_intra_z)-PG2_z(i, :)+CT_z_temp+SD2_z(i, :)+MC_z_temp, SurfIntra{3}, case_num_pat{pat_idx(i)}); cameramenu; axis off; material dull;
%     SurfStatColLim([-10 10]);
%     exportfigbo(f,[OUTPATH '/' img_contrast{1} '/' case_num_pat{pat_idx(i)} '_zscore_Total2.png'], 'png', 6); close(f);
%    
%     f = figure; CSFSurfStatViewDataINC(transpose(squeeze(z_score_database(i, :, :))), { SurfIntra{3}, SurfIntra{3}, SurfIntra{3}, SurfIntra{3}, SurfIntra{3} }, case_num_pat{pat_idx(i)});
%     CSFSurfStatViewColLim([-5 5]);
%     exportfigbo(f,[OUTPATH '/' img_contrast{1} '/' case_num_pat{pat_idx(i)} '_zscore_features.png'], 'png', 6); close(f);
%     
    f = figure; SurfStatView(MC_z_temp, SurfIntra{3});
    SurfStatColLim([-5 5]);    
    exportfigbo(f,[OUTPATH '/' img_contrast{1} '/' case_num_pat{pat_idx(i)} '_MC_zscore_features.png'], 'png', 6); close(f);
end