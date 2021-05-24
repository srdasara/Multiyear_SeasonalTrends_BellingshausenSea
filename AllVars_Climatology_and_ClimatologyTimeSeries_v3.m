% Here I plot out climatologies both spatially and temporally. Spatial
% climatologies will show seasonal patterns, and temporal climatologies
% will show monthly patterns. All variables of interest are analyzed
% within, including MAOD, chl-a, Ice, and wind speed. 

%%


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication
load('Total_timetable_amsrmf.mat')
load('Total_timetable_Depol_Ratio.mat')
load('Total_timetable_CMOD_Surface.mat')


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/New_mat_files_CHL_A_DATA_MONTHLY

load('Aqua.mat') 
load('Aqua_Longitude_Subset_BellingshausenSea.mat') 
load('Aqua_Latitude_Subset_BellingshausenSea.mat')
load('Master_chlor_a_monthly_full_res.mat') 

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication

t1 = datetime(2006,06,01);
t2 = datetime(2018,12,31);
times = t1:calmonths(1):t2; 
clear t1 t2
OPTIONS.MaxIter = 1000;

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Revisions_Figures

%%

Aqua_Lat = (linspace(Latitude_subset(1), Latitude_subset(end), 15)) ; 
Aqua_Lon = (linspace(Longitude_subset(1), Longitude_subset(end), 46)); 
Aqua_Lat = flip(Aqua_Lat); % this is for plotting purposes because pcolor flips image on y axis. 

%%

nAqua_Lat = Latitude_subset(1) : -1 : Latitude_subset(end) ; 
nAqua_Lon = Longitude_subset(1) : 1 : Longitude_subset(end); 

nAqua_Lat = flip(nAqua_Lat) ; 

%%
% I can index out the chlorophyll master array which is 1104 x 360 x 151
% into separate seasons for Winter, Spring, Summer, & Fall 

% After that, the climatology would just be the average of those seasons. 

% For Winter: 

t1 = datetime(2006,06,01);
t2 = datetime(2018,12,31);
times = t1:calmonths(1):t2; 

times = times';

[winter_x, winter_y] = find(times.Month >= 6 & times.Month <= 8);
% Can check if this worked by typing 'times(winter_x)' in command window
times(winter_x)
Chl_a_winter = Master_chl_a; 
Chl_a_winter = Chl_a_winter(:,:,winter_x);

% For Spring: 
[spring_x, spring_y] = find(times.Month >= 9 & times.Month <= 11);
times(spring_x)
Chl_a_spring = Master_chl_a; 
Chl_a_spring = Chl_a_spring(:,:, spring_x); 


% For Summer: 

[summer_x, summer_y] = find(times.Month >= 1 & times.Month <=2 | times.Month == 12);

times(summer_x)
Chl_a_summer = Master_chl_a; 
Chl_a_summer = Chl_a_summer(:,:, summer_x); 

% For Fall: 

[fall_x, fall_y] = find(times.Month >= 3 & times.Month <= 5); 
times(fall_x) 
Chl_a_fall = Master_chl_a; 
Chl_a_fall = Chl_a_fall(:,:, fall_x); 

% So now for the raw and smoothn climatologies, you can average all of
% these and plot... (for raw) or average all of them and then smoothn on
% the 2D. 

Chl_a_winter_mean = mean(Chl_a_winter,3 ,'omitnan');
Chl_a_spring_mean = mean(Chl_a_spring, 3 , 'omitnan'); 
Chl_a_summer_mean = mean(Chl_a_summer, 3, 'omitnan'); 
Chl_a_fall_mean = mean(Chl_a_fall, 3, 'omitnan'); 

%% Adjust resolution 
% 
% step_Lat = Latitude_subset(1) : -1 : Latitude_subset(end) ; 
% step_Lon = Longitude_subset(1) : 1 : Longitude_subset(end); 
% step_Lat = flip(step_Lat); 

step_Lon = Aqua_Lon; 
step_Lat = Aqua_Lat;

fun = @(block_struct) nanmean(block_struct.data, [ 1 2]); 

Chl_a_winter_mean_one_degree_res = blockproc(Chl_a_winter_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_spring_mean_one_degree_res = blockproc(Chl_a_spring_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_summer_mean_one_degree_res = blockproc(Chl_a_summer_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_fall_mean_one_degree_res = blockproc(Chl_a_fall_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 
% 
%  Chl_a_winter_mean_one_degree_res = BlockMean(Chl_a_winter_mean,...
%         (length(Longitude_subset) ./ length(step_Lon)) ,...
%         length(Latitude_subset)  ./ length(step_Lat)); 
% 
%  Chl_a_spring_mean_one_degree_res = BlockMean(Chl_a_spring_mean,...
%      (length(Longitude_subset) ./ length(step_Lon)) ,...
%         length(Latitude_subset)  ./ length(step_Lat));
%     
%  Chl_a_summer_mean_one_degree_res = BlockMean(Chl_a_summer_mean,...
%      (length(Longitude_subset) ./ length(step_Lon)) ,...
%         length(Latitude_subset)  ./ length(step_Lat));
% 
%  Chl_a_fall_mean_one_degree_res = BlockMean(Chl_a_fall_mean,...
%      (length(Longitude_subset) ./ length(step_Lon)) ,...
%         length(Latitude_subset)  ./ length(step_Lat));

% climatology chl raw full degree Res 

Chl_a_summer_mean_rot = rot90(Chl_a_summer_mean); 
Chl_a_spring_mean_rot = rot90(Chl_a_spring_mean); 
Chl_a_fall_mean_rot   = rot90(Chl_a_fall_mean); 
Chl_a_winter_mean_rot = rot90(Chl_a_winter_mean); 

% climatology chl raw one degree res

Chl_a_summer_mean_rot = rot90(Chl_a_summer_mean_one_degree_res); 
Chl_a_spring_mean_rot = rot90(Chl_a_spring_mean_one_degree_res); 
Chl_a_fall_mean_rot   = rot90(Chl_a_fall_mean_one_degree_res); 
Chl_a_winter_mean_rot = rot90(Chl_a_winter_mean_one_degree_res); 

%% plot climatology chl raw figures full resolution:

%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

subplot(3,1,1)
plot_BSea_figure_subplot(Chl_a_spring_mean_rot, ...
    step_Lon,...
    step_Lat,...
    cmocean('algae'), ...
    [0 1],...
    'Spring'); 

subplot(3,1,2)
plot_BSea_figure_subplot(Chl_a_summer_mean_rot, ...
    step_Lon,...
    step_Lat,...
    cmocean('algae'), ...
    [0 1],...
    'Summer'); 

subplot(3,1,3)
plot_BSea_figure_subplot(Chl_a_fall_mean_rot, ...
    step_Lon,...
    step_Lat,...
    cmocean('algae'), ...
    [0 1],...
    'Fall'); 
    
% subplot(2,2,1)
% plot_BSea_figure_subplot(Chl_a_winter_mean_rot, ...
%     step_Lon,...
%     step_Lat,...
%     'algae', ...
%     [0 1.5],...
%     'Winter'); 
    
% saveas(gcf, 'Chlorophyll_one_deg_winter_smoothn.png'); 

hp4 = get(subplot(3,1, 3),'Position');
h = colorbar('Position', [hp4(1)+(hp4(3)-0.2)  0.33  0.02  0.4]);
h.FontWeight = 'bold';
h.FontSize = 15;
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
hy = ylabel(h, sprintf('mg m^{-3}'), 'FontSize', 18);
% hy.FontSize = 18;

%%

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'Updated_CHLA_Climatology_1_1','-dpng','-r300');       %  *// 300 dpi


%%

% %% trying to subplot seasonal climatologies so that they are all in the same figure. 
%    
% % now smoothn full res files, and  adjust the resolution 
% 
% Chl_a_winter_mean_smoothn = smoothn(Chl_a_winter_mean, 'robust', OPTIONS); 
% Chl_a_spring_mean_smoothn = smoothn(Chl_a_spring_mean, 'robust', OPTIONS); 
% Chl_a_summer_mean_smoothn = smoothn(Chl_a_summer_mean, 'robust', OPTIONS); 
% Chl_a_fall_mean_smoothn = smoothn(Chl_a_fall_mean, 'robust', OPTIONS); 
% 
%  Chl_a_winter_mean_one_degree_smoothn = BlockMean(Chl_a_winter_mean_smoothn,...
%         (length(Longitude_subset) ./ length(step_Lon)) ,...
%         length(Latitude_subset)  ./ length(step_Lat)); 
% 
%  Chl_a_spring_mean_one_degree_smoothn = BlockMean(Chl_a_spring_mean_smoothn,...
%      (length(Longitude_subset) ./ length(step_Lon)) ,...
%         length(Latitude_subset)  ./ length(step_Lat));
%     
%  Chl_a_summer_mean_one_degree_smoothn = BlockMean(Chl_a_summer_mean_smoothn,...
%      (length(Longitude_subset) ./ length(step_Lon)) ,...
%         length(Latitude_subset)  ./ length(step_Lat));
% 
%  Chl_a_fall_mean_one_degree_smoothn = BlockMean(Chl_a_fall_mean_smoothn,...
%      (length(Longitude_subset) ./ length(step_Lon)) ,...
%         length(Latitude_subset)  ./ length(step_Lat));
% 
% % chl vars: gaps from full res interpolated before constructing one degree
% % res figures. 
% 
% % one degree Res vars: 
% Chl_a_winter_mean_smoothn_rot = rot90(Chl_a_winter_mean_one_degree_smoothn); 
% Chl_a_spring_mean_smoothn_rot = rot90(Chl_a_spring_mean_one_degree_smoothn); 
% Chl_a_summer_mean_smoothn_rot = rot90(Chl_a_summer_mean_one_degree_smoothn); 
% Chl_a_fall_mean_smoothn_rot   = rot90(Chl_a_fall_mean_one_degree_smoothn); 


%%
% If you want the image to have land pixels masked, load 'WAP_landmask_updated.mat', below. I felt this was
% unecessary for the purpose of visualizaiton. 

% cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/LandPoints
% load('WAP_landmask_updated.mat') 

%%
% winter = WAP_landmask_updated .* Chl_a_winter_mean_smoothn_rot; 
% spring = WAP_landmask_updated .* Chl_a_spring_mean_smoothn_rot; 
% summer = WAP_landmask_updated.*Chl_a_summer_mean_smoothn_rot;
% fall = WAP_landmask_updated .* Chl_a_fall_mean_smoothn_rot;

% % one degree res... 
% winter = Chl_a_winter_mean_smoothn_rot; 
% spring = Chl_a_spring_mean_smoothn_rot; 
% summer = Chl_a_summer_mean_smoothn_rot;
% fall =  Chl_a_fall_mean_smoothn_rot;
% 
% %%
% 
% % here i'm plotting out the full res chl-climatologies that have also been
% % interpolated 
% 
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
% if ~make_it_tight,  clear subplot;  end
% 
% fig = figure,clf;
% 
% subplot(2,2,3)
% plot_BSea_figure_subplot(summer, ...
%     Aqua_Lon,...
%     Aqua_Lat,...
%     'algae', ...
%     [0 1.5],...
%     'Summer'); 
% 
% subplot(2,2,2)
% plot_BSea_figure_subplot(spring, ...
%     Aqua_Lon,...
%     Aqua_Lat,...
%     'algae', ...
%     [0 1.5],...
%     'Spring'); 
% 
% subplot(2,2,4)
% plot_BSea_figure_subplot(fall, ...
%     Aqua_Lon,...
%     Aqua_Lat,...
%     'algae', ...
%     [0 1.5],...
%     'Fall'); 
%     
% 
% subplot(2,2,1)
% plot_BSea_figure_subplot(winter, ...
%     Aqua_Lon,...
%     Aqua_Lat,...
%     'algae', ...
%     [0 1.5],...
%     'Winter'); 
%     
% hp4 = get(subplot(2,2,4),'Position');
% h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
% h.FontWeight = 'bold';
% h.FontSize = 15;
% hy = ylabel(h, sprintf('mg m^{-3}'), 'FontSize', 18);
% % hy.FontSize = 18;

%% CMOD, Winds, & Ice Now 

   % CALIPSO didn't start having data until June of 2006. This winter season will be shorter. 
   
    TR_winter_2006 = timerange('2006-06-01', '2006-09-01'); % Winter: June 1st to September 1st
    TR_winter_2007 = timerange('2007-06-01', '2007-09-01'); 
    TR_winter_2008 = timerange('2008-06-01', '2008-09-01');
    TR_winter_2009 = timerange('2009-06-01', '2009-09-01');
    TR_winter_2010 = timerange('2010-06-01', '2010-09-01'); 
    TR_winter_2011 = timerange('2011-06-01', '2011-09-01');
    TR_winter_2012 = timerange('2012-06-01', '2012-09-01');
    TR_winter_2013 = timerange('2013-06-01', '2013-09-01');
    TR_winter_2014 = timerange('2014-06-01', '2014-09-01');
    TR_winter_2015 = timerange('2015-06-01', '2015-09-01');
    TR_winter_2016 = timerange('2016-06-01', '2016-09-01');
    TR_winter_2017 = timerange('2017-06-01', '2017-09-01');
    TR_winter_2018 = timerange('2018-06-01', '2018-09-01'); 
    
    TR_spring_2006 = timerange('2006-09-01', '2006-12-01'); % Spring: Sept 1st to Dec 1st
    TR_spring_2007 = timerange('2007-09-01', '2007-12-01');
    TR_spring_2008 = timerange('2008-09-01', '2008-12-01'); 
    TR_spring_2009 = timerange('2009-09-01', '2009-12-01'); 
    TR_spring_2010 = timerange('2010-09-01', '2010-12-01'); 
    TR_spring_2011 = timerange('2011-09-01', '2011-12-01');
    TR_spring_2012 = timerange('2012-09-01', '2012-12-01'); 
    TR_spring_2013 = timerange('2013-09-01', '2013-12-01'); 
    TR_spring_2014 = timerange('2014-09-01', '2014-12-01');
    TR_spring_2015 = timerange('2015-09-01', '2015-12-01'); 
    TR_spring_2016 = timerange('2016-09-01', '2016-12-01'); 
    TR_spring_2017 = timerange('2017-09-01', '2017-12-01');
    TR_spring_2018 = timerange('2018-09-01', '2018-12-01'); 
    
    TR_summer_2006 = timerange('2006-12-01', '2007-03-01'); % Summer: Dec 1st to March 1st
    TR_summer_2007 = timerange('2007-12-01', '2008-03-01');
    TR_summer_2008 = timerange('2008-12-01', '2009-03-01'); 
    TR_summer_2009 = timerange('2009-12-01', '2010-03-01'); 
    TR_summer_2010 = timerange('2010-12-01', '2011-03-01'); 
    TR_summer_2011 = timerange('2011-12-01', '2012-03-01');
    TR_summer_2012 = timerange('2012-12-01', '2013-03-01'); 
    TR_summer_2013 = timerange('2013-12-01', '2014-03-01'); 
    TR_summer_2014 = timerange('2014-12-01', '2015-03-01');
    TR_summer_2015 = timerange('2015-12-01', '2016-03-01'); 
    TR_summer_2016 = timerange('2016-12-01', '2017-03-01'); 
    TR_summer_2017 = timerange('2017-12-01', '2018-03-01');
    TR_summer_2018 = timerange('2018-12-01', '2019-01-01'); 
    
    %     TR_fall_2006   = timerange('2006-03-01', '2006-06-01'); % there
    %     is no CALIPSO DATA in fall of 2006
    TR_fall_2007   = timerange('2007-03-01', '2007-06-01'); %% Fall: March 1st to June 1st
    TR_fall_2008   = timerange('2008-03-01', '2008-06-01');
    TR_fall_2009   = timerange('2009-03-01', '2009-06-01');
    TR_fall_2010   = timerange('2010-03-01', '2010-06-01');
    TR_fall_2011   = timerange('2011-03-01', '2011-06-01');
    TR_fall_2012   = timerange('2012-03-01', '2012-06-01');
    TR_fall_2013   = timerange('2013-03-01', '2013-06-01');
    TR_fall_2014   = timerange('2014-03-01', '2014-06-01');
    TR_fall_2015   = timerange('2015-03-01', '2015-06-01');
    TR_fall_2016   = timerange('2016-03-01', '2016-06-01');
    TR_fall_2017   = timerange('2017-03-01', '2017-06-01');
    TR_fall_2018   = timerange('2018-03-01', '2018-06-01');

%%
% Initializing of my variables for the loop below. 

Total_CMOD_winter = [];
Total_Lat_CMOD_winter = [];
Total_Lon_CMOD_winter = [];

Total_CMOD_spring = [];
Total_Lat_CMOD_spring = [];
Total_Lon_CMOD_spring = [];

Total_CMOD_summer = [];
Total_Lat_CMOD_summer = [];
Total_Lon_CMOD_summer = [];

Total_CMOD_fall = [];
Total_Lat_CMOD_fall = [];
Total_Lon_CMOD_fall = [];


Total_Ice_winter = [];
Total_Lat_Ice_winter = [];
Total_Lon_Ice_winter = [];

Total_Ice_spring = [];
Total_Lat_Ice_spring = [];
Total_Lon_Ice_spring = [];

Total_Ice_summer = [];
Total_Lat_Ice_summer = [];
Total_Lon_Ice_summer = [];

Total_Ice_fall = [];
Total_Lat_Ice_fall = [];
Total_Lon_Ice_fall = [];

Total_Wind_winter = [];
Total_Lat_Wind_winter = [];
Total_Lon_Wind_winter = [];

Total_Wind_spring = [];
Total_Lat_Wind_spring = [];
Total_Lon_Wind_spring = [];

Total_Wind_summer = [];
Total_Lat_Wind_summer = [];
Total_Lon_Wind_summer = [];

Total_Wind_fall = [];
Total_Lat_Wind_fall = [];
Total_Lon_Wind_fall = [];

%%
for i = 2006:2018
    
    disp(i)
    
    SEASON = {'winter', 'spring', 'summer', 'fall'};
    
    for j = 1:length(SEASON)
        
        disp(j)
        
        % [eval(sprintf('CMOD_Bin_%d_winter', i)), eval(sprintf('night_%d_winter_occ', i)),...
        %     eval(sprintf('night_%d_winter_std')), eval(sprintf('night_%d_winter_err', i))] = ...
        
        if exist(sprintf('TR_%s_%d', SEASON{j}, i),'var') % only keep looping if the variable actually exists
            % TR_fall_2006 does not exist, which is why I need this
            % statement in here
            
            eval(sprintf('Lat_CMOD = Total_timetable_CMOD_Surface(TR_%s_%d,:).Total_Latitude_Surface;', SEASON{j},i))
            eval(sprintf('Lon_CMOD = Total_timetable_CMOD_Surface(TR_%s_%d, :).Total_Longitude_Surface;',SEASON{j}, i))
            eval(sprintf('OD       = Total_timetable_CMOD_Surface(TR_%s_%d, :).CMOD_Surface;', SEASON{j}, i))
            
            eval(sprintf('Lat_Ice = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Latitude_Ice;', SEASON{j},i))
            eval(sprintf('Lon_Ice = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Longitude_Ice;', SEASON{j},i))
            eval(sprintf('Ice     = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Surface_532_Integrated_Depolarization_Ratio;', SEASON{j},i)) 
            
            
            % I dont need to filter this out anymore since the
            % 'Saving_out_CMOD_Ice_Wind.m' file took care of this. 
%             
%             bad_Ice_values = Ice <= -0.2 | Ice > 1.2; 
%             Ice(bad_Ice_values) = NaN;
%             
%             nan_ice        = isnan(Ice(:,1)); 
%             Ice      = Ice(~nan_ice) ;
%             Lat_Ice  = Lat_Ice(~nan_ice); 
%             Lon_Ice  = Lon_Ice(~nan_ice); 
%             
            
            
            
            eval(sprintf('Lat_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Latitude_Wind;', SEASON{j}, i))
            eval(sprintf('Lon_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Longitude_Wind;', SEASON{j}, i))
            eval(sprintf('Wind     = Total_timetable_amsrmf(TR_%s_%d,:).Total_windamsrMF;', SEASON{j}, i))
         
            
            if j == 1
                
                Total_CMOD_winter = vertcat(Total_CMOD_winter, OD);                 
                Total_Lat_CMOD_winter = vertcat(Total_Lat_CMOD_winter, Lat_CMOD);
                Total_Lon_CMOD_winter = vertcat(Total_Lon_CMOD_winter, Lon_CMOD);
                
                Total_Ice_winter = vertcat(Total_Ice_winter, Ice);
                Total_Lat_Ice_winter = vertcat(Total_Lat_Ice_winter, Lat_Ice); 
                Total_Lon_Ice_winter = vertcat(Total_Lon_Ice_winter, Lon_Ice); 
                                
                Total_Wind_winter = vertcat(Total_Wind_winter, Wind); 
                Total_Lat_Wind_winter = vertcat(Total_Lat_Wind_winter, Lat_Wind); 
                Total_Lon_Wind_winter = vertcat(Total_Lon_Wind_winter, Lon_Wind); 
                
            elseif j == 2
                
                Total_CMOD_spring = vertcat(Total_CMOD_spring, OD);
                Total_Lat_CMOD_spring = vertcat(Total_Lat_CMOD_spring, Lat_CMOD);
                Total_Lon_CMOD_spring = vertcat(Total_Lon_CMOD_spring, Lon_CMOD);
                
                Total_Ice_spring = vertcat(Total_Ice_spring, Ice);
                Total_Lat_Ice_spring = vertcat(Total_Lat_Ice_spring, Lat_Ice); 
                Total_Lon_Ice_spring = vertcat(Total_Lon_Ice_spring, Lon_Ice); 
                                
                Total_Wind_spring = vertcat(Total_Wind_spring, Wind); 
                Total_Lat_Wind_spring = vertcat(Total_Lat_Wind_spring, Lat_Wind); 
                Total_Lon_Wind_spring = vertcat(Total_Lon_Wind_spring, Lon_Wind); 
                
            elseif j == 3
                
                Total_CMOD_summer = vertcat(Total_CMOD_summer, OD);               
                Total_Lat_CMOD_summer = vertcat(Total_Lat_CMOD_summer, Lat_CMOD);
                Total_Lon_CMOD_summer = vertcat(Total_Lon_CMOD_summer, Lon_CMOD);
                
                Total_Ice_summer = vertcat(Total_Ice_summer, Ice);
                Total_Lat_Ice_summer = vertcat(Total_Lat_Ice_summer, Lat_Ice); 
                Total_Lon_Ice_summer = vertcat(Total_Lon_Ice_summer, Lon_Ice); 
                                
                Total_Wind_summer = vertcat(Total_Wind_summer, Wind); 
                Total_Lat_Wind_summer = vertcat(Total_Lat_Wind_summer, Lat_Wind); 
                Total_Lon_Wind_summer = vertcat(Total_Lon_Wind_summer, Lon_Wind); 
                
            elseif j == 4
                
                Total_CMOD_fall = vertcat(Total_CMOD_fall, OD);                
                Total_Lat_CMOD_fall = vertcat(Total_Lat_CMOD_fall, Lat_CMOD);
                Total_Lon_CMOD_fall = vertcat(Total_Lon_CMOD_fall, Lon_CMOD);
                
                Total_Ice_fall = vertcat(Total_Ice_fall, Ice);
                Total_Lat_Ice_fall = vertcat(Total_Lat_Ice_fall, Lat_Ice); 
                Total_Lon_Ice_fall = vertcat(Total_Lon_Ice_fall, Lon_Ice); 
                                
                Total_Wind_fall = vertcat(Total_Wind_fall, Wind); 
                Total_Lat_Wind_fall = vertcat(Total_Lat_Wind_fall, Lat_Wind); 
                Total_Lon_Wind_fall = vertcat(Total_Lon_Wind_fall, Lon_Wind); 
                
            end
     
            
        end
        
        clear CMOD CMOD_OCC CMOD_STD CMOD_ERR 
        
    end
end

%%

[CMOD_one_degree_winter, CMOD_OCC_winter, CMOD_STD_winter, CMOD_ERR_winter]  = hist_wt_occ_tot(Total_Lat_CMOD_winter, Total_Lon_CMOD_winter, Total_CMOD_winter, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_spring, CMOD_OCC_spring, CMOD_STD_spring, CMOD_ERR_spring]  = hist_wt_occ_tot(Total_Lat_CMOD_spring, Total_Lon_CMOD_spring, Total_CMOD_spring, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_summer, CMOD_OCC_summer, CMOD_STD_summer, CMOD_ERR_summer]  = hist_wt_occ_tot(Total_Lat_CMOD_summer, Total_Lon_CMOD_summer, Total_CMOD_summer, Aqua_Lat', Aqua_Lon');
[CMOD_one_degree_fall, CMOD_OCC_fall, CMOD_STD_fall, CMOD_ERR_fall]          = hist_wt_occ_tot(Total_Lat_CMOD_fall, Total_Lon_CMOD_fall, Total_CMOD_fall, Aqua_Lat', Aqua_Lon');



[Ice_one_degree_winter, Ice_OCC_winter, Ice_STD_winter, Ice_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Ice_winter, Total_Lon_Ice_winter, Total_Ice_winter, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_spring, Ice_OCC_spring, Ice_STD_spring, Ice_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Ice_spring, Total_Lon_Ice_spring, Total_Ice_spring, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_summer, Ice_OCC_summer, Ice_STD_summer, Ice_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Ice_summer, Total_Lon_Ice_summer, Total_Ice_summer, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_fall, Ice_OCC_fall, Ice_STD_fall, Ice_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Ice_fall, Total_Lon_Ice_fall, Total_Ice_fall, Aqua_Lat', Aqua_Lon');




[Wind_one_degree_winter, Wind_OCC_winter, Wind_STD_winter, Wind_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Wind_winter, Total_Lon_Wind_winter, Total_Wind_winter, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_spring, Wind_OCC_spring, Wind_STD_spring, Wind_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Wind_spring, Total_Lon_Wind_spring, Total_Wind_spring, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_summer, Wind_OCC_summer, Wind_STD_summer, Wind_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Wind_summer, Total_Lon_Wind_summer, Total_Wind_summer, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_fall, Wind_OCC_fall, Wind_STD_fall, Wind_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Wind_fall, Total_Lon_Wind_fall, Total_Wind_fall, Aqua_Lat', Aqua_Lon');


%%
winter = Wind_one_degree_winter; 
spring = Wind_one_degree_spring; 
summer = Wind_one_degree_summer;
fall   =  Wind_one_degree_fall;

%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure,clf;

subplot(2,2,3)
plot_BSea_figure_subplot(summer, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('amp'), ...
    [0 13],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(spring, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('amp'), ...
    [0 13],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(fall, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('amp'), ...
    [0 13],...
    'Fall'); 
    
subplot(2,2,1)
plot_BSea_figure_subplot(winter, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('amp'), ...
    [0 13],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, sprintf('m s^{-1}'), 'FontSize', 18);
% hy.FontSize = 18;
%%

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'Updated_WindSpeed_Climatology_1_1','-dpng','-r300');       %  *// 300 dpi


%%

winter = Ice_one_degree_winter; 
spring = Ice_one_degree_spring; 
summer = Ice_one_degree_summer;
fall   = Ice_one_degree_fall;


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

subplot(2,2,3)
plot_BSea_figure_subplot(summer, ...
    step_Lon,...
    step_Lat,...
    cmocean('Ice'), ...
    [0 0.8],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(spring, ...
    step_Lon,...
    step_Lat,...
    cmocean('Ice'), ...
    [0 0.8],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(fall, ...
    step_Lon,...
    step_Lat,...
    cmocean('Ice'), ...
    [0 0.8],...
    'Fall'); 
    

subplot(2,2,1)
plot_BSea_figure_subplot(winter, ...
    step_Lon,...
    step_Lat,...
    cmocean('Ice'), ...
    [0 0.8],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, 'Depolarization Ratio \delta', 'FontSize', 16);
% hy.FontSize = 18;
%%

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'Updated_Ice_Climatology_1_1','-dpng','-r300');       %  *// 300 dpi


%%
% Smoothn Climatologies 
% 
% CMOD_winter_smoothn = smoothn(CMOD_one_degree_winter);
% CMOD_spring_smoothn = smoothn(CMOD_one_degree_spring); 
% CMOD_summer_smoothn = smoothn(CMOD_one_degree_summer); 
% CMOD_fall_smoothn = smoothn(CMOD_one_degree_fall); 
% 
% 
% %%
% winter = CMOD_winter_smoothn; 
% spring = CMOD_spring_smoothn; 
% summer = CMOD_summer_smoothn;
% fall   = CMOD_fall_smoothn;

%%

winter = CMOD_one_degree_winter; 
spring = CMOD_one_degree_spring;
summer = CMOD_one_degree_summer; 
fall   = CMOD_one_degree_fall;




make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

subplot(2,2,3)
plot_BSea_figure_subplot(summer, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('tempo'), ...
    [0 .08],...
    'Summer'); 

subplot(2,2,2)
plot_BSea_figure_subplot(spring, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('tempo'), ...
    [0 .08],...
    'Spring'); 

subplot(2,2,4)
plot_BSea_figure_subplot(fall, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('tempo'), ...
    [0 .08],...
    'Fall'); 
    

subplot(2,2,1)
plot_BSea_figure_subplot(winter, ...
    Aqua_Lon,...
    Aqua_Lat,...
    cmocean('tempo'), ...
    [0 .08],...
    'Winter'); 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, 'MAOD', 'FontSize', 18);
% hy.FontSize = 18;

%%

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'Updated_MAOD_Climatology_1_1','-dpng','-r300');       %  *// 300 dpi

%%


% %%
% Wind_winter_smoothn = smoothn(Wind_one_degree_winter);
% Wind_spring_smoothn = smoothn(Wind_one_degree_spring); 
% Wind_summer_smoothn = smoothn(Wind_one_degree_summer); 
% Wind_fall_smoothn = smoothn(Wind_one_degree_fall); 
% 
% %%
% winter = Wind_winter_smoothn; 
% spring = Wind_spring_smoothn; 
% summer = Wind_summer_smoothn;
% fall   = Wind_fall_smoothn;
% % 
% 
% %%
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
% if ~make_it_tight,  clear subplot;  end
% 
% fig = figure,clf;
% 
% subplot(2,2,3)
% plot_BSea_figure_subplot(summer, ...
%     step_Lon,...
%     step_Lat,...
%     'amp', ...
%     [0 20],...
%     'Summer'); 
% 
% subplot(2,2,2)
% plot_BSea_figure_subplot(spring, ...
%     step_Lon,...
%     step_Lat,...
%     'amp', ...
%     [0 20],...
%     'Spring'); 
% 
% subplot(2,2,4)
% plot_BSea_figure_subplot(fall, ...
%     step_Lon,...
%     step_Lat,...
%     'amp', ...
%     [0 20],...
%     'Fall'); 
%     
% 
% subplot(2,2,1)
% plot_BSea_figure_subplot(winter, ...
%     step_Lon,...
%     step_Lat,...
%     'amp', ...
%     [0 20],...
%     'Winter'); 
%     
% hp4 = get(subplot(2,2,4),'Position');
% h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
% h.FontWeight = 'bold';
% h.FontSize = 15;
% hy = ylabel(h, sprintf('m s^{-1}'), 'FontSize', 18);
% % hy.FontSize = 18;

%% Time Series Climatologies...


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication

load('Depol_Ratio_Monthly_avg_Vars.mat')
load('amsrmf_Monthly_avg_Vars.mat')
load('CMOD_Monthly_avg_Vars_Surface.mat') 


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Revisions_Figures

t1 = datetime(2006,06,01);
t2 = datetime(2018,12,31);
times = t1:calmonths(1):t2; 



t1 = datetime(2007,01,01);
t2 = datetime(2007,12,31);
times_months = t1:calmonths(1):t2; 
times_months_num = month(times_months);

times = times'; 
Ice_Vector = double(Depol_Ratio_Monthly_avg(:,1));
Ice_Vector(117,:) = [];
CMOD_Vector = CMOD_Monthly_avg_Surface(:,1);
CMOD_Vector(117,:) = [];
Wind_Vector = amsrmf_Monthly_avg;
times_without_feb = times;
times_without_feb(117,:) = []; %
times_for_wind = times;
times_for_wind(117,:) =[];
times_for_wind(64:73,:) = [];
Wind_Vector(117,:) = [];
Wind_Vector(64:73,:)=[];



[Ac_CMOD, tc_CMOD] = climatology(CMOD_Vector,times_without_feb, 'monthly'); 
[Ac_Ice, tc_Ice] = climatology(Ice_Vector, times_without_feb, 'monthly'); 
[Ac_Wind, tc_Wind] = climatology(Wind_Vector, times_for_wind, 'monthly'); 

Master_chl_a(isnan(Master_chl_a)) = 0; 
[Ac_Chl, tc_Chl] = climatology(Master_chl_a, times, 'monthly'); 

for i = 1:12
Ac_Chl_monthly(i) = mean2(Ac_Chl(:,:,i));
end

Ac_Chl_monthly(Ac_Chl_monthly <= 0) = NaN;


%%
% this is for standard deviation calculation, check cdt monthly to see what
% I did 

jan_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[1], @mad);
feb_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[2], @mad);
mar_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[3], @mad);
apr_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[4],  @mad);
may_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[5],  @mad);
jun_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[6],  @mad);
jul_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[7],  @mad);
aug_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[8],  @mad);
sep_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[9],  @mad);
oct_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[10],  @mad);
nov_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[11],  @mad);
dec_CMOD_mad = monthly(CMOD_Vector,times_without_feb,[12],  @mad);

[Ac_Ice, tc_Ice] = climatology(Ice_Vector, times_without_feb, 'monthly'); 


jan_Ice_mad = monthly(Ice_Vector,times_without_feb,[1],  @mad);
feb_Ice_mad = monthly(Ice_Vector,times_without_feb,[2],  @mad);
mar_Ice_mad = monthly(Ice_Vector,times_without_feb,[3],  @mad);
apr_Ice_mad = monthly(Ice_Vector,times_without_feb,[4],  @mad);
may_Ice_mad = monthly(Ice_Vector,times_without_feb,[5],  @mad);
jun_Ice_mad = monthly(Ice_Vector,times_without_feb,[6],  @mad);
jul_Ice_mad = monthly(Ice_Vector,times_without_feb,[7],  @mad);
aug_Ice_mad = monthly(Ice_Vector,times_without_feb,[8],  @mad);
sep_Ice_mad = monthly(Ice_Vector,times_without_feb,[9],  @mad);
oct_Ice_mad = monthly(Ice_Vector,times_without_feb,[10],  @mad);
nov_Ice_mad = monthly(Ice_Vector,times_without_feb,[11],  @mad);
dec_Ice_mad = monthly(Ice_Vector,times_without_feb,[12],  @mad);


[Ac_Wind, tc_Wind] = climatology(Wind_Vector, times_for_wind, 'monthly'); 

jan_Wind_mad = monthly(Wind_Vector,times_for_wind,[1],  @mad);
feb_Wind_mad = monthly(Wind_Vector,times_for_wind,[2],  @mad);
mar_Wind_mad = monthly(Wind_Vector,times_for_wind,[3],  @mad);
apr_Wind_mad = monthly(Wind_Vector,times_for_wind,[4],  @mad);
may_Wind_mad = monthly(Wind_Vector,times_for_wind,[5],  @mad);
jun_Wind_mad = monthly(Wind_Vector,times_for_wind,[6],  @mad);
jul_Wind_mad = monthly(Wind_Vector,times_for_wind,[7],  @mad);
aug_Wind_mad = monthly(Wind_Vector,times_for_wind,[8],  @mad);
sep_Wind_mad = monthly(Wind_Vector,times_for_wind,[9],  @mad);
oct_Wind_mad = monthly(Wind_Vector,times_for_wind,[10],  @mad);
nov_Wind_mad = monthly(Wind_Vector,times_for_wind,[11],  @mad);
dec_Wind_mad = monthly(Wind_Vector,times_for_wind,[12],  @mad);

% chl and standard deviation calculation
Master_chl_a(isnan(Master_chl_a)) = 0; 
[Ac_Chl, tc_Chl] = climatology(Master_chl_a, times, 'monthly'); 

for i = 1:12
Ac_Chl_monthly(i) = mean2(Ac_Chl(:,:,i));
end

Ac_Chl_monthly(Ac_Chl_monthly <= 0) = NaN;

Ac_Chl_monthly = Ac_Chl_monthly';

for i = 1:12 
    disp(i)
    mad_chl_monthly_2d(:,:,i) = monthly(Master_chl_a, times, [i], 'dim', 3, @mad); 
end

    mad_chl_monthly = nanmean(mad_chl_monthly_2d, [1 2]);
    mad_chl_monthly = squeeze(mad_chl_monthly); 
    
% jan_chl_std = monthly(Master_chl_a,times,[1], 'omitnan', @std);
% jan_chl_std = 
% feb_chl_std = monthly(Master_chl_a,times,[2], 'omitnan', @std);
% mar_chl_std = monthly(Master_chl_a,times,[3], 'omitnan', @std);
% apr_chl_std = monthly(Master_chl_a,times,[4], 'omitnan', @std);
% may_chl_std = monthly(Master_chl_a,times,[5], 'omitnan', @std);
% jun_chl_std = monthly(Master_chl_a,times,[6], 'omitnan', @std);
% jul_chl_std = monthly(Master_chl_a,times,[7], 'omitnan', @std);
% aug_chl_std = monthly(Master_chl_a,times,[8], 'omitnan', @std);
% sep_chl_std = monthly(Master_chl_a,times,[9], 'omitnan', @std);
% oct_chl_std = monthly(Master_chl_a,times,[10], 'omitnan', @std);
% nov_chl_std = monthly(Master_chl_a,times,[11], 'omitnan', @std);
% dec_chl_std = monthly(Master_chl_a,times,[12], 'omitnan', @std);





%%
sienna = rgb('sienna'); 
gray = rgb('gray'); 
black = rgb('black'); 
grass_green = rgb('grass green');
ice_blue = rgb('sky blue'); 
wind_blue = rgb('scarlet');

% Create textbox
% pos = [x-start y-start x-width y-height]


%%

make_it_tight = true;
 subplot = @(m,n,p) subtightplot (m, n, p, [0.035 0.05], [0.1 0.05], [0.1 0.05]);
if ~make_it_tight,  clear subplot;  end


fig = figure(1), clf;

ax(1) = subplot(4,1,1);
grid on
aa_splot(1:12, Ac_CMOD, 'x-',...
    'linewidth', 2, ...
    'Color', black,...
    'MarkerSize', 9,...
    'MarkerFaceColor', black,...
    'MarkerEdgeColor', black);

ylabel('MAOD')
%         xlabel('Months')
xlim([1,12])
ylim([0.01 0.09])
set(gca, 'ytick',0: 0.02: 0.08); 
AX=findall(0,'type','axes');
set(AX, 'FontSize', 18)
%                   xtickangle(20)


ax(2) = subplot(4,1,2);
grid on
aa_splot(1:12, Ac_Chl_monthly, 'v-',...
    'linewidth', 2, ...
    'Color', grass_green,...
    'MarkerSize', 7,...
    'MarkerFaceColor', grass_green,...
    'MarkerEdgeColor', grass_green);
% set(gca, 'xtick', 1:12,...
%     'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Chl-{\ita} Climatology')
ylabel(sprintf('mg m^{-3}')),...
    %         xlabel('Months')
xlim([1,12])
ylim([0.0001 0.3])
set(gca, 'ytick',0: 0.1: 0.3); 


AX=findall(0,'type','axes');
set(AX, 'FontSize', 18)
% xtickangle(45)




ax(3) = subplot(4,1,3); 
grid on
aa_splot(1:12, Ac_Ice, 'o-',...
    'linewidth', 2, ...
    'Color', ice_blue,...
    'MarkerSize', 7,...
    'MarkerFaceColor', ice_blue,...
    'MarkerEdgeColor', ice_blue);
% set(gca, 'xtick', 1:12,...
%     'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Ice Climatology')
ylabel('\delta',...
   'FontName','Helvetica Neue');%         xlabel('Months')
xlim([1,12])
ylim([0.4 0.7])
set(gca, 'ytick',0.4: 0.1: 0.7); 


AX=findall(0,'type','axes');
set(AX, 'FontSize', 18)
% xtickangle(45)




ax(4) = subplot(4,1,4);
grid on
aa_splot(1:12, Ac_Wind, 'd-',...
    'linewidth', 2, ...
    'Color', wind_blue,...
    'MarkerSize', 7,...
    'MarkerFaceColor', wind_blue,...
    'MarkerEdgeColor', wind_blue);
set(gca, 'xtick', 1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% title('Wind Speed Climatology')
ylabel(sprintf('m s^{-1}'))
xlim([1,12])
ylim([6 13])

AX=findall(0,'type','axes');
set(AX, 'FontSize', 20)
 xtickangle(45)


[ax(1:3).XTickLabel] = deal([]);
ax(1).YTick(1) = []; 
ax(2).YTick(1) = []; 
ax(3).YTick(1) = []; 
ax(4).YTick(1) = [];
ax(4).FontSize = 23;
% ax(1) to see properties 

% from below: an(1).Position(2) = 0.87
% ax(1).Position(2) = 0.7673

% (0.7673 - should give the correct position for first subplot 

an(1) = annotation(gcf,'textbox',... % I drew this on in the figure 
  [0.11729695024077 0.915 0.16 0.0349],...
    'String','MAOD',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(2) = annotation(gcf,'textbox',...
    [an(1).Position(1),...
    (ax(2).Position(2) + (an(1).Position(2) - ax(1).Position(2))),... % so this is the second axes, y position  + (y position of text box - y position of first axes)
    0.16 0.0349],...
    'String',{'Chl-{\ita}'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(3) = annotation(gcf, 'textbox', ...
    [an(1).Position(1),...
    (ax(3).Position(2) + (an(2).Position(2) - ax(2).Position(2))),...
    0.16 0.0349],...
    'String',{'Ice'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

an(4) = annotation(gcf, 'textbox', ...
    [an(1).Position(1),...
    (ax(4).Position(2) + (an(3).Position(2) - ax(3).Position(2))),...
    0.16 0.0349],...
    'String',{'Wind'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');


%%


% set(gcf,'PaperPositionMode','auto')
%  set(gcf,'PaperPosition','fillpage') 
orient(fig,'landscape')

print(fig,'Updated_Climatology_TimeSeries','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju


%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end
%% Upper and Lower Subplots with Titles
income = [3.2,4.1,5.0,5.6];
outgo = [2.5,4.0,3.35,4.9];
subplot(2,1,1); plot(income)
title('Income')
subplot(2,1,2); plot(outgo)
title('Outgo')

%% Subplots in Quadrants
figure
subplot(2,2,1)
text(.5,.5,{'subplot(2,2,1)';'or subplot 221'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,2)
text(.5,.5,{'subplot(2,2,2)';'or subplot 222'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,3)
text(.5,.5,{'subplot(2,2,3)';'or subplot 223'},...
    'FontSize',14,'HorizontalAlignment','center')
subplot(2,2,4)
text(.5,.5,{'subplot(2,2,4)';'or subplot 224'},...
    'FontSize',14,'HorizontalAlignment','center')
