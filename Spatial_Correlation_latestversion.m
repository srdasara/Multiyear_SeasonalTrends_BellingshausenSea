%% Now, taking my seasonal grids of CMOD, Ice, Winds, and CHL-A, I will run a correlation through all of them and 
% Producing a variant of the surf plot 


 cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/

%%

load Total_timetable_amsrmf.mat
load Total_timetable_CMOD_Surface.mat
load Total_timetable_Depol_Ratio.mat


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/New_mat_files_CHL_A_DATA_MONTHLY/
load('Master_chlor_a_monthly_full_res.mat') 
load('Aqua_Longitude_Subset_BellingshausenSea.mat')
load('Aqua_Latitude_Subset_BellingshausenSea.mat') 

% for i = 1:151 
% Total_chl_monthly_2(i) = nanmean(Master_chl_a(:,:,i), [1 2]); 
% end
% 
% Total_chl_a_monthly = Total_chl_monthly_2; 

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Revisions_Figures

%%
step_Lat = (linspace(Latitude_subset(1), Latitude_subset(end), 15)) ; 
step_Lon = (linspace(Longitude_subset(1), Longitude_subset(end), 46)); 
step_Lat = flip(step_Lat); % this is for plotting purposes because pcolor flips image on y axis. 



% At the moment this is just to check how many steps they are and if they
% correspond...
nstep_Lat = Latitude_subset(1) : -1 : Latitude_subset(end) ; 
nstep_Lon = Longitude_subset(1) : 1 : Longitude_subset(end); 
nstep_Lat = flip(nstep_Lat); 


%% Here I am averaging data from each individual season for the temporal range of record.
clear Chl_a_seasons
count = 0;
for i = 1:3:148
    count = count+1; 
    disp(count)
    Chl_a_seasons(:,:,count) = mean(Master_chl_a(:,:, i : i+2), 3, 'omitnan'); 
    
end

Chl_a_seasons = cat(3,Chl_a_seasons, Master_chl_a(:,:,151)); % concatenating the last month since it's the remainder 

%% This is a better way to adjust resolution of chl-a, that way, you can use nanmean for the missing values.

fun = @(block_struct) nanmean(block_struct.data, [ 1 2]); 


CHL_SEASONS = blockproc(Chl_a_seasons,...
    [length(Chl_a_seasons(:,1, :)) ./ length(step_Lon) length(Chl_a_seasons(1,:, 1)) ./ length(step_Lat)],...
     fun); 






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
    % 13
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
    % 13
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
    % 13
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
% 12

% 51 seasons in total across June 2006 to Dec 2018 


%%

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



clear CMOD_SEASONS SeaIce_SEASONS Winds_SEASONS
for i = 2006:2018
    
    disp(i)
    
    SEASON = {'fall','winter', 'spring', 'summer'};
    
    for j = 1:length(SEASON)
        
        disp(j)
        
        % [eval(sprintf('CMOD_Bin_%d_winter', i)), eval(sprintf('night_%d_winter_occ', i)),...
        %     eval(sprintf('night_%d_winter_std')), eval(sprintf('night_%d_winter_err', i))] = ...
        
        if exist(sprintf('TR_%s_%d', SEASON{j}, i),'var') % only keep looping if the variable actually exists
            % TR_fall_2006 does not exist, which is why I need this
            % statement in here. It will skip fall 2006, which doesn't
            % exist, and start with winter 2006, proceeding onwards until
            % end of the timeframe.. 
            
            eval(sprintf('Lat_CMOD = Total_timetable_CMOD_Surface(TR_%s_%d,:).Total_Latitude_Surface;', SEASON{j},i))
            eval(sprintf('Lon_CMOD = Total_timetable_CMOD_Surface(TR_%s_%d, :).Total_Longitude_Surface;',SEASON{j}, i))
            eval(sprintf('OD       = Total_timetable_CMOD_Surface(TR_%s_%d, :).CMOD_Surface;', SEASON{j}, i))
            
            eval(sprintf('Lat_Ice = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Latitude_Ice;', SEASON{j},i))
            eval(sprintf('Lon_Ice = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Longitude_Ice;', SEASON{j},i))
            eval(sprintf('Ice     = Total_timetable_Depol_Ratio(TR_%s_%d,:).Total_Surface_532_Integrated_Depolarization_Ratio;', SEASON{j},i)) 
            
            eval(sprintf('Lat_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Latitude_Wind;', SEASON{j}, i))
            eval(sprintf('Lon_Wind = Total_timetable_amsrmf(TR_%s_%d,:).Total_Longitude_Wind;', SEASON{j}, i))
            eval(sprintf('Wind     = Total_timetable_amsrmf(TR_%s_%d,:).Total_windamsrMF;', SEASON{j}, i))
                                                
            [CMOD, CMOD_OCC, CMOD_STD, CMOD_ERR]         = hist_wt_occ_tot(Lat_CMOD, Lon_CMOD, OD, step_Lat', step_Lon');
            [SeaIce, SeaIce_OCC, SeaIce_STD, SeaIce_ERR] = hist_wt_occ_tot(Lat_Ice, Lon_Ice, Ice, step_Lat', step_Lon');
            [Winds, Winds_OCC, Winds_STD, Winds_ERR]     = hist_wt_occ_tot(Lat_Wind, Lon_Wind, Wind, step_Lat', step_Lon');
            
            
            
            
            
            % CMOD 3D Seasons Now
            if ~exist('CMOD_SEASONS', 'var')
                
                CMOD_SEASONS     = CMOD;
                SeaIce_SEASONS   = SeaIce;
                Winds_SEASONS    = Winds;
                
                CMOD_OCC_SEASONS = CMOD_OCC;
                CMOD_STD_SEASONS = CMOD_STD;
                CMOD_ERR_SEASONS = CMOD_ERR;
                
                SeaIce_OCC_SEASONS = SeaIce_OCC;
                SeaIce_STD_SEASONS = SeaIce_STD;
                SeaIce_ERR_SEASONS = SeaIce_ERR;
                
                Winds_OCC_SEASONS = Winds_OCC;
                Winds_STD_SEASONS = Winds_STD;
                Winds_ERR_SEASONS = Winds_ERR;
                
            else
                
                CMOD_SEASONS(:,:, end + 1) = cat(3, CMOD);
                SeaIce_SEASONS(:,:, end + 1) = cat(3, SeaIce);
                Winds_SEASONS(:,:, end + 1) = cat(3, Winds);
                
                CMOD_OCC_SEASONS(:,:, end + 1) = cat(3, CMOD_OCC);
                CMOD_STD_SEASONS(:,:, end + 1) = cat(3, CMOD_STD);
                CMOD_ERR_SEASONS(:,:, end + 1) = cat(3, CMOD_ERR);
                
                SeaIce_OCC_SEASONS(:,:, end + 1) = cat(3, SeaIce_OCC);
                SeaIce_STD_SEASONS(:,:, end + 1) = cat(3, SeaIce_STD);
                SeaIce_ERR_SEASONS(:,:, end + 1) = cat(3, SeaIce_ERR);
                
                Winds_OCC_SEASONS(:,:, end + 1) = cat(3, Winds_OCC);
                Winds_STD_SEASONS(:,:, end + 1) = cat(3, Winds_STD);
                Winds_ERR_SEASONS(:,:, end + 1) = cat(3, Winds_ERR);
                
            end           
     
            
        end
        
        clear CMOD CMOD_OCC CMOD_STD CMOD_ERR SeaIce SeaIce_OCC SeaIce_STD SeaIce_ERR Winds Winds_OCC Winds_STD Winds_ERR 
        
    end
end

%%


%%
%%%%%% This was the procedure used to get seasonal ice contour lines %%%%%%

% I basically used climatological averages as the contour lines. 

% Initializing of my variables for the loop below. 


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


%
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
                
                Total_Ice_winter = vertcat(Total_Ice_winter, Ice);
                Total_Lat_Ice_winter = vertcat(Total_Lat_Ice_winter, Lat_Ice); 
                Total_Lon_Ice_winter = vertcat(Total_Lon_Ice_winter, Lon_Ice); 
                                 
                
            elseif j == 2
                
                
                Total_Ice_spring = vertcat(Total_Ice_spring, Ice);
                Total_Lat_Ice_spring = vertcat(Total_Lat_Ice_spring, Lat_Ice); 
                Total_Lon_Ice_spring = vertcat(Total_Lon_Ice_spring, Lon_Ice); 
                                
                
                
            elseif j == 3
                
                
                Total_Ice_summer = vertcat(Total_Ice_summer, Ice);
                Total_Lat_Ice_summer = vertcat(Total_Lat_Ice_summer, Lat_Ice); 
                Total_Lon_Ice_summer = vertcat(Total_Lon_Ice_summer, Lon_Ice); 
                                
               
                
            elseif j == 4
                
                
                Total_Ice_fall = vertcat(Total_Ice_fall, Ice);
                Total_Lat_Ice_fall = vertcat(Total_Lat_Ice_fall, Lat_Ice); 
                Total_Lon_Ice_fall = vertcat(Total_Lon_Ice_fall, Lon_Ice); 
                                
                
            end
     
            
        end
        
        clear CMOD CMOD_OCC CMOD_STD CMOD_ERR 
        
    end
end

%


[Ice_one_degree_winter, Ice_OCC_winter, Ice_STD_winter, Ice_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Ice_winter, Total_Lon_Ice_winter, Total_Ice_winter, step_Lat', step_Lon');
[Ice_one_degree_spring, Ice_OCC_spring, Ice_STD_spring, Ice_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Ice_spring, Total_Lon_Ice_spring, Total_Ice_spring, step_Lat', step_Lon');
[Ice_one_degree_summer, Ice_OCC_summer, Ice_STD_summer, Ice_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Ice_summer, Total_Lon_Ice_summer, Total_Ice_summer, step_Lat', step_Lon');
[Ice_one_degree_fall, Ice_OCC_fall, Ice_STD_fall, Ice_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Ice_fall, Total_Lon_Ice_fall, Total_Ice_fall, step_Lat', step_Lon');

%

winter_contour = Ice_one_degree_winter; 
spring_contour = Ice_one_degree_spring; 
summer_contour = Ice_one_degree_summer;
fall_contour   = Ice_one_degree_fall;

%%
%%%%% Keeping track of my seasons here %%%%%%%%%%%%

t1 = datetime(2006,06,01);
t2 = datetime(2018,12,31);
times = t1:calmonths(1):t2; 
clear t1 t2

times_seasons = times(1) : calmonths(3) : times(end); 


%%

% You are not using this function in the updated version of this
% manuscript, so commenting it out: 
% 
% CMOD_SEASONS_SMOOTHN = zeros(size(CMOD_SEASONS)); 

for i = 1:50
    CMOD_SEASONS_SMOOTHN(:,:,i) = smoothn(CMOD_SEASONS(:,:,i), 'robust') ;
end


CMOD_SEASONS_SMOOTHN(:,:,39) = CMOD_SEASONS(:,:,39); 
CMOD_SEASONS_SMOOTHN(:,:,11) = CMOD_SEASONS(:,:,11); 
CMOD_SEASONS_SMOOTHN(:,:,31) = CMOD_SEASONS(:,:,31); 
CMOD_SEASONS_SMOOTHN(:,:,35) = CMOD_SEASONS(:,:,35); 
CMOD_SEASONS_SMOOTHN(:,:,19) = CMOD_SEASONS(:,:,19); 
CMOD_SEASONS_SMOOTHN(:,:,47) = CMOD_SEASONS(:,:,47);

%%

% Throw out last month of December because it's not a full season: 

CMOD_SEASONS         = CMOD_SEASONS(:,:, 1:50);
% CMOD_SEASONS_SMOOTHN = CMOD_SEASONS_SMOOTHN(:,:, 1:50); 

Winds_SEASONS        = Winds_SEASONS(:,:,1:50); 
SeaIce_SEASONS       = SeaIce_SEASONS(:,:,1:50); 

CHL_SEASONS         = rot90(CHL_SEASONS); % this should be the correct orientation! 
CHL_SEASONS         = CHL_SEASONS(:,:, 1:50); 


%%

% I need to find the mean of each season and then plot out how it's changing from that mean across time. 


%% SeaIce Test with masking values below 0.15

SeaIce_SEASONS_masked = SeaIce_SEASONS; 

SeaIce_SEASONS_masked(SeaIce_SEASONS_masked < 0.15) = NaN;



 
%%

%%%%%%%%% Here I start calculating the spatial correlation between seasons %%%%%%%%% 

% Matrix preallocation 
CMOD_CHL_RHO_MATRIX  = zeros(15,46,4) ; 
CMOD_CHL_PVAL_MATRIX = zeros(15,46,4) ; 

CMOD_ICE_RHO_MATRIX  = zeros(15,46,4) ; 
CMOD_ICE_PVAL_MATRIX = zeros(15,46,4) ; 

CMOD_WIND_RHO_MATRIX  = zeros(15,46,4) ; 
CMOD_WIND_PVAL_MATRIX = zeros(15,46,4) ;  

CHL_ICE_RHO_MATRIX = zeros(15,46,4) ; 
CHL_ICE_PVAL_MATRIX = zeros(15,46,4) ;  

CHL_WIND_RHO_MATRIX = zeros(15,46,4) ; 
CHL_WIND_PVAL_MATRIX = zeros(15,46,4) ; 

WIND_ICE_RHO_MATRIX = zeros(15,46,4) ; 
WIND_ICE_PVAL_MATRIX = zeros(15,46,4) ; 
    
%%
clear x y z k 

for i = 1 : length(step_Lat(1,:))
    % disp(i)
    
    for j = 1 : length(step_Lon(1,:))
        %  disp(j)
        
        for s = 1:4
            
            x = squeeze(CMOD_SEASONS(i,j,s:4:end));
            y = squeeze(CHL_SEASONS(i,j,s:4:end));
            z = squeeze(SeaIce_SEASONS_masked(i,j,s:4:end));
            k = squeeze(Winds_SEASONS(i,j,s:4:end));
            
            % CMOD & CHL
            if sum(isnan(x)) || sum(isnan(y)) <= 3
                [CMOD_CHL_RHO_MATRIX(i, j,s), CMOD_CHL_PVAL_MATRIX(i, j, s)]  = corr(x, y,'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(x)) || sum(isnan(y)) > 3
                isnan([CMOD_CHL_RHO_MATRIX(i, j,s)  , CMOD_CHL_PVAL_MATRIX(i, j, s)]);
            end
            
            
            % CMOD & Ice
            if sum(isnan(x)) || sum(isnan(z)) <= 3
                [CMOD_ICE_RHO_MATRIX(i, j,s), CMOD_ICE_PVAL_MATRIX(i, j,s)]   = corr(x, z, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(x)) || sum(isnan(z)) > 3
                isnan([CMOD_ICE_RHO_MATRIX(i, j,s)  , CMOD_ICE_PVAL_MATRIX(i, j,s)]); 
            end
            
            
            % CMOD & Wind
            if sum(isnan(x)) || sum(isnan(k)) <= 3
                [CMOD_WIND_RHO_MATRIX(i, j,s), CMOD_WIND_PVAL_MATRIX(i, j,s)] = corr(x, k, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(x)) || sum(isnan(k)) > 3
                isnan([CMOD_WIND_RHO_MATRIX(i, j,s) , CMOD_WIND_PVAL_MATRIX(i, j,s)]);
            end
            
            
            % CHL & Ice
            if sum(isnan(z)) || sum(isnan(y)) <= 3
                [CHL_ICE_RHO_MATRIX(i, j,s), CHL_ICE_PVAL_MATRIX(i, j,s)]     = corr(z, y, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(z)) || sum(isnan(y)) > 3
                isnan([CHL_ICE_RHO_MATRIX(i, j,s)   , CHL_ICE_PVAL_MATRIX(i, j,s)]);
            end
            
            
            % CHL & Wind
            if sum(isnan(k)) || sum(isnan(y)) <= 3
                [CHL_WIND_RHO_MATRIX(i,j,s), CHL_WIND_PVAL_MATRIX(i,j,s)]     = corr(k, y, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(k)) || sum(isnan(y)) > 3
                isnan([CHL_WIND_RHO_MATRIX(i,j,s)   , CHL_WIND_PVAL_MATRIX(i,j,s)]);
            end
            
            
            % Wind & Ice
            if sum(isnan(k)) || sum(isnan(z)) <= 3
                [WIND_ICE_RHO_MATRIX(i,j,s), WIND_ICE_PVAL_MATRIX(i,j,s)]     = corr(k, z, 'Type', 'Pearson', 'rows', 'pairwise');
            elseif sum(isnan(k)) || sum(isnan(z)) > 3
                isnan([WIND_ICE_RHO_MATRIX(i,j,s)   , WIND_ICE_PVAL_MATRIX(i,j,s)]);
            end
            
        end
        
    end
    
end










%%%%%%%%% PLOTTING OF SPATIAL CORRELATION FIGURES BELOW %%%%%%%%%%%%%%%%

%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(2,2,1)
plot_BSea_figure_subplot(CMOD_ICE_RHO_MATRIX(:,:,1), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Winter'); 

    hold on;


    [C, h] = m_contour(step_Lon, step_Lat, winter_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2,2,2)
plot_BSea_figure_subplot(CMOD_ICE_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Spring'); 

    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,3)
plot_BSea_figure_subplot(CMOD_ICE_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue,...    
    [-1 1],...
    'Summer'); 
    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,4)
plot_BSea_figure_subplot(CMOD_ICE_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue,...    
    [-1 1],...
    'Fall'); 

    hold on;
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'MAOD & Ice seasonal correlation coefficient', 'FontSize', 18);

%%


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Revisions_Figures

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_MAOD_Ice_Correlation_1_1_newcolormap_masked','-dpdf','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju

%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(2,2,1)
plot_BSea_figure_subplot(CMOD_ICE_PVAL_MATRIX(:,:,1), ...
    step_Lon,...
    step_Lat,...
    [0 0.2],...
    'Winter'); 

    hold on;


    [C, h] = m_contour(step_Lon, step_Lat, winter_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]);  
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2,2,2)
plot_BSea_figure_subplot(CMOD_ICE_PVAL_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    [0 0.2],...
    'Spring'); 

    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,3)
plot_BSea_figure_subplot(CMOD_ICE_PVAL_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    [0 0.2],...
    'Summer'); 
    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor',[ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,4)
plot_BSea_figure_subplot(CMOD_ICE_PVAL_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    [0 0.2],...
    'Fall'); 

    hold on;
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'CMOD & Ice seasonal p-values', 'FontSize', 18);

%%


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end


fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(3, 1, 1 )
plot_BSea_figure_subplot(CMOD_CHL_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Spring'); 

    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1, 2)
plot_BSea_figure_subplot(CMOD_CHL_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue,...   
    [-1 1],...
    'Summer'); 
    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1, 3)
plot_BSea_figure_subplot(CMOD_CHL_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Fall'); 

    hold on;
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 

hp4 = get(subplot(3,1, 3),'Position');
h = colorbar('Position', [hp4(1)+(hp4(3)-0.2)  0.33  0.02  0.4]);
h.FontWeight = 'bold';
h.FontSize = 15;
 hy = ylabel(h, 'MAOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
% hy = ylabel(h, 'CMOD & Ice seasonal correlation coefficient', 'FontSize', 18);

%%



set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_MAOD_CHLA_Correlation_1_1_newcolormap','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju


%%
fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(2,2,1)
plot_BSea_figure_subplot(CMOD_WIND_RHO_MATRIX(:,:,1), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Winter'); 
    hold on; 
    [C, h] = m_contour(step_Lon, step_Lat, winter_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13

subplot(2,2,2)
plot_BSea_figure_subplot(CMOD_WIND_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Spring'); 
    hold on; 


    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]);  
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13

subplot(2,2,3)
plot_BSea_figure_subplot(CMOD_WIND_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Summer'); 
    hold on; 

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13

subplot(2,2,4)
plot_BSea_figure_subplot(CMOD_WIND_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Fall'); 
    hold on; 


    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor',[ 0 0 0]); 
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13

hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
% hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'MAOD & Wind seasonal correlation coefficient', 'FontSize', 18);
% hy = ylabel(h, 'CMOD & Ice seasonal correlation coefficient', 'FontSize', 18);
 

%%



set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_MAOD_Wind_Correlation_1_1_newcolormap.png','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju



%%


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end


fig = figure;clf;


subplot(3, 1, 1)
plot_BSea_figure_subplot(CHL_ICE_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Spring'); 

    hold on ; 

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3, 1, 2)
plot_BSea_figure_subplot(CHL_ICE_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue,...    
    [-1 1],...
    'Summer'); 

    hold on ; 

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3, 1, 3)
plot_BSea_figure_subplot(CHL_ICE_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Fall'); 
    hold on ; 
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

hp4 = get(subplot(3, 1 ,3),'Position');
h = colorbar('Position', [hp4(1)+(hp4(3)-0.25)  0.2  0.02  0.6]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, 'Ice & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);


%%




set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_CHLA_Ice_Correlation_newcolormap_masked','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju



%%


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end


fig = figure;clf;


subplot(2,2,1)
plot_BSea_figure_subplot(WIND_ICE_RHO_MATRIX(:,:,1), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Winter'); 

    hold on;


    [C, h] = m_contour(step_Lon, step_Lat, winter_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]);  
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2, 2, 2)
plot_BSea_figure_subplot(WIND_ICE_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Spring'); 

    hold on ; 

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2, 2, 3)
plot_BSea_figure_subplot(WIND_ICE_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue,...    
    [-1 1],...
    'Summer'); 

    hold on ; 

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2, 2, 4)
plot_BSea_figure_subplot(WIND_ICE_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Fall'); 
    hold on ; 
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
hy = ylabel(h, 'Ice & Wind Speed seasonal correlation coefficient', 'FontSize', 18);


%%




set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Updated_Wind_Ice_Correlation_newcolormap_masked_2.png','-dpng','-r96');       %  *// 300 dpi
%  print(gcf, 'Ju




%% STATISTICAL CALCULATIONS %%


% Matrix preallocation 
randn_CMOD_CHL_RHO_MATRIX  = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 
randn_CMOD_CHL_PVAL_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 

randn_CMOD_ICE_RHO_MATRIX  = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 
randn_CMOD_ICE_PVAL_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 

randn_CMOD_WIND_RHO_MATRIX  = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 
randn_CMOD_WIND_PVAL_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 

randn_CHL_ICE_RHO_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 
randn_CHL_ICE_PVAL_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 

randn_CHL_WIND_RHO_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 
randn_CHL_WIND_PVAL_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 

randn_WIND_ICE_RHO_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 
randn_WIND_ICE_PVAL_MATRIX = zeros([length(CMOD_SEASONS(:, 1, 1)), length(CMOD_SEASONS(1, :, 1))]) ; 
    

%%

CMOD_randn_SEASONS = randn(size(CMOD_SEASONS)); 
CHL_randn_SEASONS = randn(size(CHL_SEASONS)); 
SeaIce_randn_SEASONS = randn(size(SeaIce_SEASONS_masked)); 
Winds_randn_SEASONS = randn(size(Winds_SEASONS)); 


CMOD_randn_SEASONS   = spatialPattern(size(CMOD_SEASONS),-1);
CHL_randn_SEASONS    = spatialPattern(size(CMOD_SEASONS),-1);
SeaIce_randn_SEASONS = spatialPattern(size(CMOD_SEASONS),-1);
Winds_randn_SEASONS  = spatialPattern(size(CMOD_SEASONS),-1);

%%
   
for i = 1 : length(nstep_Lat(1,:))
    % disp(i)
    
    for j = 1 : length(nstep_Lon(1,:))
        %  disp(j)
        
        for s = 1:4
            
            [randn_CMOD_CHL_RHO_MATRIX(i, j,s), randn_CMOD_CHL_PVAL_MATRIX(i, j, s)] = corr(squeeze(CMOD_randn_SEASONS(i,j,s:4:end)),...
                squeeze(CHL_randn_SEASONS(i,j,s:4:end)),'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_CMOD_ICE_RHO_MATRIX(i, j,s), randn_CMOD_ICE_PVAL_MATRIX(i, j,s)] = corr(squeeze(CMOD_randn_SEASONS(i,j,s:4:end)),...
                squeeze(SeaIce_randn_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_CMOD_WIND_RHO_MATRIX(i, j,s), randn_CMOD_WIND_PVAL_MATRIX(i, j,s)] = corr(squeeze(CMOD_randn_SEASONS(i,j,s:4:end)), ...
                squeeze(Winds_randn_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_CHL_ICE_RHO_MATRIX(i, j,s), randn_CHL_ICE_PVAL_MATRIX(i, j,s)] = corr(squeeze(SeaIce_randn_SEASONS(i,j,s:4:end)),...
                squeeze(CHL_randn_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_CHL_WIND_RHO_MATRIX(i,j,s), randn_CHL_WIND_PVAL_MATRIX(i,j,s)] = corr(squeeze(Winds_randn_SEASONS(i,j,s:4:end)),...
                squeeze(CHL_randn_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
            
            [randn_WIND_ICE_RHO_MATRIX(i,j,s), randn_WIND_ICE_PVAL_MATRIX(i,j,s)] = corr(squeeze(Winds_randn_SEASONS(i,j,s:4:end)),...
                squeeze(SeaIce_randn_SEASONS(i,j,s:4:end)), 'Type', 'Pearson', 'rows', 'pairwise');
            
        end
        
    end
    
end


%%
% rand matrix creation with NaNs in same places as NaNs in data matrices

% CMOD & ICe
NAN_randn_CMOD_ICE_RHO_MATRIX = nan(size(CMOD_ICE_RHO_MATRIX));
NAN_randn_CMOD_ICE_RHO_MATRIX(~isnan(CMOD_ICE_RHO_MATRIX)) = randn_CMOD_ICE_RHO_MATRIX(~isnan(CMOD_ICE_RHO_MATRIX)); 


% CHL & ICE
test_CHL_ICE_RHO_MATRIX = CHL_ICE_RHO_MATRIX; 
test_CHL_ICE_RHO_MATRIX(test_CHL_ICE_RHO_MATRIX==0) = NaN; 

% test_CHL_ICE_RHO_MATRIX(:,:,1) = NaN(size(randn_CHL_ICE_RHO_MATRIX(:,:,1))); 

NAN_randn_CHL_ICE_RHO_MATRIX = nan(size(randn_CHL_ICE_RHO_MATRIX)); 
NAN_randn_CHL_ICE_RHO_MATRIX(~isnan(test_CHL_ICE_RHO_MATRIX)) = randn_CHL_ICE_RHO_MATRIX(~isnan(test_CHL_ICE_RHO_MATRIX)); 


% CMOD & CHL
test_CMOD_CHL_RHO_MATRIX = CMOD_CHL_RHO_MATRIX; 
test_CMOD_CHL_RHO_MATRIX(:,:,1) = NaN(size(CMOD_CHL_RHO_MATRIX(:,:,1))); 

NAN_randn_CMOD_CHL_RHO_MATRIX = nan(size(test_CMOD_CHL_RHO_MATRIX)); 
NAN_randn_CMOD_CHL_RHO_MATRIX(~isnan(test_CMOD_CHL_RHO_MATRIX)) = randn_CMOD_CHL_RHO_MATRIX(~isnan(test_CMOD_CHL_RHO_MATRIX)); 


% CMOD & WIND

NAN_randn_CMOD_WIND_RHO_MATRIX = nan(size(CMOD_WIND_RHO_MATRIX)); 
NAN_randn_CMOD_WIND_RHO_MATRIX(~isnan(CMOD_WIND_RHO_MATRIX)) = randn_CMOD_WIND_RHO_MATRIX(~isnan(CMOD_WIND_RHO_MATRIX));


% trying to keep NaNs in same position as where ice was removed so that
% random matrices statistics can be computed
nx = length(CMOD_ICE_RHO_MATRIX(:,1,1));
ny = length(CMOD_ICE_RHO_MATRIX(1,:,1));
nz = length(CMOD_ICE_RHO_MATRIX(1,1,:));

CMOD_ICE_DATA_vector = reshape(CMOD_ICE_RHO_MATRIX, [nx*ny*nz 1]);
CMOD_CHLA_DATA_vector = reshape(CMOD_CHL_RHO_MATRIX, [nx*ny*nz 1]);
CMOD_WIND_DATA_vector = reshape(CMOD_WIND_RHO_MATRIX, [nx*ny*nz 1]); 
ICE_CHLA_DATA_vector = reshape(test_CHL_ICE_RHO_MATRIX, [nx*ny*nz 1]); 

% RAND_vector1 = reshape(randn_CMOD_CHL_RHO_MATRIX, [nx*ny*nz 1] ) ;
RAND_vector1 = reshape(NAN_randn_CMOD_ICE_RHO_MATRIX, [nx*ny*nz 1]); 
RAND_vector2 = reshape(NAN_randn_CHL_ICE_RHO_MATRIX, [nx*ny*nz 1]); 
RAND_vector3 = reshape(NAN_randn_CMOD_CHL_RHO_MATRIX, [nx*ny*nz 1]); 
RAND_vector4 = reshape(NAN_randn_CMOD_WIND_RHO_MATRIX, [nx*ny*nz 1]); 

% Check assumptions before conducting ANOVA / Kruskal - Wallis

% If h = 1, this indicates the rejection of the null hypothesis at the Alpha significance level.
% If h = 0, this indicates a failure to reject the null hypothesis at the Alpha significance level.

% returns a test decision for the null hypothesis that 
% the data in vector x is from a population with a normal distribution, 
% using the Anderson-Darling test

% For ANOVA, samples need to come from a normal distribution. if h = 0, you
% fail to reject the null that the samples come from a normal distribution,
% and can continue with the ANOVA. 

h_rand1     = adtest(RAND_vector1); 
h_rand2     = adtest(RAND_vector2); 
h_rand3     = adtest(RAND_vector3);
h_rand4     = adtest(RAND_vector4);
h_CMOD_wind = adtest(CMOD_WIND_DATA_vector); 
h_CMOD_ice  = adtest(CMOD_ICE_DATA_vector); 
h_CMOD_chla = adtest(CMOD_CHLA_DATA_vector); 
h_ice_chla  = adtest(ICE_CHLA_DATA_vector); 

%%

DATA = cat(2,RAND_vector1, RAND_vector2, RAND_vector3, RAND_vector4, CMOD_ICE_DATA_vector,...
    CMOD_CHLA_DATA_vector, CMOD_WIND_DATA_vector, ICE_CHLA_DATA_vector);


groupnames = {'random MAOD & ice', 'random ice & chl-a', 'random MAOD & chl-a','random MAOD & wind',...
    'MAOD & ice', 'MAOD & chl-a', 'MAOD & wind', 'ice & chl-a'}; 


[p_anova,tbl_anova, stats_anova]                         = anova1(DATA, groupnames);

% p_kruskalwallis = kruskalwallis(DATA)

% 
% Fstat  = tbl{2,5};
% pvalue = tbl{2,6};

%%

fig = figure; clf; 
[c_anova, m_anova, h_anova, nms_anova]  = multcompare(stats_anova); 
set(gca, 'FontSize', 16);
title('Multiple Comparison Plot, ANOVA'); 
grid off

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 


%%
print(fig,'Multiple_Comparison_Plot_ANOVA.png','-dpng','-r96');       %  *// 300 dpi
writematrix(c_anova,'c_anova.xls')

%%

fig = figure;clf; 
[p_kruskalwallis, tbl_kruskalwallis, stats_kruskalwallis] = kruskalwallis(DATA, groupnames);
set(gca, 'FontSize', 22);
xtickangle(30)

fig = figure; clf;
[c_kruskalwallis, m_kruskalwallis,...
    h_kruskalwallis, nms_kruskalwallis] = multcompare(stats_kruskalwallis); 

set(gca, 'FontSize', 22);
title('Multiple Comparison Plot, Kruskal-Wallis'); 
xlabel('Testing whether mean ranks of groups are significantly different')
% ylabel('Groups')
% grid on
%%
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'Multiple_Comparison_Plot_KruskalWallis','-dpng','-r300');       %  *// 300 dpi


%%
% writematrix(c_kruskalwallis,'c_kruskalwallis.xls')
% c returns a matrix of the pairwise comparison results from a multiple
% comparison test using information contained in stats structure


%%
% 
% 
% DATA = cat(2,RAND_vector1, RAND_vector2, RAND_vector3, RAND_vector4, CMOD_ICE_DATA_vector,...
%     CMOD_CHLA_DATA_vector, CMOD_WIND_DATA_vector, ICE_CHLA_DATA_vector);
% 
% groupnames = {'random MOD & ice', 'random ice & chl-a', 'random MOD & chl-a','random MOD & wind', 'MOD & ice', 'MOD & chl-a', 'MOD & wind', 'ice & chl-a'}; 


[p,tbl,stats] = anova1(DATA, groupnames);
[p,tbl,stats] = kruskalwallis(DATA, groupnames);

Fstat  = tbl{2,5};
pvalue = tbl{2,6};

[c, m, h, nms] = multcompare(stats); 


%% % Plotting of spatial statistical distributions given random numbers 


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(2,2,1)
plot_BSea_figure_subplot(NAN_randn_CMOD_ICE_RHO_MATRIX(:,:,1), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Winter'); 

    hold on;


    [C, h] = m_contour(step_Lon, step_Lat, winter_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65 ; 
%     clabel(C, h,'fontsize',13);
%     clegendm(C,h,'m')

subplot(2,2,2)
plot_BSea_figure_subplot(NAN_randn_CMOD_ICE_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Spring'); 

    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,3)
plot_BSea_figure_subplot(NAN_randn_CMOD_ICE_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Summer'); 
    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,4)
plot_BSea_figure_subplot(NAN_randn_CMOD_ICE_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue,...
    [-1 1],...
    'Fall'); 

    hold on;
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 
    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
%  hy = ylabel(h, 'CMOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
 hy = ylabel(h, 'Random number MAOD & Ice seasonal correlation coefficient', 'FontSize', 18);





%%

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'RandomNumber_Seasonal_Spatial_Correlation_MAOD_ICE_1_1_masked_newcolormap.png','-dpng','-r96');       %  *// 300 dpi
%  print(gcf, 'Ju





%%


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end


fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(3, 1, 1 )
plot_BSea_figure_subplot(NAN_randn_CHL_ICE_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Spring'); 

    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1, 2)
plot_BSea_figure_subplot(NAN_randn_CHL_ICE_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Summer'); 
    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1, 3)
plot_BSea_figure_subplot(NAN_randn_CHL_ICE_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Fall'); 

    hold on;
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 

hp4 = get(subplot(3,1, 3),'Position');
h = colorbar('Position', [hp4(1)+(hp4(3)-0.2)  0.33  0.02  0.4]);
h.FontWeight = 'bold';
h.FontSize = 15;
 hy = ylabel(h, 'Random Number ice & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
% hy = ylabel(h, 'CMOD & Ice seasonal correlation coefficient', 'FontSize', 18);



%%

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'RandomNumber_Seasonal_Spatial_Correlation_ice_chla_1_1_masked_newcolormap','-dpng','-r300');       %  *// 300 dpi
%  



%%


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end


fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(3, 1, 1 )
plot_BSea_figure_subplot(NAN_randn_CMOD_CHL_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Spring'); 

    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1, 2)
plot_BSea_figure_subplot(NAN_randn_CMOD_CHL_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Summer'); 
    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(3,1, 3)
plot_BSea_figure_subplot(NAN_randn_CMOD_CHL_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Fall'); 

    hold on;
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 

hp4 = get(subplot(3,1, 3),'Position');
h = colorbar('Position', [hp4(1)+(hp4(3)-0.2)  0.33  0.02  0.4]);
h.FontWeight = 'bold';
h.FontSize = 15;
 hy = ylabel(h, 'Random Number MAOD & chl-{\ita} seasonal correlation coefficient', 'FontSize', 18);
% hy = ylabel(h, 'CMOD & Ice seasonal correlation coefficient', 'FontSize', 18);





%%

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'RandomNumber_Seasonal_Spatial_Correlation_MAOD_chla_1_1_masked_newcolormap','-dpng','-r300');       %  *// 300 dpi
%  


%%


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.08 0.08], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end


fig = figure;clf;

% chl_season = rot90(Chl_a_seasons_one_degree_res); 

% test_2 = smoothn(test, 'robust'); 
subplot(2,2,1)
plot_BSea_figure_subplot(NAN_randn_CMOD_WIND_RHO_MATRIX(:,:,1), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Winter'); 

    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, winter_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,2)
plot_BSea_figure_subplot(NAN_randn_CMOD_WIND_RHO_MATRIX(:,:,2), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Spring'); 
    hold on;

    [C, h] = m_contour(step_Lon, step_Lat, spring_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

subplot(2,2,3)
plot_BSea_figure_subplot(NAN_randn_CMOD_WIND_RHO_MATRIX(:,:,3), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Summer'); 

    hold on;
    [C, h] = m_contour(step_Lon, step_Lat, summer_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 subplot(2,2,4)
plot_BSea_figure_subplot(NAN_randn_CMOD_WIND_RHO_MATRIX(:,:,4), ...
    step_Lon,...
    step_Lat,...
    redblue, ...
    [-1 1],...
    'Fall'); 

    hold on;
    [C, h] = m_contour(step_Lon, step_Lat, fall_contour, 'linewi', 5.5, 'LineColor', [ 0 0 0]); 
    h.LevelList = 0.65  ; 
%     clabel(C, h,'fontsize',13);

 

    
hp4 = get(subplot(2,2,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.15  0.02  hp4(2)+hp4(3) * 1.3]);
h.FontWeight = 'bold';
h.FontSize = 15;
 hy = ylabel(h, 'Random Number MAOD & Wind Speed seasonal correlation coefficient', 'FontSize', 18);
% hy = ylabel(h, 'CMOD & Ice seasonal correlation coefficient', 'FontSize', 18);





%%

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(fig,'RandomNumber_Seasonal_Spatial_Correlation_MAOD_wind_1_1_masked_newcolormap','-dpng','-r300');       %  *// 300 dpi
%  

