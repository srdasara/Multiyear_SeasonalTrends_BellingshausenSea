%% Read MAOD, Construct into timetable, save out in 
% /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/CALIPSO_LID_L2_05kmAPro-Standard-V4-20/UPDATED_Calipso_2006_2018_CMOD_Night_CloudFree_AdjustedAltRange

% can these push these into monthly averaged time series
% climatologies
% spearman's rho 
%%

load('Total_Surface_Type.mat')
load('Calipso_2006_2018_Surface_Elevation_Statistics.mat')


%%

Total_Profile_Time_Surface = Total_Profile_Time ; 


Total_Profile_Time_New_Surface =datetime(Total_Profile_Time_Surface,'ConvertFrom',...
        'epochtime','Epoch','1993-01-01');



    Total_table_surface = table(Total_Profile_Time_New_Surface,Total_Surface_Elevation_Statistics); 

    Total_table_surface                = sortrows(Total_table_surface,'Total_Profile_Time_New_Surface','ascend'); % sort values with increasing time duration
    Total_timetable_surface            = table2timetable(Total_table_surface); % make table into a timetable
    Total_Surface_Elevation_Statistics = Total_timetable_surface.Total_Surface_Elevation_Statistics;
   
    clear Total_table_surface Total_timetable_surface Total_Profile_Time Total_Profile_Time_New_Surface Total_Profile_Time Total_Profile_Time_Surface
    
%     isequal(Total_Profile_Time_New_Surface,Total_Profile_Time_New )
% test to see if these arrays are the same;
% Looks like microseconds are off. 

%% 0.0977KM TO 2.0137KM if not passing any filters, otherwise 
% 0.03779 km to 2.0137 km, if surface elev statistics permit it. 

% these correspond to altitude bins 358 to 391. 

% first load these
cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication

load('Total_EC_532.mat')
load('Total_Day_Night_Flag.mat')
load('Total_COD_Cloud.mat')
load('Total_altitudes.mat')
load('Total_Longitude.mat')
load('Total_Latitude.mat') 
load('Total_Profile_Time_New.mat')
load('Total_Surface_532_Integrated_Depolarization_Ratio.mat')
% load('Total_Profile_Time.mat') 
load('Total_windamsrMF.mat') 

   
% cd /Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication


%% Filter to only include Extinction Coefficients that are cloud free and nighttime only. 
% Total_COD_Cloud == 0 , DayNightFlag == 1

Cloud_Free                                          = Total_COD_Cloud(:, 1) == 0; 
Total_Latitude_Cloud_Free                           = Total_Latitude(Cloud_Free); 
Total_Longitude_Cloud_Free                          = Total_Longitude(Cloud_Free); 
Total_Profile_Time_New_Cloud_Free                   = Total_Profile_Time_New(Cloud_Free);
Total_Surface_Elevation_Statistics_Cloud_Free       = Total_Surface_Elevation_Statistics(Cloud_Free);
Total_Surface_Type_Cloud_Free                       = Total_Surface_Type(Cloud_Free); 
Total_EC_532_Cloud_Free                             = Total_EC_532(Cloud_Free, :) ; 


Total_Day_Night_Flag_Cloud_Free                     = Total_Day_Night_Flag(Cloud_Free,:); 
Night_Cloud_Free                                    = Total_Day_Night_Flag_Cloud_Free(:,1) == 1 ; 
Total_EC_532_Cloud_Free_Night                       = Total_EC_532_Cloud_Free(Night_Cloud_Free,:); 
Total_Latitude_Night_Cloud_Free                     = Total_Latitude_Cloud_Free(Night_Cloud_Free); 
Total_Longitude_Night_Cloud_Free                    = Total_Longitude_Cloud_Free(Night_Cloud_Free); 
Total_Profile_Time_New_Night_Cloud_Free             = Total_Profile_Time_New_Cloud_Free(Night_Cloud_Free); 
Total_EC_532_Night_Cloud_Free                       = Total_EC_532_Cloud_Free(Night_Cloud_Free, :); 
Total_Surface_Elevation_Statistics_Night_Cloud_Free = Total_Surface_Elevation_Statistics_Cloud_Free(Night_Cloud_Free);
Total_Surface_Type_Night_Cloud_Free                 = Total_Surface_Type_Cloud_Free(Night_Cloud_Free);




% 0.03779 km to 2.0137 km, if surface elev statistics & Surface type permit it. 

% Good Surface Type and Surface Elev Stats:

Total_Surface_Good                  = Total_Surface_Type_Night_Cloud_Free(:,1) == 17 & Total_Surface_Elevation_Statistics_Night_Cloud_Free(:,1) == 0 ; 
Total_EC_532_Surface_Good           = Total_EC_532_Night_Cloud_Free(Total_Surface_Good,:); 
Total_Profile_Time_New_Surface_Good = Total_Profile_Time_New_Night_Cloud_Free(Total_Surface_Good); 
Total_Latitude_Surface_Good         = Total_Latitude_Night_Cloud_Free(Total_Surface_Good); 
Total_Longitude_Surface_Good        = Total_Longitude_Night_Cloud_Free(Total_Surface_Good);

Total_EC_532_Surface_Bad            = Total_EC_532_Night_Cloud_Free(~Total_Surface_Good, :); 
Total_Profile_Time_New_Surface_Bad  = Total_Profile_Time_New_Night_Cloud_Free(~Total_Surface_Good);
Total_Latitude_Surface_Bad          = Total_Latitude_Night_Cloud_Free(~Total_Surface_Good);
Total_Longitude_Surface_Bad         = Total_Longitude_Night_Cloud_Free(~Total_Surface_Good); 


Total_EC_532_Surface_Good_adjusted_alt = Total_EC_532_Surface_Good(:, 358:391) ; % 0.03779 km to 2.0137 km, if surface elev statistics permit it. 
Total_adjusted_alt_Surface_Good        = Total_altitudes(358:391, :) ; 

Total_EC_532_Surface_Bad_adjusted_alt = Total_EC_532_Surface_Bad(:, 358:390); 
Total_adjusted_alt_Surface_Bad        = Total_altitudes(358:390, :);
% Total_Night_Flag_Cloud_Free           = Total_Day_Night_Flag_Cloud_Free(Night_Cloud_Free); 

     
%%                              
    
% Convert NaNs to 0 for trapz function:

 Total_EC_532_Surface_Good_adjusted_alt(isnan(Total_EC_532_Surface_Good_adjusted_alt)) = 0 ; % 0 is clear air, NaN has been filtered out by quality screening
 Total_EC_532_Surface_Bad_adjusted_alt(isnan(Total_EC_532_Surface_Bad_adjusted_alt)) = 0; 
 %  Total_EC_532_Night_Cloud_Free(Total_EC_532_Night_Cloud_Free == 0) = NaN ; % Converting all zeros in sigma to NaNs. 

 % to keep NaN or not keep NaN?
 %%
clear CMOD_Surface_Good
CMOD_Surface_Good = zeros(length(Total_EC_532_Surface_Good_adjusted_alt(:,1)), 1); 

for i = 1:length(Total_EC_532_Surface_Good_adjusted_alt(:,1))
%     disp(i)
    CMOD_Surface_Good(i) = -1 .* (trapz(Total_adjusted_alt_Surface_Good, Total_EC_532_Surface_Good_adjusted_alt(i,:))) ; 
    % -1 in equation above was to flip in consideration of the fact that altitudes start from 2.0137 km and end at 0.0977km
end

% save('CMOD.mat', 'CMOD') ; 

clear CMOD_Surface_Bad
CMOD_Surface_Bad = zeros(length(Total_EC_532_Surface_Bad_adjusted_alt(:,1)), 1); 

for i = 1:length(Total_EC_532_Surface_Bad_adjusted_alt(:,1)) 
   
    CMOD_Surface_Bad(i) = -1 .* (trapz(Total_adjusted_alt_Surface_Bad, Total_EC_532_Surface_Bad_adjusted_alt(i,:)));
    % -1 in equation above was to flip in consideration of the fact that altitudes start from 2.0137 km and end at 0.0977km

end

    %% Make a time table with all of these values, 3 separate ones: CMOD, Ice, and Winds
    Total_Profile_Time_New_Surface = Total_Profile_Time_New_Surface_Good;
    Total_Latitude_Surface = Total_Latitude_Surface_Good; 
    Total_Longitude_Surface = Total_Longitude_Surface_Good; 
    CMOD_Surface = CMOD_Surface_Good;
    
    Total_table_CMOD_Surface_Good = table(Total_Profile_Time_New_Surface,...
        Total_Latitude_Surface,...
        Total_Longitude_Surface,...
        CMOD_Surface); 

    Total_table_CMOD_Surface_Good             = sortrows(Total_table_CMOD_Surface_Good,'Total_Profile_Time_New_Surface','ascend'); % sort values with increasing time duration
    Total_timetable_CMOD_Surface_Good        = table2timetable(Total_table_CMOD_Surface_Good); % make table into a timetable

    
    
    
    Total_Profile_Time_New_Surface = Total_Profile_Time_New_Surface_Bad; 
    Total_Latitude_Surface = Total_Latitude_Surface_Bad;
    Total_Longitude_Surface = Total_Longitude_Surface_Bad; 
    CMOD_Surface = CMOD_Surface_Bad; 
    
    Total_table_CMOD_Surface_Bad = table(Total_Profile_Time_New_Surface,...
        Total_Latitude_Surface,...
        Total_Longitude_Surface,...
        CMOD_Surface); 
    
    Total_table_CMOD_Surface_Bad = sortrows(Total_table_CMOD_Surface_Bad, 'Total_Profile_Time_New_Surface', 'ascend');
    Total_timetable_CMOD_Surface_Bad       = table2timetable(Total_table_CMOD_Surface_Bad); % make table into a timetable

    Total_timetable_CMOD_Surface = [Total_timetable_CMOD_Surface_Good ; Total_timetable_CMOD_Surface_Bad] ; 
    Total_timetable_CMOD_Surface = sortrows(Total_timetable_CMOD_Surface, 'Total_Profile_Time_New_Surface', 'ascend');
    
    
   
     save('Total_timetable_CMOD_Surface.mat', 'Total_timetable_CMOD_Surface', '-v7.3')
    %% I first have to filter for only good values of Ice 
    
    bad_Ice_values = Total_Surface_532_Integrated_Depolarization_Ratio <= -0.2 | Total_Surface_532_Integrated_Depolarization_Ratio > 1.2;
    Total_Surface_532_Integrated_Depolarization_Ratio(bad_Ice_values) = NaN; % I set these bad values to NaNs so I can easily index and remove them
    
    nan_ice        = isnan(Total_Surface_532_Integrated_Depolarization_Ratio(:,1));
    Total_Surface_532_Integrated_Depolarization_Ratio      = Total_Surface_532_Integrated_Depolarization_Ratio(~nan_ice) ;
    Total_Latitude_Ice  = Total_Latitude(~nan_ice);
    Total_Longitude_Ice  = Total_Longitude(~nan_ice);
    Total_Profile_Time_New_Ice = Total_Profile_Time_New(~nan_ice); 
    
    Total_table_Depol_Ratio = table(Total_Profile_Time_New_Ice,...
        Total_Latitude_Ice,...
        Total_Longitude_Ice,...
        Total_Surface_532_Integrated_Depolarization_Ratio);
    
    Total_table_Depol_Ratio = sortrows(Total_table_Depol_Ratio, 'Total_Profile_Time_New_Ice', 'ascend');
    
    Total_timetable_Depol_Ratio = table2timetable(Total_table_Depol_Ratio);
    save('Total_timetable_Depol_Ratio.mat', 'Total_timetable_Depol_Ratio', '-v7.3') 
    
    %%
    
    bad_Wind_values = Total_windamsrMF <= 0| Total_windamsrMF > 50;
    Total_windamsrMF(bad_Wind_values) = NaN; % I set these bad values to NaNs so I can easily index and remove them

    nan_wind        = isnan(Total_windamsrMF(:,1));
    Total_windamsrMF      = Total_windamsrMF(~nan_wind) ;
    Total_Latitude_Wind  = Total_Latitude(~nan_wind);
    Total_Longitude_Wind  = Total_Longitude(~nan_wind);
    Total_Profile_Time_New_Wind = Total_Profile_Time_New(~nan_wind); 
    
    
    Total_table_amsrmf = table(Total_Profile_Time_New_Wind,... 
        Total_Latitude_Wind,...
        Total_Longitude_Wind,... 
        Total_windamsrMF);
    
    Total_table_amsrmf = sortrows(Total_table_amsrmf, 'Total_Profile_Time_New_Wind', 'ascend'); 
    Total_timetable_amsrmf = table2timetable(Total_table_amsrmf); 
    save('Total_timetable_amsrmf.mat', 'Total_timetable_amsrmf', '-v7.3') 
    
    %%
    
    timetable_CMOD_monthly_avg    = retime(Total_timetable_CMOD_Surface, 'monthly', @nanmean); 
    CMOD_Monthly_avg_Surface = timetable_CMOD_monthly_avg.CMOD_Surface; 
    CMOD_Time_Months_Surface = timetable_CMOD_monthly_avg.Total_Profile_Time_New_Surface; 
    CMOD_Lat_Months_Surface  = timetable_CMOD_monthly_avg.Total_Latitude_Surface;
    CMOD_Lon_Months_Surface  = timetable_CMOD_monthly_avg.Total_Longitude_Surface;
        
    timetable_Depol_Ratio_monthly_avg = retime(Total_timetable_Depol_Ratio, 'monthly', @nanmean); 
    Depol_Ratio_Monthly_avg = timetable_Depol_Ratio_monthly_avg.Total_Surface_532_Integrated_Depolarization_Ratio;
    Depol_Ratio_Time_Months = timetable_Depol_Ratio_monthly_avg.Total_Profile_Time_New_Ice; 
    Depol_Ratio_Lat_Months = timetable_Depol_Ratio_monthly_avg.Total_Latitude_Ice; 
    Depol_Ratio_Lon_Months = timetable_Depol_Ratio_monthly_avg.Total_Longitude_Ice; 
    
    timetable_amsrmf_monthly_avg = retime(Total_timetable_amsrmf, 'monthly', @nanmean); 
    amsrmf_Monthly_avg = timetable_amsrmf_monthly_avg.Total_windamsrMF; 
    amsrmf_Time_Months = timetable_amsrmf_monthly_avg.Total_Profile_Time_New_Wind; 
    amsrmf_Lat_Months = timetable_amsrmf_monthly_avg.Total_Latitude_Wind;
    amsrmf_Lon_Months = timetable_amsrmf_monthly_avg.Total_Longitude_Wind; 
    
    save('CMOD_Monthly_avg_Vars_Surface.mat', ...
        'CMOD_Monthly_avg_Surface',...
        'CMOD_Time_Months_Surface',...
        'CMOD_Lat_Months_Surface',...
        'CMOD_Lon_Months_Surface',...
        '-v7.3') 
    
    save('Depol_Ratio_Monthly_avg_Vars.mat',...
        'Depol_Ratio_Monthly_avg',...
        'Depol_Ratio_Time_Months',...
        'Depol_Ratio_Lat_Months',...
        'Depol_Ratio_Lon_Months',...
        '-v7.3') 
    
    save('amsrmf_Monthly_avg_Vars.mat',...
        'amsrmf_Monthly_avg',...
        'amsrmf_Time_Months',...
        'amsrmf_Lat_Months',...
        'amsrmf_Lon_Months',...
        '-v7.3')
    
    %%
    
    % I went back in to calculation the standard deviation
    % loaded the timetables from
    % '/Users/srishtidasarathy/Documents/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication'

    %      I created this function in a .m file
%     function y = std_timetable(x) 
%     y = std(x, 'omitnan') ;
%     end

%%

% I've constructed all of the standard deviation values into timetables. 
% Draft code also inside
        CMOD_no_zeros = Total_timetable_CMOD_Surface.CMOD_Surface; 
        CMOD_no_zeros(CMOD_no_zeros==0) = nan;
        Total_timetable_CMOD_test = addvars(Total_timetable_CMOD_Surface,CMOD_no_zeros); 
        
        CMOD_std = retime(Total_timetable_CMOD_test, 'monthly', @nanstd);
        
        
        
        CMOD_std = retime(Total_timetable_CMOD_test, 'monthly', @nanstd); 
        CMOD_mean = retime(Total_timetable_CMOD_test, 'monthly', @nanmean); 
        CMOD_mean_absol_dev = retime(Total_timetable_CMOD_test, 'monthly', @mad); 
        CMOD_median_absol_dev = retime(Total_timetable_CMOD_test, 'monthly', @mad_median); 
        
        
        Ice_median_absol_dev = retime(Total_timetable_Depol_Ratio, 'monthly', @mad); 
        Ice_mean = retime(Total_timetable_Depol_Ratio, 'monthly', @nanmean); 
        
        
%         CMOD_std = CMOD_std.CMOD;
        Wind_std = retime(Total_timetable_amsrmf, 'monthly', @std_timetable); 
%         Wind_std = Wind_std.Total_windamsrMF; 
        Ice_std = retime(Total_timetable_Depol_Ratio, 'monthly', @std_timetable); 
%         Ice_std = Ice_std.Total_Surface_532_Integrated_Depolarization_Ratio; 
    
% 
%         Chl_std_test = std(Master_chl_a, 0, [1 2], 'omitnan'); 
%         Chl_std = squeeze(Chl_std_test);  
%         
%         Total_table_chl_a_std = table(times, Chl_std);    
%         timetable_chl_a_std = sortrows(Total_table_chl_a_std, 'times', 'ascend'); 
    
        %%
        % I did this after loading Master_chlor_a_monthly_full_res.mat
        % Here is the mean absolute deviation across every page of
        % Master_chl_a 
        
        
        t1 = datetime(2006,06,01);
        t2 = datetime(2018,12,31);
        times = t1:calmonths(1):t2; 
        times = times';
    
    
        
        
        for i = 1:151 
            Total_chl_monthly_2(i) = nanmean(Master_chl_a(:,:,i), [1 2]); 
        end

        Total_chl_a_monthly = Total_chl_monthly_2'; 


        Total_chl_a_mad = mad(Master_chl_a, 0, [1 2]);
        
        Total_chl_a_std = nanstd(Master_chl_a, 0, [ 1 2]); 
        
        Total_chl_a_mad = squeeze(Total_chl_a_mad); 
        Total_chl_a_std = squeeze(Total_chl_a_std); 
        
        Total_timetable_chl_a_monthly_plus_mad = timetable(times, Total_chl_a_monthly,...
            Total_chl_a_mad, Total_chl_a_std);

        
        %%
        
        % Standard deviation test: 
         A = [4 -5 1 2 3 5 -9 1 7];
         std_A = std(A); 
         mean_A = mean(A); 
        
        
        % June 2006 example
        
        S = timerange('06/01/2006','07/01/2006');
        timetable_test = Total_timetable_CMOD(S,:); 
        
        CMOD_max_five = maxk(timetable_test.CMOD, 10); 
        edges = [0 0:0.0001:0.01 0.01];
        histogram(timetable_test.CMOD, edges) 
    
        
        
        
        
        
        
        
        
        
    