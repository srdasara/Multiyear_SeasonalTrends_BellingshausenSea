
%% Attempting to calculate the number of extinction coefficients that pass quality control: 
cd '/Users/srishtidasarathy/Documents/Documents - Srishti’s MacBook Pro /Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication'
load vars_of_interest.mat % this is me taking extinction coefficient, COD_Cloud, and Day_Night_Flag arrays
% This file is all of data from 2006 to end of 2018.

% Parse Total_EC_532 to correspond to the same length as earlier
% Bellingshausen Sea subset, i.e., only until end of 2018 in length. Also
% true for COD and Day Night Flag 

load('Total_EC_532.mat')
load('Total_COD_Cloud.mat')
load('Total_Day_Night_Flag.mat')
load('Total_Profile_Time_New.mat')


%%

% I have variables 'Total_EC_532' and 'Total_COD_Cloud' and
% 'Total_Day_Night_Flag'
% where NaNs in Total_EC_532 are those extinction coefficients 
% that did not pass quality control:

% I can first filter my Total_EC_532 to only include columns that are cloud
% free & Nighttime

Night_Cloud_Free              = (Total_COD_Cloud(:, 1) == 0) & (Total_Day_Night_Flag(:,1) == 1); % where there is no cloud presence and corresponds to nighttime. 
Total_EC_532_Night_Cloud_Free = Total_EC_532(Night_Cloud_Free, :); 

% Now I'm already left with 334453 columns out of the total 6386088

% Of this remaining array, I can count the number of NaNs, where NaNs are those extinction coefficients 
% that did not pass quality control:

Total_NaN = sum(isnan(Total_EC_532_Night_Cloud_Free), 'all') ;
disp(Total_NaN);

% Now I can count those that are good values, where they only correspond to qc'd clean
% marine extinction coefficients (clean air extinction coefficients are
% also excluded).
good_values = sum(~isnan(Total_EC_532_Night_Cloud_Free) & Total_EC_532_Night_Cloud_Free ~= 0, 'all');
disp(good_values); 

percent_cleanmarine = ( good_values / (Total_NaN + good_values)) * 100;
disp( percent_cleanmarine )
% so of nighttime cloudfree extinction coefficients, ~44.5% of them
% correspond to clean marine values. This does not include out of the
% fraction that are also reported as clear air. 

%% 
% But when we consider the columns that were removed because of nighttime
% and cloud free (this section is more just a guess)

nighttime_cloudfree_columns = length(Total_EC_532_Night_Cloud_Free(:,1));
daytime_cloudy_columns      = length(Total_EC_532(:,1)) - nighttime_cloudfree_columns;
% 6051635

daytime_cloudy_ECs = daytime_cloudy_columns * 399; % full size of matrix of possible extinction coefficients
percent_cleanmarine_night_cloudfree = (good_values / (daytime_cloudy_ECs + Total_NaN + good_values))  * 100; 

% 0.0815 percent (!!!) of extinction coefficients make the cut. This is a lower
% estimate since not all altitude "blocks" in daytime_cloudy_ECs will
% actually have reportings of Extinction coefficients. 




