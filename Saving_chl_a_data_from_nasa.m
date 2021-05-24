


%% Here is how I read chl-a data and saved out the variables needed! 


% function [stat] = Read_hdf_from_server(year_run, month_run)
% 
% Year          = {'2006', '2007', '2008', '2009', '2010', '2011', '2012','2013', '2014', '2015','2016','2017', '2018', '2019'};
% Month         =  {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10','11', '12'};
cd /Users/srishtidasarathy/Documents/Bowman/New_download_CHL_A_Data_03_2020

Chl_data_dir = '/Users/srishtidasarathy/Documents/Bowman/New_download_CHL_A_Data_03_2020/' ;


nc_files     = dir([Chl_data_dir '*4km.nc']);

% file_with_data_in_range = 0;

% Bellingshausen Sea Spatial Subset
 LatLims  = [-75  -60];
 LonLims  = [-100 -54];

%%
  Master_count_n   = [];
  Master_filename  = [];
        
% data
  Master_Latitude  = [];
  Master_Longitude = [];
  Master_chl_a     = [];
  inRange_Lon      = [];
    
  %%
  for n = 1:length(nc_files)
      
      disp(n)
      myFileName   = nc_files(n,:).name;
      fullFileName = fullfile(Chl_data_dir, myFileName);
      
      product_name        = 'lat';
      Latitude            = ncread( fullFileName , product_name );
      
      inRange_Lat = Latitude >= LatLims(1) & Latitude < LatLims(end); 
      
      product_name        = 'lon';
      Longitude           = ncread( fullFileName , product_name );
      
      inRange_Lon = Longitude >= LonLims(1) & Longitude < LonLims(end) ; 

      % sum(inRange) >= 1 so we have data
      %tracking number of file with data in range
%       file_with_data_in_range         = file_with_data_in_range +1;
      
      product_name        = 'chlor_a';
      chlor_a             = ncread( fullFileName, product_name) ;
      
      
      Latitude_subset   = Latitude(inRange_Lat) ;
      Longitude_subset  = Longitude(inRange_Lon) ;
      chl_a_subset      = chlor_a(inRange_Lon, inRange_Lat);
      time              = nc_files(n).name(2:15);  % filename has year followed by number of days from first day of year 
                      % to year then end of number of days from first day of year
     
      if  n == 1          
          % counter
          Master_count_n   = n ;
          Master_filename  = fullFileName;
          Master_time      = time ; 
          Master_chl_a     = chl_a_subset ;          
      else
          % counter          
          Master_count_n   = cat(1, Master_count_n, n);
          Master_filename  = cat(1, Master_filename, fullFileName);
          Master_time      = cat(1, Master_time, time) ; 
          Master_chl_a     = cat(3, Master_chl_a, chl_a_subset);    %% concatenated along the third dimension
      end    
      
  end
  
    %%
% Master_chl_a_mean = mean(Master_chl_a, 3, 'omitnan'); 
% Master_chl_a_median = median(Master_chl_a, 3, 'omitnan'); 
 
%%

%%% CHANGE DIRECTORY!!!!!

cd /Users/srishtidasarathy/Documents/Bowman/New_mat_files_CHL_A_DATA_MONTHLY

%%
save('Aqua_Latitude_Subset_BellingshausenSea.mat', 'Latitude_subset', '-v7.3') 
save('Aqua_Longitude_Subset_BellingshausenSea.mat', 'Longitude_subset', '-v7.3') 
save('Master_chlor_a_monthly_full_res.mat', 'Master_chl_a', '-v7.3') 

save('Aqua.mat',...
    'Master_time',...
    'Master_count_n', ...
    'Master_filename',...
    'Master_chl_a', '-v7.3');
%     'Master_chl_a_mean',...
%     'Master_chl_a_median', '-v7.3');
%       
  