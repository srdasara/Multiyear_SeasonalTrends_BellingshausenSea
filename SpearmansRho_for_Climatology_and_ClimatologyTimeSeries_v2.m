%%

% Here, I'm looking at the temporal correlation between CMOD, Ice, chl-a, and wind speed and plotting a figure to visualize trends. I have calculated rho with a lag, as seen in this script. 


%%


black = rgb('black'); 
grass_green = rgb('true green');
ice_blue = rgb('lightish blue'); 
ice_blue = rgb('dodger blue');
wind_blue = rgb('scarlet');

%%
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication

load('Depol_Ratio_Monthly_avg_Vars.mat')
load('CMOD_Monthly_avg_Vars_Surface.mat')
load('amsrmf_Monthly_avg_Vars.mat') 


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/New_mat_files_CHL_A_DATA_MONTHLY/

load('Master_chlor_a_monthly_full_res.mat') 
for i = 1:151 
Total_chl_monthly_2(i) = nanmean(Master_chl_a(:,:,i), [1 2]); 
end

Total_chl_a_monthly = Total_chl_monthly_2'; 
% 
%

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Revisions_Figures
%%

 Total_chl_a_monthly(Total_chl_a_monthly == 0) = NaN; 

% ^ Use this if you need to make chlorophyll 0 values into NaNs


% First make all chl_a_values that are zero for no data into NaNs...

%%
% NOTES on how Spearman Rank Correlation works...
% a = [1 4 6 3 4 6 7 8]; 
% b = [34 56 34 56 79 23 48 28];
% [RHO,PVAL] = corr(a',b','Type','Spearman');

%%
[rho_no_lag_Ice_CMOD, pval_no_lag_Ice_CMOD]   = corr(CMOD_Monthly_avg_Surface, Depol_Ratio_Monthly_avg, 'Type', 'Spearman', 'Rows', 'pairwise');
[rho_no_lag_Wind_CMOD, pval_no_lag_Wind_CMOD] = corr(CMOD_Monthly_avg_Surface, amsrmf_Monthly_avg, 'Type', 'Spearman','Rows', 'pairwise');
[rho_no_lag_Chl_a_CMOD, pval_no_lag_Chl_a_CMOD] = corr(CMOD_Monthly_avg_Surface, Total_chl_a_monthly, 'Type', 'Spearman','Rows', 'pairwise');
  
%% are ranks tied or not when performing spearman's? here's a brief check to see if values are the same. 

sorted_CMOD = sort(CMOD_Monthly_avg_Surface); 
sorted_Ice = sort(Depol_Ratio_Monthly_avg); 
sorted_Chl = sort(Total_chl_a_monthly); 

% yes, I do have tied ranks! Keep in mind the adjusted formula for
% calculating rho 


%%
for iShift=1:6
    
    disp(iShift) 
    
    CMOD_negative_shift = [CMOD_Monthly_avg_Surface(iShift + 1 : end)  ; nan(iShift, 1)];
    
          
        % This is the negative lag. CMOD Night monthly average is shifted
        % to the left by each iteration. Here I am choosing a 6 period iteration, or half a
        % year, to find where the corr coefs line up best for all my data. 
         
        
    CMOD_positive_shift = [nan(iShift, 1) ; CMOD_Monthly_avg_Surface(1 : end - iShift)];
    
        % This is the positive lag. CMOD Night monthly average is shifted
        % to the right by each iteration. Here I am choosing a 6 period iteration, or half a
        % year, to find where the corr coefs line up best for all my data. 
        
    [rho_pos_Ice_CMOD, pval_pos_Ice_CMOD] = corr(CMOD_positive_shift, Depol_Ratio_Monthly_avg, 'Type', 'Spearman', 'Rows', 'pairwise');
    [rho_neg_Ice_CMOD, pval_neg_Ice_CMOD] = corr(CMOD_negative_shift, Depol_Ratio_Monthly_avg, 'Type', 'Spearman', 'Rows', 'pairwise');
    
    rho_positive_lag_Ice_CMOD(iShift,1) = rho_pos_Ice_CMOD;
    rho_negative_lag_Ice_CMOD(iShift,1) = rho_neg_Ice_CMOD;
    
    pval_positive_lag_Ice_CMOD(iShift, 1) = pval_pos_Ice_CMOD; 
    pval_negative_lag_Ice_CMOD(iShift, 1) = pval_neg_Ice_CMOD;
    
    [rho_pos_Wind_CMOD, pval_pos_Wind_CMOD] = corr(CMOD_positive_shift, amsrmf_Monthly_avg, 'Type', 'Spearman', 'Rows', 'pairwise');
    [rho_neg_Wind_CMOD, pval_neg_Wind_CMOD] = corr(CMOD_negative_shift, amsrmf_Monthly_avg, 'Type', 'Spearman', 'Rows', 'pairwise');
    
    rho_positive_lag_Wind_CMOD(iShift,1) = rho_pos_Wind_CMOD;
    rho_negative_lag_Wind_CMOD(iShift,1) = rho_neg_Wind_CMOD;
    
    pval_positive_lag_Wind_CMOD(iShift, 1) = pval_pos_Wind_CMOD; 
    pval_negative_lag_Wind_CMOD(iShift,1)  = pval_neg_Wind_CMOD;
    
    [rho_pos_Chl_a_CMOD, pval_pos_Chl_a_CMOD] = corr(CMOD_positive_shift, Total_chl_a_monthly, 'Type', 'Spearman', 'Rows', 'pairwise');
    [rho_neg_Chl_a_CMOD, pval_neg_Chl_a_CMOD] = corr(CMOD_negative_shift, Total_chl_a_monthly, 'Type', 'Spearman', 'Rows', 'pairwise');
    
    rho_positive_lag_Chl_a_CMOD(iShift,1) = rho_pos_Chl_a_CMOD;
    rho_negative_lag_Chl_a_CMOD(iShift,1) = rho_neg_Chl_a_CMOD;
    
    pval_positive_lag_Chl_a_CMOD(iShift, 1) = pval_pos_Chl_a_CMOD;
    pval_negative_lag_Chl_a_CMOD(iShift, 1) = pval_neg_Chl_a_CMOD; 
    
end

%%

% I am flipping it because in the above loop, the iteration goes
% sequentially from 1 to 6. HOWEVER, when I want to look at all of my
% correlation coefficients from - 6 to 6, I need to flip the negative lag
% vector because it goes from -1:-6 versus the other way around, necessary
% for concatenating... 

rho_negative_lag_Chl_a_CMOD_flipped = flipud(rho_negative_lag_Chl_a_CMOD) ; 
rho_negative_lag_Wind_CMOD_flipped  = flipud(rho_negative_lag_Wind_CMOD); 
rho_negative_lag_Ice_CMOD_flipped   = flipud(rho_negative_lag_Ice_CMOD);

pval_negative_lag_Chl_a_CMOD_flipped = flipud(pval_negative_lag_Chl_a_CMOD);
pval_negative_lag_Wind_CMOD_flipped  = flipud(pval_negative_lag_Wind_CMOD);
pval_negative_lag_Ice_CMOD_flipped   = flipud(pval_negative_lag_Ice_CMOD); 


%%
spearman_rho_Chl_a_CMOD = [rho_negative_lag_Chl_a_CMOD_flipped ; rho_no_lag_Chl_a_CMOD ; rho_positive_lag_Chl_a_CMOD]; 

spearman_rho_Wind_CMOD  = [rho_negative_lag_Wind_CMOD_flipped ; rho_no_lag_Wind_CMOD ; rho_positive_lag_Wind_CMOD]; 

spearman_rho_Ice_CMOD   = [rho_negative_lag_Ice_CMOD_flipped ; rho_no_lag_Ice_CMOD ; rho_positive_lag_Ice_CMOD ] ;


pval_Chl_a_CMOD = [pval_negative_lag_Chl_a_CMOD_flipped ; pval_no_lag_Chl_a_CMOD ; pval_positive_lag_Chl_a_CMOD];

pval_Wind_CMOD = [pval_negative_lag_Wind_CMOD_flipped ; pval_no_lag_Wind_CMOD ; pval_positive_lag_Wind_CMOD ];

pval_Ice_CMOD = [pval_negative_lag_Ice_CMOD_flipped ; pval_no_lag_Ice_CMOD ; pval_positive_lag_Ice_CMOD ] ; 




%%
%
fig = figure(1); clf;  

% subplot(3,1,1);
plot(-6:6, spearman_rho_Wind_CMOD,...
    'LineWidth',2,...
    'Color', wind_blue) 
grid on

AX=findobj(fig,'Type','Axes'); 
set(AX, 'FontSize', 16)    

set(AX,'YAxisLocation','origin')
set(AX, 'Ylim', [-1 1]); %ax.YLim = [-2 2];
set(AX, 'YTick', -1:0.2:1 )       

set(AX, 'XAxisLocation', 'origin')
set(AX, 'XTick', -6:1:6)

[ax,h1] = suplabel('Spearman''s rho', 'y');  %#ok<ASGLU>
[ax,h2] = suplabel('Lag Period, Monthly','x');   %#ok<ASGLU>

set(h1,'FontSize', 18)
set(h2, 'FontSize', 18)

title('CMOD and Wind')

%%
saveas(gcf, 'Spearmans_rho_CMOD_Wind', 'fig') ; 
saveas(gcf, 'Spearmans_rho_CMOD_Wind', 'png') ; 

%%
fig = figure(2); clf;

plot(-6:6, spearman_rho_Ice_CMOD, ... 
    'LineWidth', 2,...
    'Color', ice_blue)

grid on

AX=findobj(fig,'Type','Axes'); 
set(AX, 'FontSize', 16)    

set(AX,'YAxisLocation','origin')
set(AX, 'Ylim', [-1 1]); 
set(AX, 'YTick', -1 : 0.2 : 1 )       

set(AX, 'XAxisLocation', 'origin')
set(AX, 'XTick', -6:1:6)

[ax,h1] = suplabel('Spearman''s rho', 'y');  %#ok<ASGLU>
[ax,h2] = suplabel('Lag Period, Monthly','x');   %#ok<ASGLU>

set(h1,'FontSize', 18)
set(h2, 'FontSize', 18)

title('CMOD and Ice') 

%%
saveas(gcf, 'Spearmans_rho_CMOD_Ice', 'fig') ; 
saveas(gcf, 'Spearmans_rho_CMOD_Ice', 'png') ; 

%%
% subplot(3,1,3);

fig = figure(3); clf; 

plot(-6:6, spearman_rho_Chl_a_CMOD, ...
    'LineWidth', 2,...
    'Color', grass_green)
grid on


AX=findobj(fig,'Type','Axes'); 
set(AX, 'FontSize', 16)    

set(AX,'YAxisLocation','origin')
set(AX, 'Ylim', [-1 1]); %ax.YLim = [-2 2];
set(AX, 'YTick', -1:0.2:1 )       

set(AX, 'XAxisLocation', 'origin')
set(AX, 'XTick', -6:1:6)

[ax,h1] = suplabel('Spearman''s rho', 'y');  %#ok<ASGLU>
[ax,h2] = suplabel('Lag Period, Monthly','x');   %#ok<ASGLU>

set(h1,'FontSize', 18)
set(h2, 'FontSize', 18)

%%


saveas(gcf, 'Spearmans_rho_CMOD_chl_a', 'fig') ; 
saveas(gcf, 'Spearmans_rho_CMOD_chl_a', 'png') ; 








% %%
% saveas(gcf, 'spearmans_rho_6mthlag_CMOD_Wind_Ice_Chl_a', 'fig') ; 
% 
% saveas(gcf, 'spearmans_rho_6mthlag_CMOD_Wind_Ice_Chl_a.png') ; 

%%




%%
%
fig = figure; clf;  

subplot(3,1,1);
plot(-6:6, spearman_rho_Wind_CMOD,'-kd',...
    'LineWidth',2,...
    'MarkerSize', 8,...
    'MarkerFaceColor', wind_blue) 
grid on

title('CMOD and Wind')

subplot(3,1,2); 
plot(-6:6, spearman_rho_Ice_CMOD, '-kd', ... 
    'LineWidth', 2,...
    'MarkerSize', 8,...
    'MarkerFaceColor', ice_blue)
grid on

title('CMOD and Ice') 

subplot(3,1,3);
plot(-6:6, spearman_rho_Chl_a_CMOD, '-kd', ...
    'LineWidth', 2,...
    'MarkerSize', 8,...
    'MarkerFaceColor', grass_green)
grid on

title('CMOD and Chlorophyll') 

AX=findobj(fig,'Type','Axes'); 
set(AX, 'FontSize', 13)    
set(AX,'YAxisLocation','origin')
set(AX, 'XAxisLocation', 'origin')
set(AX, 'XTick', -6:1:6)
set(AX, 'YTick', -1: 0.1 : 1 )

[ax,h1] = suplabel('Spearman''s rho', 'y');  %#ok<ASGLU>
[ax,h2] = suplabel('Lag Period','x');  

set(h1, 'FontWeight', 'bold','FontSize', 20)
set(h2, 'FontWeight', 'bold', 'FontSize', 20)


%%
saveas(gcf, 'spearmans_rho_6mthlag_CMOD_Wind_Ice_Chl_a', 'fig') ; 

saveas(gcf, 'spearmans_rho_6mthlag_CMOD_Wind_Ice_Chl_a.png') ; 



%%

make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.5], [0.125 0.125], [0.125 0.125]);
 subplot = @(m,n,p) subtightplot (m, n, p, [0.035 0.05], [0.1 0.05], [0.1 0.05]);

if ~make_it_tight,  clear subplot;  end


clear plot
fig = figure(4);clf; 



ax(1) = subplot(3,1,1);

plot(-6:6, spearman_rho_Chl_a_CMOD, '-v', ...
    'LineWidth', 2,...
    'MarkerSize', 6,...
    'MarkerFaceColor', grass_green,...
    'Color', grass_green)
ylim([-0.8 0.8])
set(gca, 'xtick', -6:1:6,'ytick',  -0.6:0.2:0.6)
legend('MAOD & chl-{\ita}')
legend boxoff
grid on 

% title('CMOD and Wind')

ax(2) = subplot(3,1,2); 
plot(-6:6, spearman_rho_Ice_CMOD, '-o', ... 
    'LineWidth', 2,...
    'MarkerSize', 6,...
    'MarkerFaceColor', ice_blue,...
    'Color', ice_blue)
ylim([-0.8 0.8])
% ax(2).XTickLabel = [];

set(gca, 'xtick', -6:1:6,'ytick',  -0.6:0.2:0.6)
legend('MAOD & Ice')
legend boxoff

grid on

% title('CMOD and Ice') 

ax(3) = subplot(3,1,3);

plot(-6:6, spearman_rho_Wind_CMOD,'-d',...
    'LineWidth',2,...    
    'MarkerSize', 6,...
    'MarkerFaceColor', wind_blue,...
    'Color', wind_blue)
ylim([-0.8 0.8])
% ax(1).XTickLabel = [];
set(gca, 'xtick', -6:1:6, 'ytick', -0.6:0.2:0.6)
legend('MAOD & Wind Speed')
legend boxoff

ylim([-0.8 0.8])

grid on

[ax(1:2).XTickLabel] = deal([]);

% 
 [s,h1] = suplabel('rho', 'y');  %#ok<ASGLU>
 [s,h2] = suplabel('Lag Period, Monthly','x');  
% 

AX=findall(0,'type','axes');
set(AX, 'FontSize', 19)
% set(AX, 'XTick', -6:1:6)

%     [0.617957746478874 0.799442896935933 0.248239436619719 0.0668523676880223],...
% 
% an(1) = annotation(fig,'textbox',...
%     [0.46 0.73 0.248239436619719 0.0668523676880223],...
%     'String','CMOD & Chl-{\ita} Temporal Correlation',...
%     'LineStyle','none',...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FontName','Helvetica-Narrow',...
%     'FitBoxToText','off');
% 
% an(2) = annotation(fig,'textbox',...
%     [an(1).Position(1),...
%     (ax(2).Position(2) + (an(1).Position(2) - ax(1).Position(2))),... % so this is the second axes, y position  + (y position of text box - y position of first axes)
%     an(1).Position(3) an(1).Position(4)],...
%     'String',{'CMOD & Ice Temporal Correlation'},...
%     'LineStyle','none',...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FontName','Helvetica-Narrow',...
%     'FitBoxToText','off');
% 
% an(3) = annotation(fig, 'textbox', ...
%     [an(1).Position(1),...
%     (ax(3).Position(2) + (an(2).Position(2) - (ax(2).Position(2)))),...
%     an(1).Position(3) an(1).Position(4)],...
%     'String',{'CMOD & Wind Temporal Correlation'},...
%     'LineStyle','none',...
%     'FontWeight','bold',...
%     'FontSize',14,...
%     'FontName','Helvetica-Narrow',...
%     'FitBoxToText','off');
% 
% 
% 

 c = newline;

an(4) = annotation(fig,'textbox',...
    [0.77 0.87 0.4 0.0376044568245124],...
    'String','lag +1, rho = 0.6181',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica-Narrow',...
    'FitBoxToText','off');



an(5) = annotation(fig,'textbox',...
    [0.76,...
    (ax(2).Position(2) + (an(4).Position(2) - ax(1).Position(2))),... % so this is the second axes, y position  + (y position of text box - y position of first axes)
    an(4).Position(3) an(4).Position(4)],...
    'String','lag +3, rho = -0.5873',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica-Narrow',...
    'FitBoxToText','off');

an(6) = annotation(fig, 'textbox', ...
    [0.76,...
    (ax(3).Position(2) + (an(5).Position(2) - (ax(2).Position(2)))),...
    an(4).Position(3) an(4).Position(4)],...
    'String','lag +1, rho = -0.5986',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica-Narrow',...
    'FitBoxToText','off');




%%
% orient(fig,'landscape')

 set(gcf,'PaperPositionMode','auto')


print(fig,'Updated_SpearmansRho_Fig.png','-dpng','-r300');       %  *// 300 dpi
%  print(gcf, 'Ju













