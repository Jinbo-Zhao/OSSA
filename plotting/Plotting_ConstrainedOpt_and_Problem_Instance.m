clc
clear
close all

folders = strcat({'..\'}, {'data','src'});
for k = 1:length(folders)
    addpath(folders{k})
end
%% Loading data
% loading the experiments data
load('Exp_Env_Data_Constrained.mat')
load('Exp_data_BootCamp_median_beta_0_05.mat')
%load('Plotting_Data.mat')
exp_result_Boot=returnedtruesettingbeta005;
frquency_Boot=sum(exp_result_Boot);
return_amount_Boot=sum(exp_result_Boot,2);
true_mean(:,1)=1+true_mean(:,1);%Change the negative utility rate to idle rate

fesible_systems=true_mean(true_mean(:,2)<feasible,:);
[optimal_primary,optimal_index]=min(fesible_systems(:,1));
optimal_system=fesible_systems(optimal_index,:);
%% Gnerating Used Plots
%Two figures were generated om thi sub block, which is the Fig6 in the
%first version manuscipt
figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_half,'filled')
hold on
p1=scatter(optimal_system(:,1),optimal_system(:,2),[],'r',LineWidth=1.5);
hold on
p2=yline(feasible,'r',LineWidth=1.5);
xlim([0.26 0.32])
ylim([6 9])
xlabel('Expected Average Idle Rate')
ylabel('Expected Average Waiting Time')
cb=colorbar();
ylabel(cb,'Frequency of Being Returned','FontSize',15,'Rotation',270,Position=[3.760952495393299,490.7676050786943,0])
cb.Ticks = 0:100:1000;
cb.TickLabels = num2cell(0:0.1:1) ;
h1=legend([p1 p2],{'Acceptable System','Constraint'});
set( h1, 'Box', 'on','Position',[0.148809449995557,0.185317463084818,0.367857134608286,0.120238091974032]) ;
set(gca,"LineWidth",1.5,"FontSize",15)
f = gcf;
savefig('Constrained_Opt\Frequency_Half_Constrian_Local.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Half_Constrian_Local.png','Resolution',600)

figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_Boot,'filled')
hold on
scatter(optimal_system(:,1),optimal_system(:,2),[],'r',LineWidth=1.5)
hold on
yline(feasible,'r',LineWidth=1.5)
xlabel('Expected Average Idle Rate')
ylabel('Expected Average Waiting Time')
cb=colorbar();
caxis([0,1000])
cb.Ticks = 0:100:1000;
cb.TickLabels = num2cell(0:0.1:1) ;
ylabel(cb,'Frequency of Being Returned','FontSize',15,'Rotation',270,Position=[3.760952495393299,490.7676050786943,0])
xlim([0.26 0.32])
ylim([6 9])
set(gca,"LineWidth",1.5,"FontSize",15)
f = gcf;
savefig('Constrained_Opt\Frequency_Boot_Constrian_Local.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Boot_Constrian_Local.png','Resolution',600)

%% Gnerating full plots
%In this subsection, we generated 4 colored scatter plots corresponding to
%3 different confidence regions and BOOTCOMP
optimal_system=fesible_systems(optimal_index,:);

true_mean=flip(true_mean);
frquency_ellipse=flip(frquency_ellipse);
frquency_box=flip(frquency_box);
frquency_half=flip(frquency_half);
frquency_Boot=flip(frquency_Boot);
figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_ellipse,'filled')
%title('Returned Frequency Ellipse')
hold on
p1=scatter(optimal_system(:,1),optimal_system(:,2),[],'r');
hold on
p2=yline(feasible,'r');
xlabel('Expected Average Idle Rate')
ylabel('Expected Average Waiting Time')
cb=colorbar();
ylabel(cb,'Propotion Being Returned','FontSize',11,'Rotation',270,Position=[3.38,503,0])
cb.Ticks = 0:100:1000;
cb.TickLabels = num2cell(0:0.1:1) ;
h1=legend([p1 p2],{'Acceptable System','Constraint'});
set( h1, 'Box', 'on', 'Color', [0.8,0.8,0.8], 'EdgeColor', get( h1, 'Color' )) ;
f = gcf;
savefig('Constrained_Opt\Frequency_Ellipe_full.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Ellipe_full.png','Resolution',600)

figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_box,'filled')
colorbar
%title('Returned Frequency Box')
hold on
scatter(optimal_system(:,1),optimal_system(:,2),[],'r')
hold on
yline(feasible,'r')
xlabel('Expected Average Idle Rate')
ylabel('Expected Average Waiting Time')
cb=colorbar();
ylabel(cb,'Propotion Being Returned','FontSize',11,'Rotation',270,Position=[3.38,503,0])
cb.Ticks = 0:100:1000;
cb.TickLabels = num2cell(0:0.1:1) ;
f = gcf;
savefig('Constrained_Opt\Frequency_Box_full.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Box_full.png','Resolution',600)

figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_half,'filled')%Scatter
%title('Returned Frequency Ellipse')
hold on
p1=scatter(optimal_system(:,1),optimal_system(:,2),[],'r');%Optimal system
hold on
p2=yline(feasible,'r');%Constraint
xlabel('Expected Average Idle Rate')%Labels
ylabel('Expected Average Waiting Time')
cb=colorbar();
ylabel(cb,'Propotion Being Returned','FontSize',11,'Rotation',270,Position=[3.38,503,0])
cb.Ticks = 0:100:1000;%Colorbar ticks
cb.TickLabels = num2cell(0:0.1:1) ;
h1=legend([p1 p2],{'Acceptable System','Constraint'});%Legend
set( h1, 'Box', 'on', 'Color', [0.8,0.8,0.8], 'EdgeColor', get( h1, 'Color' )) ;
rectangle('Position',[0.26 6 0.06 3])
f = gcf;%Saving graph
savefig('Constrained_Opt\Frequency_Half_full.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Half_full.png','Resolution',600)

figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_Boot,'filled')
%title('Returned Frequency Ellipse')
hold on
p1=scatter(optimal_system(:,1),optimal_system(:,2),[],'r');
hold on
p2=yline(feasible,'r');
xlabel('Expected Average Idle Rate')
ylabel('Expected Average Waiting Time')
cb=colorbar();
ylabel(cb,'Propotion Being Returned','FontSize',11,'Rotation',270,Position=[3.38,503,0])
caxis([0,1000])
cb.Ticks = [0:100:1000];
cb.TickLabels = num2cell(0:0.1:1) ;
rectangle('Position',[0.26 6 0.06 3])
f = gcf;
savefig('Constrained_Opt\Frequency_Boot_full.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Boot_full.png','Resolution',600)

%% Gnerating zoom-in plots
%In this subsection, we generated 4 local behavior of 
% colored scatter plots corresponding to 3 different confidence regions and
% BOOTCOMP. These were for zoom-in plots, whcih is not used.
figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_ellipse,'filled')
%title('Returned Frequency Ellipse')
hold on
p1=scatter(optimal_system(:,1),optimal_system(:,2),[],'r');
hold on
p2=yline(feasible,'r');
xlim([0.26 0.32])
ylim([6 9])
f = gcf;
savefig('Constrained_Opt\Frequency_Ellipsoid_local.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Ellipsoid_local.png','Resolution',600)

figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_box,'filled')
%title('Returned Frequency Box')
hold on
scatter(optimal_system(:,1),optimal_system(:,2),[],'r')
hold on
yline(feasible,'r')
xlim([0.26 0.32])
ylim([6 9])
f = gcf;
savefig('Constrained_Opt\Frequency_Box_local.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Box_local.png','Resolution',600)

figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_half,'filled')
%title('Returned Frequency Half')
hold on
scatter(optimal_system(:,1),optimal_system(:,2),[],'r')
hold on
yline(feasible,'r')
xlim([0.26 0.32])
ylim([6 9])
f = gcf;
savefig('Constrained_Opt\Frequency_Half_local.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Half_local.png','Resolution',600)

figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_Boot,'filled')
caxis([0 1000])
%title('Returned Frequency Boot')
hold on
scatter(optimal_system(:,1),optimal_system(:,2),[],'r')
hold on
yline(feasible,'r')
xlim([0.26 0.32])
ylim([6 9])
f = gcf;
savefig('Constrained_Opt\Frequency_Boot_local.fig')
exportgraphics(f,'Constrained_Opt\Frequency_Boot_local.png','Resolution',600)
%% Reporting Type I error and Generating Returned Amount Histo
disp('Type I error of FOSSA Ellipsoid is')
frquency_ellipse(find(true_mean(:,1) ==optimal_system(:,1)))
disp('Type I error of FOSSA Box is')
frquency_box(find(true_mean(:,1) ==optimal_system(:,1)))
disp('Type I error of FOSSA HalfBox is')
frquency_half(find(true_mean(:,1) ==optimal_system(:,1)))
disp('Type I error of BootComp is')
frquency_Boot(find(true_mean(:,1) ==optimal_system(:,1)))

BW=1;
figure
histogram(return_amount_ellipse,BinWidth=BW)
hold on
histogram(return_amount_box,BinWidth=BW)
hold on
histogram(return_amount_half,BinWidth=BW)
hold on
histogram(return_amount_Boot,BinWidth=BW)
h1=legend({'Ellipsoid','Box','Half-box','BootComp'});
xlabel('Subset Size')
ylabel('Frequency')
set( h1,'Box', 'off','Position',[0.434523796511722,0.711507941049243,0.192857139823692,0.165476185934884]);
f = gcf;
savefig('Constrained_Opt\SubsetSize_Constrained.fig')
exportgraphics(f,'Constrained_Opt\SubsetSize_Constrained.png','Resolution',600)
%title('Size of Returned Subset')
%% Standardize the Performance
true_mean=flip(true_mean);
utility=true_mean(:,1);
waiting_time=true_mean(:,2);
utility_min=min(utility);
utility_max=max(utility);
utility_range=utility_max-utility_min;
waiting_time_min=min(waiting_time);
waiting_time_max=max(waiting_time);
waiting_time_range=waiting_time_max-waiting_time_min;

standardized_mean(:,1)=(true_mean(:,1)-utility_min)/utility_range;
standardized_mean(:,2)=(true_mean(:,2)-waiting_time_min)/waiting_time_range;
standardized_optimal_primary=(optimal_primary-utility_min)/utility_range;
standardized_feasible=(feasible-waiting_time_min)/waiting_time_range;

standard_point=[standardized_optimal_primary,standardized_feasible];

quality_all= quality_Constrained(standardized_mean,standard_point);
%quality_returned= quality_all(logical(return_index_ellipse),:);
%%
quality_storage=zeros(4,macro_replications);
for rep=1:macro_replications
    temp=exp_result_ellipse(rep,:);
    quality_returned=quality_all(logical(temp),:);
    quality_storage(1,rep)=mean(quality_returned);
    
    temp=exp_result_box(rep,:);
    quality_returned=quality_all(logical(temp),:);
    quality_storage(2,rep)=mean(quality_returned);

    temp=exp_result_half(rep,:);
    quality_returned=quality_all(logical(temp),:);
    quality_storage(3,rep)=mean(quality_returned);

    temp=exp_result_Boot(rep,:);
    quality_returned=quality_all(logical(temp),:);
    quality_storage(4,rep)=mean(quality_returned);
end


figure
average_quality_ellipse=quality_storage(1,:);
histogram(average_quality_ellipse,BinWidth=0.001)
hold on
average_quality_box=quality_storage(2,:);
histogram(average_quality_box,BinWidth=0.001)
hold on
average_quality_half=quality_storage(3,:);
histogram(average_quality_half,BinWidth=0.001)
hold on
average_quality_Boot=quality_storage(4,:);
histogram(average_quality_Boot,BinWidth=0.001)

h1=legend({'Ellipsoid','Box','Half-box','BootComp'});
xlabel('Average Quality')
ylabel('Frequency')
set( h1,'Box', 'off','Position',[0.161309510797436,0.704365083906386,0.192857139823692,0.165476185934884]);
%title('Avergae Quality of Returned Subset')
f = gcf;
savefig('Constrained_Opt\Average Quality_Constrained.fig')
exportgraphics(f,'Constrained_Opt\Average Quality_Constrained.png','Resolution',600)

%%
figure
p1=histogram(quality_all,BinWidth=0.05,FaceColor=[200 200 200]/255);
set(gca,'Box','off')
hold on
temp=exp_result_half(rep,:);
quality_returned=quality_all(logical(temp),:);
p2=histogram(quality_returned,BinWidth=0.05,FaceColor=	"#EDB120",FaceAlpha=1);
set(gca,'Box','off')
h1=legend([p2,p1],{'$\mathcal{S}^{\mathrm{FOSSA}}$','Other systems'},'Interpreter','latex');
set( h1,'Box', 'off','Position',[0.06202381798554,0.757960320383784,0.374999991538269,0.137999996423722]);
xlabel('Distance to Acceptability')
ylabel('Number of Systems')
f = gcf;
set(gcf, 'Position', [100 100 600 250]); 
savefig('Constrained_Opt\Quality of Returned Systems_Constrained.fig')
exportgraphics(f,'Constrained_Opt\Quality of Returned Systems_Constrained.png','Resolution',600)
%% 
load('Exp_data_time_constrained.mat')

disp('The 10%quantile mean and 90%quantile of running time of FOSSA Elipsoid are')
quantile(time_result_ellipse,0.1)
mean(time_result_ellipse)
quantile(time_result_ellipse,0.9)
disp('The 10%quantile mean and 90%quantile of running time of FOSSA Box are')
quantile(time_result_box,0.1)
mean(time_result_box)
quantile(time_result_box,0.9)
disp('The 10%quantile mean and 90%quantile of running time of FOSSA HalfBox are')
quantile(time_result_half,0.1)
mean(time_result_half)
quantile(time_result_half,0.9)
disp('The 10%quantile mean and 90%quantile of running time of BootComp are')
quantile(time_result_Boot,0.1)
mean(time_result_Boot)
quantile(time_result_Boot,0.9)

%%
disp('The 10%quantile mean and 90%quantile of average quanlity of returned subset(FOSSA Elipsoid) are')
quantile(average_quality_ellipse,0.1)
mean(average_quality_ellipse)
quantile(average_quality_ellipse,0.9)
disp('The 10%quantile mean and 90%quantile of average quanlity of returned subset(FOSSA Box) are')
quantile(average_quality_box,0.1)
mean(average_quality_box)
quantile(average_quality_box,0.9)
disp('The 10%quantile mean and 90%quantile of average quanlity of returned subset(FOSSA HalfBox) are')
quantile(average_quality_half,0.1)
mean(average_quality_half)
quantile(average_quality_half,0.9)
disp('The 10%quantile mean and 90%quantile of average quanlity of returned subset(BootComp) are')
quantile(average_quality_Boot,0.1)
mean(average_quality_Boot)
quantile(average_quality_Boot,0.9)
%%
disp('The 10%quantile mean and 90%quantile of size of returned subset(FOSSA Elipsoid) are')
quantile(return_amount_ellipse,0.1)
mean(return_amount_ellipse)
quantile(return_amount_ellipse,0.9)
disp('The 10%quantile mean and 90%quantile of size of returned subset(FOSSA Box) are')
quantile(return_amount_box,0.1)
mean(return_amount_box)
quantile(return_amount_box,0.9)
disp('The 10%quantile mean and 90%quantile of size of returned subset(FOSSA HalfBox) are')
quantile(return_amount_half,0.1)
mean(return_amount_half)
quantile(return_amount_half,0.9)
disp('The 10%quantile mean and 90%quantile of size of returned subset(FOSSA BootComp) are')
quantile(return_amount_Boot,0.1)
mean(return_amount_Boot)
quantile(return_amount_Boot,0.9)
%% Scatter instance
figure
scatter(true_mean(:,1),true_mean(:,2),'black','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
set(gca,'linewidth',1.5,'fontsize',15)
%title('Returned Frequency Ellipse')
%hold on
%p2=yline(5,'r');
xlabel('Expected Average Idle Rate')
ylabel('Expected Average Waiting Time')
f = gcf;
savefig('Problem Instance\Scatter_All_Systems.fig')
exportgraphics(f,'Problem Instance\Scatter_All_Systems.png','Resolution',600)
%% Scatter of samples form one specific system
load('Simulation_Data.mat')

sample_size=50000;
transparency=0.008;
system_index=1000;%Change here to see different systems
temp=record(system_index,:,1:sample_size);
temp=reshape(temp,2,sample_size);
temp=temp';
temp(:,1)=1-temp(:,1);

figure
scatter(temp(1:1000,1),temp(1:1000,2),'black')
xlabel('Average Idle Rate')
ylabel('Average Waiting Time')
f = gcf;
savefig('Problem Instance\Scatter_one_System_Dots_1000.fig')
exportgraphics(f,'Problem Instance\Scatter_one_System_Dots_1000.png','Resolution',600)

figure
scatter(temp(:,1),temp(:,2),'black','filled','MarkerFaceAlpha',transparency,'MarkerEdgeAlpha',transparency)
set(gca,'linewidth',1.5,'fontsize',15)
xlabel('Average Idle Rate')
ylabel('Average Waiting Time')
f = gcf;
savefig('Problem Instance\Scatter_one_System_cloud_10000.fig')
exportgraphics(f,'Problem Instance\Scatter_one_System_cloud_10000.png','Resolution',600)
% This a figure of probability cloud of samples from the normal distn
% taking the sample mean of sample covaraince of slected system as the true parameter
%It is printed for showing the normality of slected system
NV = mvnrnd(mean(temp),cov(temp),sample_size);
%NV = mvnrnd(mean(temp),cov(temp),10000);
figure
scatter(NV(:,1),NV(:,2),'black','filled','MarkerFaceAlpha',transparency,'MarkerEdgeAlpha',transparency)
ylim=[22,34];
xlabel('Average Idle Rate')
ylabel('Average Waiting Time')
f = gcf;
savefig('Problem Instance\Scatter_Normal_cloud_10000.fig')
exportgraphics(f,'Problem Instance\Scatter_Normal_cloud_10000.png','Resolution',600)

%% The histogram of correlation coefficent of all systems
load('Simulation_Data.mat')
record(:,1,:)=1-record(:,1,:);
corelation_stoarage=zeros(1001,1);
for rep =1:1001
    temp=record(rep,:,:);
    temp=reshape(temp,2,50000);
    temp=cov(temp');
    corelation_stoarage(rep,1)=temp(2,1)/sqrt(temp(1,1)*temp(2,2));
end
figure
histogram(corelation_stoarage)
xlabel('Correlation Coefficient')
ylabel('Frequency')
f = gcf;
savefig('Problem Instance\Histo_Correlation.fig')
exportgraphics(f,'Problem Instance\Histo_Correlation.png','Resolution',600)

%%
for k = 1:length(folders)
    rmpath(folders{k})
end