clc
clear
close all

folders = strcat({'..\'}, {'data','src'});
for k = 1:length(folders)
    addpath(folders{k})
end
%%
load('Exp_Env_Data_Pareto_full.mat')
load('Exp_Result_GSPRT.mat')

frquency_ellipse=sum(exp_result_ellipse);
frquency_box=sum(exp_result_box);
frquency_half=sum(exp_result_half);
frquency_GSPRT=sum(exp_result_GSPRT);
%true_mean(:,1)=1+true_mean(:,1);
[membership, member_value]=find_pareto_frontier(true_mean);
optimal_system=member_value;

%% Confidence Region For Pareto Front
Plotting_confidence_box_updated(system_info,sample_size_vec,0.05,'#FFFF00');
hold on
%scatter(true_mean(:,1),true_mean(:,2),20,'black','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%hold on
scatter(true_mean(:,1),true_mean(:,2),10,'black','filled')
hold on
scatter(optimal_system(:,1),optimal_system(:,2),10,'filled',MarkerFaceColor='#D92525',MarkerEdgeColor='none')
hold on
histogram(-2+rand(10,1),'FaceColor','#FFFF00',EdgeColor='none');%Stealing its legend
h1=legend({'All systems','Acceptable Systems','Confidence Region'});
%set( h1, 'Box', 'on', 'Color', [0.8,0.8,0.8], 'EdgeColor', get( h1, 'Color' )) ;
set( h1, 'Box', 'on') ;
%set( h1,'Box', 'off')
xlabel('Expected Average Idle Rate')
ylabel('Expected Average Waiting Time')
xlim([0.2,0.55])
ylim([0,40])

f = gcf;
savefig('Pareto Opt\Confidence_Region_Box_Illustration_new.fig')
exportgraphics(f,'Pareto Opt\Confidence_Region_Box_Illustration_new.png','Resolution',600)

%% Checking PAS
count=0;
for rep=1:1000
    
    temp=exp_result_ellipse(rep,:);
    temp=sum(temp(membership));
    if temp~=sum(membership)
        count=count+1;
        temp
    end
    
end
disp('Ellipsoid Method missing Optimal system:')
count

count=0;
for rep=1:1000
    
    temp=exp_result_box(rep,:);
    temp=sum(temp(membership));
    if temp~=sum(membership)
        count=count+1;
        temp
    end
    
end
disp('Box Method missing Optimal system:')
count

count=0;
for rep=1:1000
    
    temp=exp_result_half(rep,:);
    temp=sum(temp(membership));
    if temp~=sum(membership)
        count=count+1;
        temp
    end
    
end
disp('Half-box Method missing Optimal system:')
count

count=0;
for rep=1:1000
    
    temp=exp_result_GSPRT(rep,:);
    temp=sum(temp(membership));
    if temp~=sum(membership)
        count=count+1;
        temp
    end
    
end
disp('GSPRT Method missing Optimal system:')
count

%% Checking Prescreening
count=0;
for rep=1:1000
    
    temp=exp_result_ellipse(rep,:);
    logic_GSPRT=exp_result_GSPRT(rep,:);
    temp=sum(temp(logic_GSPRT));
    if temp~=sum(logic_GSPRT)
        count=count+1;
        temp
    end
    
end
disp('Ellipsoid Method missing PreScreening:')
count

count=0;
for rep=1:1000
    
    temp=exp_result_box(rep,:);
    logic_GSPRT=exp_result_GSPRT(rep,:);
    temp=sum(temp(logic_GSPRT));
    if temp~=sum(logic_GSPRT)
        count=count+1;
        temp
    end
    
end
disp('Box Method missing PreScreening:')
count

count=0;
for rep=1:1000
    
    temp=exp_result_half(rep,:);
    logic_GSPRT=exp_result_GSPRT(rep,:);
    temp=sum(temp(logic_GSPRT));
    if temp~=sum(logic_GSPRT)
        count=count+1;
        temp
    end
    
end
disp('Half Method missing PreScreening:')
count
%% full plots
% The figure is the colored scatter plot for FOSSA (BOX) 
%figure
Plotting_confidence_box_updated(system_info,sample_size_vec,0.05,[140 140 140]/256);
hold on
scatter(true_mean(:,1),true_mean(:,2),[],frquency_ellipse,'filled')

%title('Returned Frequency Ellipse')
hold on
p2=scatter(optimal_system(:,1),optimal_system(:,2),[],'r');
hold on
p3=histogram(-2+rand(10,1),'FaceColor',[140 140 140]/256,EdgeColor='none');%Stealing its legend
%h1=legend([p1],{'Acceptable Systems','Confidence Region'});
h1=legend([p2,p3],{'Acceptable Systems','Confidence Region'});
xlabel('Expected Average Idle Rate','FontSize',15)
ylabel('Expected Average Waiting Time','FontSize',15)
xlim([0.2,0.6])
ylim([0,40])
cb=colorbar();
ylabel(cb,'Frequency of Being Returned','FontSize',15,'Rotation',270,Position=[4.332381066821871,499.9419170052998,0])
cb.Ticks = [0:100:1000];
cb.TickLabels = num2cell(0:0.1:1) ;
%set( h1, 'Box', 'on', 'Color', [0.8,0.8,0.8], 'EdgeColor', get( h1, 'Color' )) ;
set( h1, 'Box', 'on','Position',[0.395238127625313,0.701984130248191,0.383928562700748,0.120238091974031]) ;
set(gca,'LineWidth',1.5,'FontSize',15);
f = gcf;
savefig('Pareto Opt\Frequency_Box_Pareto.fig')
exportgraphics(f,'Pareto Opt\Frequency_Box_Pareto.png','Resolution',600)

%% The figure is the colored scatter plot for FOSSA (GSPRT)
% hold on
% histogram(-2+rand(10,1),'FaceColor','#FFFF00',EdgeColor='none');%Stealing its legend
% h1=legend({'All systems','Acceptable Systems','Confidence Region'});
% %set( h1, 'Box', 'on', 'Color', [0.8,0.8,0.8], 'EdgeColor', get( h1, 'Color' )) ;
% set( h1, 'Box', 'on') ;
% %set( h1,'Box', 'off')
% xlabel('Expected Average Idle Rate')
% ylabel('Expected Average Waiting Time')
% xlim([0.2,0.55])
% ylim([0,40])

figure
scatter(true_mean(:,1),true_mean(:,2),[],frquency_GSPRT,'filled')
%title('Returned Frequency Ellipse')
hold on
p1=scatter(optimal_system(:,1),optimal_system(:,2),[],'r');
xlabel('Expected Average Idle Rate')
ylabel('Expected Average Waiting Time')
cb=colorbar();
ylabel(cb,'Frequency of Being Returned','FontSize',15,'Rotation',270,Position=[4.332381066821871,499.9419170052998,0])
cb.Ticks = [0:100:1000];
cb.TickLabels = num2cell(0:0.1:1) ;
%h1=legend([p1],{'Acceptable Systems'});
%set( h1, 'Box', 'on', 'Color', [0.8,0.8,0.8], 'EdgeColor', get( h1, 'Color' )) ;
set(gca,'LineWidth',1.5,'FontSize',15);
f = gcf;
savefig('Pareto Opt\Frequency_GSPRT_Pareto.fig')
exportgraphics(f,'Pareto Opt\Frequency_GSPRT_Pareto.png','Resolution',600)

%% Quality figure
% Standardized the true mean
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

[~, member_value]=find_pareto_frontier(standardized_mean);
quality_all = quality_Pareto(standardized_mean,member_value);

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

    temp=exp_result_GSPRT(rep,:);
    quality_returned=quality_all(logical(temp),:);
    quality_storage(4,rep)=mean(quality_returned);
end



% Histogram of averaged quality
figure
Ellipse=quality_storage(1,:);
histogram(Ellipse,BinWidth=0.0001)
hold on
Box=quality_storage(2,:);
histogram(Box,BinWidth=0.0001)
hold on
Half=quality_storage(3,:);
histogram(Half,BinWidth=0.0001)
hold on
GSPRT=quality_storage(4,:);
histogram(GSPRT,BinWidth=0.0001)
title('Averaged Quality of Returned Subset')

%% Quality Plot of oen random returned subset
temp=exp_result_box(1000,:);
quality_returned_box=quality_all(logical(temp),:);
temp=exp_result_GSPRT(1000,:);
quality_returned_GSPRT=quality_all(logical(temp),:);

BW=0.01;
figure
histogram(quality_all,BinWidth=BW,FaceColor=[96 96 96]/255)
set(gca, 'YScale', 'log')
hold on
histogram(quality_returned_box,BinWidth=BW,FaceColor=	"#D95319")
h1=legend({'All systems','Systems returned by Box'});
set( h1,'Box', 'off')
%set( h1,'Box', 'off','Position',[0.14202381798554,0.729960320383784,0.374999991538269,0.137999996423722]);
xlabel('System Quality')
ylabel('Number of Systems')
f = gcf;
%set(gcf, 'Position', [100 100 600 250]); 
ylim([0.9,1000])
%savefig('Pareto Opt\Quality of Returned Systems_Pareto_GSPRT.fig')
%exportgraphics(f,'Pareto Opt\Quality of Returned Systems_Pareto_Box.png','Resolution',600)

figure
histogram(quality_all,BinWidth=BW,FaceColor=[96 96 96]/255)
set(gca, 'YScale', 'log')
hold on
histogram(quality_returned_GSPRT,BinWidth=BW,FaceColor=	"#7E2F8E")
h1=legend({'All systems','Systems returned by GSPRT'});
set( h1,'Box', 'off')
%set( h1,'Box', 'off','Position',[0.14202381798554,0.729960320383784,0.374999991538269,0.137999996423722]);
xlabel('System Quality')
ylabel('Number of Systems')
f = gcf;
%set(gcf, 'Position', [100 100 600 250]); 
ylim([0.9,1000])
%savefig('Pareto Opt\Quality of Returned Systems_Pareto_GSPRT.fig')
%exportgraphics(f,'Pareto Opt\Quality of Returned Systems_Pareto_GSPRT.png','Resolution',600)
%%
% figure
% [n, xout]=histogram(quality_all,BinWidth=BW,FaceColor=[96 96 96]/255);
% bar(xout, n, 'barwidth', 1, 'basevalue', 1);
% f = gcf;
% set(gcf,'YScale','log')
%% The figured 9 of the manuscript
%The nested histogram of quality of one random returned subset 
figure
p1=histogram(quality_all,BinWidth=BW,FaceColor=[200 200 200]/255);
set(gca, 'YScale', 'log','Box','off')
hold on
p2=histogram(quality_returned_box,BinWidth=BW,FaceColor=	"#D95319",FaceAlpha=1);
hold on
p3=histogram(quality_returned_GSPRT,BinWidth=BW,FaceColor=	"#7E2F8E",FaceAlpha=1);
h1=legend([p3,p2,p1],{'$\mathcal{S}^{\mathrm{GSPRT}}$','$\mathcal{S}^{\mathrm{FOSSA}}\backslash \mathcal{S}^{\mathrm{GSPRT}}$','Other systems'},'Interpreter','latex',FontSize=10);
set( h1,'Box', 'off')
%set( h1,'Box', 'off','Position',[0.14202381798554,0.729960320383784,0.374999991538269,0.137999996423722]);
xlabel('Distance to Acceptability')
ylabel('Number of Systems')
f = gcf;
set(gcf, 'Position', [100 100 600 250]); 
ylim([0.85,1000])
savefig('Pareto Opt\Quality of Returned Systems_Pareto.fig')
exportgraphics(f,'Pareto Opt\Quality of Returned Systems_Pareto.png','Resolution',600)
%%
disp('Ellipsoid')
quantile(Ellipse,0.1)
mean(Ellipse)
quantile(Ellipse,0.9)
disp('Box')
quantile(Box,0.1)
mean(Box)
quantile(Box,0.9)
disp('Half')
quantile(Half,0.1)
mean(Half)
quantile(Half,0.9)
disp('GSPRT')
quantile(GSPRT,0.1)
mean(GSPRT)
quantile(GSPRT,0.9)

%%
return_amount_GSPRT=sum(exp_result_GSPRT,2);

disp('Ellipsoid')
quantile(return_amount_ellipse,0.1)
mean(return_amount_ellipse)
quantile(return_amount_ellipse,0.9)
disp('Box')
quantile(return_amount_box,0.1)
mean(return_amount_box)
quantile(return_amount_box,0.9)
disp('Half')
quantile(return_amount_half,0.1)
mean(return_amount_half)
quantile(return_amount_half,0.9)
disp('GSPRT')
quantile(return_amount_GSPRT,0.1)
mean(return_amount_GSPRT)
quantile(return_amount_GSPRT,0.9)
%%

disp('Ellipsoid')
quantile(time_result_ellipse,0.1)
mean(time_result_ellipse)
quantile(time_result_ellipse,0.9)
disp('Box')
quantile(time_result_box,0.1)
mean(time_result_box)
quantile(time_result_box,0.9)
disp('Half')
quantile(time_result_half,0.1)
mean(time_result_half)
quantile(time_result_half,0.9)
%%
for k = 1:length(folders)
    rmpath(folders{k})
end
