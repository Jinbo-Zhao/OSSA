clc
clear
close all
folders = strcat({'..\'}, {'data','src'});
for k = 1:length(folders)
    addpath(folders{k})
end
load('Simulation_Data.mat')
%%

system_amount=1001;
num_responses=2;
sample_amount=50;
macro_replications=1000;
%Using the average of all data point as a true mean

record(:,1,:)=1-record(:,1,:);%Change the utility rate to idel rate
overall_data=record;
true_mean= mean(overall_data,3);

for system_index=1:system_amount
    temp=overall_data(system_index,:,:);
    temp=reshape(temp,num_responses,50000);
    sample_covariance_matrix=cov(temp');
    true_varaiance(num_responses*system_index-num_responses+1:num_responses*system_index,:)=sample_covariance_matrix;
end
system_info_temp1 = mat2cell(true_mean,ones(system_amount,1));
system_info_temp2 = mat2cell(true_varaiance,num_responses*ones(system_amount,1),[num_responses]);
true_system_info=[system_info_temp1,system_info_temp2];

%%


for macro_rep=1:macro_replications
    macro_rep
    record=overall_data(:,:,(macro_rep-1)*sample_amount+1:macro_rep*sample_amount);
    varaiance_temp=zeros(system_amount*num_responses,num_responses);
    mean_temp=zeros(system_amount,num_responses);
    for system_index=1:system_amount
        data_temp=record(system_index,:,:);
        data_temp=reshape(data_temp,2,sample_amount);
        mean_temp(system_index,:)=mean(data_temp,2);
        varaiance_temp(num_responses*system_index-num_responses+1:num_responses*system_index,:)=cov(data_temp');
    end

    sample_size_vec=sample_amount*ones(system_amount,1);
    alpha=0.05;
    system_info_temp1 = mat2cell(mean_temp,ones(system_amount,1));
    system_info_temp2 = mat2cell(varaiance_temp,num_responses*ones(system_amount,1),[num_responses]);
    system_info=[system_info_temp1,system_info_temp2];
    
    tic
    return_index_ellipse = pareto_RS_confidence_ellipse(system_info,sample_size_vec,alpha);
    time_ellipse=toc;
    
    tic
    return_index_box = pareto_RS_confidence_box(system_info,sample_size_vec,alpha);
    time_box=toc;
    
    tic
    return_index_half = pareto_RS_confidence_half_box(system_info,sample_size_vec,alpha);
    time_half=toc;
    
    exp_result_ellipse(macro_rep,:)=return_index_ellipse;
    exp_result_box(macro_rep,:)=return_index_box;
    exp_result_half(macro_rep,:)=return_index_half;
    
    time_result_ellipse(macro_rep,1)=time_ellipse;
    time_result_box(macro_rep,1)=time_box;
    time_result_half(macro_rep,1)=time_half;
end

%%
frquency_ellipse=sum(exp_result_ellipse);
frquency_box=sum(exp_result_box);
frquency_half=sum(exp_result_half);

return_amount_ellipse=sum(exp_result_ellipse,2);
return_amount_box=sum(exp_result_box,2);
return_amount_half=sum(exp_result_half,2);

[membership, member_value]=find_pareto_frontier(true_mean);

subplot(1,3,1)
scatter(true_mean(:,1),true_mean(:,2),[],frquency_ellipse,'filled')
colorbar
hold on
scatter(member_value(:,1),member_value(:,2),[],'r')
title('Returned Frequency Ellipse')

subplot(1,3,2)
scatter(true_mean(:,1),true_mean(:,2),[],frquency_box,'filled')
colorbar
title('Returned Frequency Box')
hold on
scatter(member_value(:,1),member_value(:,2),[],'r')

subplot(1,3,3)
scatter(true_mean(:,1),true_mean(:,2),[],frquency_half,'filled')
colorbar
title('Returned Frequency Half')
hold on
scatter(member_value(:,1),member_value(:,2),[],'r')

figure
histogram(return_amount_ellipse)
hold on
histogram(return_amount_box)
hold on
histogram(return_amount_half)
title('Size of Returned Subset')

%% Quality figure
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
quality_returned=quality_all(logical(return_index_ellipse),:);

figure
histogram(quality_all,BinWidth=0.01)
hold on
histogram(quality_returned,BinWidth=0.01)

%% Returned frequency figure




ellipsoid_radius=0.01;

[~, member_value]=find_pareto_frontier(true_mean);
uncertainty_scatter(true_system_info,ellipsoid_radius,frquency_ellipse)
hold on
scatter(member_value(:,1),member_value(:,2),10,'red','filled')
title('Returned Frequency Ellipse')

uncertainty_scatter(true_system_info,ellipsoid_radius,frquency_box)
hold on
scatter(member_value(:,1),member_value(:,2),10,'red','filled')
title('Returned Frequency Box')

uncertainty_scatter(true_system_info,ellipsoid_radius,frquency_half)
hold on
scatter(member_value(:,1),member_value(:,2),10,'red','filled')
title('Returned Frequency Half')
%%
for k = 1:length(folders)
    rmpath(folders{k})
end
save("Exp_Env_Data_Pareto_full.mat")
