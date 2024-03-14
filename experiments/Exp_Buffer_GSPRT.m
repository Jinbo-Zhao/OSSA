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
alpha=0.05;
%Using the average of all data point as a true mean

record(:,1,:)=1-record(:,1,:);
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
%The Following part should be run one a parallel computing enviroment. 
% It will takes approximately 10 hours to finish
crunch_cluster = parcluster;
%pool_obj = parpool(crunch_cluster);
maxNumCompThreads(48);%change it according to the number available cores

%for macro_rep=1:macro_replications
parfor (macro_rep=1:macro_replications, crunch_cluster)

    %macro_rep
    %collecting outputs for this macroreplication
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
    

    system_info_temp1 = mat2cell(mean_temp,ones(system_amount,1));
    system_info_temp2 = mat2cell(varaiance_temp,num_responses*ones(system_amount,1),[num_responses]);
    system_info=[system_info_temp1,system_info_temp2];
    
    %tic
    return_index_GSPRT = GSPRT(system_info,sample_size_vec,alpha);
    %time_ellipse=toc
    
    
    exp_result_GSPRT(macro_rep,:)=return_index_GSPRT;
   
    %time_result_GSPRT(macro_rep,1)=time_GSPRT;

end
%%
for k = 1:length(folders)
    rmpath(folders{k})
end
save("Exp_Result_GSPRT.mat",'exp_result_GSPRT')