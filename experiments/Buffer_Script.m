clc
clear
close all
%% Generating all systems
Total_Buffer=20;
node_amount=[4,2,2,1];
Max_each_satge=Total_Buffer./node_amount;

solution_list=[];
for rep1=0:1:Max_each_satge(1)
    temp_solution=zeros(1,4);
    temp_solution(1,1)=rep1;
    for rep2=0:1:Max_each_satge(2)
        temp_solution(1,2)=rep2; 
        for rep3=0:1:Max_each_satge(3)
            temp_solution(1,3)=rep3; 
            for rep4=0:1:Max_each_satge(4)
                temp_solution(1,4)=rep4; 
                if  temp_solution*node_amount'>20
                    continue
                else
                    solution_list=[solution_list;temp_solution];
                end
            end
        end
    end
end
%%
crunch_cluster = parcluster;
%pool_obj = parpool(crunch_cluster);
maxNumCompThreads(48);

runlength=50000;
num_solution=length(solution_list);
tic
record=zeros(num_solution,2,runlength);

%for rep=1:num_solution
parfor (rep=1:num_solution, crunch_cluster)
    x=solution_list(rep,:);
    x = [x(1),x(1),x(1),x(1),x(2),x(3),x(2),x(3),x(4)];
    solution = [x, 20-sum(x)]; %buffer size of the 10 nodes
    seed=321+rep;
    [AvgUtility,~,AverageWaitingTime] = Buffer_Simulation(solution, runlength, seed);
    temp=[AvgUtility';AverageWaitingTime];
    record(rep,:,:)=temp;
end

save("Simulation_Data.mat")