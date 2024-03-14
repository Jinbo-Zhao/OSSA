%This is function used to screen systems when acceptability is defined as
%Constrained Optimality. It takes each system's sampel mean vector and sample
%covariance matrix as the input and returns the index of returned sys. 
%The Common Random Number generator is assumed to be not used, otherwise,
%please update the splited_confidence.
function [return_index] = constrained_RS_confidence_box(system_info,sample_size_vec,alpha,feasible)

%initialization
[num_systems,~]=size(system_info);
[num_responses,~]=size(cell2mat(system_info(1,2)));
worst_corner_storage=zeros(num_systems,num_responses);
best_corner_storage=zeros(num_systems,num_responses);

if num_responses==2
   splited_confidence=(1-alpha)^(1/(num_systems*2));%the confidence on each system's each responses
elseif num_responses>=3
   splited_confidence=(1-alpha)^(1/(num_systems));%the confidence on each system's each responses
end


%Calculating the critical value for specifying each confidence region, See
%section 3.3 Confidence Regions and Table 1
if all(sample_size_vec==sample_size_vec(1)*ones(num_systems,1))
    sample_size=sample_size_vec(1);
    
    if num_responses==2
        closed_box_t_critical_value=tinv(1-((1-splited_confidence)/2),sample_size-1);
    elseif num_responses>=3
        closed_box_t_critical_value=tinv(1-((1-splited_confidence)/(2*num_responses)),sample_size-1);
    end

    unique_confidence_flag=1;
else
    unique_confidence_flag=0;
end
% Following for loop corresponds to the the line4-5
%l_i are the best corners
%u_i are the worst corners
for rep =1:num_systems
    sample_mean=system_info{rep,1}.';
    sample_covariance_matrix=system_info{rep,2};
    sample_size=sample_size_vec(rep);
    sample_variances = diag(sample_covariance_matrix);
    % The sample size are not same for all systems, hence, we have to
    % calculate the critical value for each system. 
    if unique_confidence_flag==0
        if num_responses==2
            closed_box_t_critical_value=tinv(1-((1-splited_confidence)/2),sample_size-1);
        elseif num_responses>=3
            closed_box_t_critical_value=tinv(1-((1-splited_confidence)/(2*num_responses)),sample_size-1);
        end
    end
    closed_box_radius=closed_box_t_critical_value*sqrt(sample_variances/sample_size);%A d-dimension vector
    closed_box_right_up_corner=sample_mean+closed_box_radius;
    closed_box_left_down_corner=sample_mean-closed_box_radius;
    worst_corner_storage(rep,:)=closed_box_right_up_corner;
    best_corner_storage(rep,:)=closed_box_left_down_corner;
end
%Following corresponds to the line 6-9
primary_lower_bound=zeros(num_systems,1);
primary_upper_bound=zeros(num_systems,1);
for rep =1:num_systems
    if sum(best_corner_storage(rep,2:end)>feasible)==0
        primary_lower_bound(rep,1)=best_corner_storage(rep,1);
    else
        primary_lower_bound(rep,1)=inf;
    end
        
    if sum(worst_corner_storage(rep,2:end)>feasible)==0
        primary_upper_bound(rep,1)=worst_corner_storage(rep,1);
    else
        primary_upper_bound(rep,1)=inf;
    end  
end
% Line 10-13
threshold=min(primary_upper_bound);
return_index=primary_lower_bound<threshold;

temp=primary_lower_bound==inf;
return_index=return_index.*(~temp);

end

