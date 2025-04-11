function [return_index]=Plotting_confidence_box_updated(system_info,sample_size_vec,alpha,color)

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


if all(sample_size_vec==sample_size_vec(1)*ones(num_systems,1))
    sample_size=sample_size_vec(1);

    if num_responses==2
        closed_box_t_critical_value=tinv(1-((1-splited_confidence)/2),sample_size-1);
    elseif num_responses>=3
        closed_box_t_critical_value=tinv(1-((1-splited_confidence)/(2*num_responses)),sample_size-1);
    end

    unique_elipsoid_flag=1;
else
    unique_elipsoid_flag=0;
end

for rep =1:num_systems
    sample_mean=system_info{rep,1}.';
    sample_covariance_matrix=system_info{rep,2};
    sample_size=sample_size_vec(rep);
    sample_variances = diag(sample_covariance_matrix);

    if unique_elipsoid_flag==0
        if num_responses==2
            closed_box_t_critical_value=tinv(1-((1-splited_confidence)/2),sample_size-1);
        elseif num_responses>=3
            closed_box_t_critical_value=tinv(1-((1-splited_confidence)/(2*num_responses)),sample_size-1);
        end   
    end
    closed_box_radius=closed_box_t_critical_value*sqrt(sample_variances/sample_size);
    closed_box_right_up_corner=sample_mean+closed_box_radius;
    closed_box_left_down_corner=sample_mean-closed_box_radius;
    worst_corner_storage(rep,:)=closed_box_right_up_corner;
    best_corner_storage(rep,:)=closed_box_left_down_corner;
end

[membership,worst_corner_pareto]=find_pareto_frontier(worst_corner_storage);
return_index=zeros(num_systems,1);


for rep =1:num_systems
    if membership(rep)==1
        return_index(rep,1)=1;%The system shold be reterned and end the loop immediately
        continue
    end

    %Using best-corner to precheck the overlap. If best corner is
    %dominated, the system must be screened out
    BHS=best_corner_storage(rep,:);
    num_pareto=sum(membership);
    return_index(rep,1)=1;%the system's best corner is not dominated by any other systems's worst corner
    for pareto_rep=1:num_pareto
        check=(BHS>=worst_corner_pareto(pareto_rep,:));
        if sum(check)==num_responses%the BHS already across the pareto front
            return_index(rep,1)=0;
            break
        end
    end
        
end

figure
for rep =1:num_systems
    if return_index(rep,1)==1
        temp1=worst_corner_storage(rep,:);
        temp2=best_corner_storage(rep,:);
        position=[temp2, temp1-temp2];
        hold on
        rectangle('Position',position,'FaceColor',color,'EdgeColor','none') 
    end
end
max_loc=max(worst_corner_storage); 
for rep =1:num_systems

    temp1=worst_corner_storage(rep,:);
    position=[temp1, max_loc-temp1];
    hold on
    rectangle('Position',position,'FaceColor','w','EdgeColor','w') 

end
 

end