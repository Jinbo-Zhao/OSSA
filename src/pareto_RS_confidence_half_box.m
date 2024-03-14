%This is function used to screen systems when acceptability is defined as
%pareto Optimality. It takes each system's sampel mean vector and sample
%covariance matrix as the input and returns the index of returned sys. 
%The Common Random Number generator is assumed to be not used, otherwise,
%please update the splited_confidence.
function [return_index] = pareto_RS_confidence_half_box(system_info,sample_size_vec,alpha)

%initialization
[num_systems,~]=size(system_info);
[num_responses,~]=size(cell2mat(system_info(1,2)));
worst_corner_storage=zeros(num_systems,num_responses);
best_corner_storage=zeros(num_systems,num_responses);

if num_responses==2
   splited_confidence=(1-alpha)^(1/(num_systems*2));%the confidence on each system's each responces
elseif num_responses>=3
   splited_confidence=(1-alpha)^(1/(num_systems));%the confidence on each system's each responces
end

%Calculating the critical value for specifying each confidence region, See
%section 3.3 Confidence Regions and Table 1
if all(sample_size_vec==sample_size_vec(1)*ones(num_systems,1))
    sample_size=sample_size_vec(1);

    if num_responses==2
        half_box_t_critical_value=tinv(splited_confidence,sample_size-1);
    elseif num_responses>=3
        half_box_t_critical_value=tinv(1-(1-splited_confidence)/(num_responses),sample_size-1);
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
    % The sample size are not same for all systems, hence, we have to
    % calculate the critical value for each system. 
    if unique_elipsoid_flag==0
        if num_responses==2
            half_box_t_critical_value=tinv(splited_confidence,sample_size-1);
        elseif num_responses>=3
            half_box_t_critical_value=tinv(1-(1-splited_confidence)/(num_responses),sample_size-1);
        end
    end
    half_box_radius=half_box_t_critical_value*sqrt(sample_variances/sample_size);
    half_box_right_up_corner=sample_mean+half_box_radius;
    half_box_left_down_corner=sample_mean-half_box_radius;
    worst_corner_storage(rep,:)=half_box_right_up_corner;
    best_corner_storage(rep,:)=half_box_left_down_corner;
end
% Above content corresponds to the line1-line6 of Algorithm 4(FOSSA (Half-)Box for Pareto Optimality 

%Following correspond to line7- line10.
%Due to transitivity of Pareto Optimality, we only need to compare the
%l_i(best_corner_storage(i,:)) with the efficient set of u_i s.
[membership,worst_corner_pareto]=find_pareto_frontier(worst_corner_storage);
return_index=zeros(num_systems,1);

for rep =1:num_systems
    if membership(rep)==1
        %Since the system's worst corner is not dominated, its best corner
        %must be more preferable than its worst corner.
        return_index(rep,1)=1;%The system shold be reterned and end the loop immediately
        continue
    end

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

end

