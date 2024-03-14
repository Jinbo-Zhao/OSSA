%This is function used to screen systems when acceptability is defined as
%pareto Optimality. It takes each system's sampel mean vector and sample
%covariance matrix as the input and returns the index of returned sys. 
%The Common Random Number generator is assumed to be not used, otherwise,
%please update the F_splited_confidence.
function [return_index] = pareto_RS_confidence_ellipse(system_info,sample_size_vec,alpha)

%initialization
[num_systems,~]=size(system_info);
[num_responses,~]=size(cell2mat(system_info(1,2)));
worst_corner_storage=zeros(num_systems,num_responses);
best_corner_storage=zeros(num_systems,num_responses);

%define the ellipse
F_splited_confidence=(1-alpha)^(1/num_systems);


%checking if the sample size are the same. If systems have the same
%sample size, the F crtical value is unique to avoid redundant computation
if all(sample_size_vec==sample_size_vec(1)*ones(num_systems,1))
    sample_size=sample_size_vec(1);
    F_crtical_value=finv(F_splited_confidence,num_responses,sample_size-num_responses);
    ellipsoid_radius=(num_responses*(sample_size-1)*F_crtical_value)/(sample_size*(sample_size-num_responses));
    unique_elipsoid_flag=1;
else
    unique_elipsoid_flag=0;
end
%Following for loop is corressponding to the line 3-5
%w_j is named worst_corner here
%Best corners will be used to pre-solve the optmization problems defined in
%Line7
for rep =1:num_systems
    sample_mean=system_info{rep,1}.';
    sample_covariance_matrix=system_info{rep,2};
    sample_size=sample_size_vec(rep);
    sample_variances = diag(sample_covariance_matrix);
    % The sample size are not same for all systems, hence, we have to
    % calculate the critical value for each system. 
    if unique_elipsoid_flag==0
        F_crtical_value=finv(F_splited_confidence,num_responses,sample_size-num_responses);
        ellipsoid_radius=(num_responses*(sample_size-1)*F_crtical_value)/(sample_size*(sample_size-num_responses));
    end
    %Each ellipse can be charicterized by the samllest box containning it
    %We name it Ellipse_Box here
   
    ellipse_box_radius=sqrt(ellipsoid_radius*sample_variances);
    ellipse_box_right_up_corner=sample_mean+ellipse_box_radius;
    ellipse_box_left_down_corner=sample_mean-ellipse_box_radius;
    %Worst corners are the Red dots in the figure 3
    worst_corner_storage(rep,:)=ellipse_box_right_up_corner;
    %These will be used to provide a relaxation of overlapping-checking
    %problem which will be solved later.
    best_corner_storage(rep,:)=ellipse_box_left_down_corner;
end
%Decompose the W(\cap{w_j})^c in to multiple open "half-boxes"
[membership,worst_corner_pareto]=find_pareto_frontier(worst_corner_storage);
phantoms = FindPhants(worst_corner_pareto);
phantoms=flip(phantoms);
return_index=zeros(num_systems,1);
[num_phantoms,~]=size(phantoms);

%Following correspond to the line6-7
for rep =1:num_systems
    
    %System's right-up corner is not dominated by any system's right-up corner, thus its
    %confidence region itself must not dominated by any other system's right-up corner
    if membership(rep)==1
        return_index(rep,1)=1;%The system shold be reterned and end the loop immediately
        continue
    end

    %Using best-corner to precheck the overlap. If best corner is
    %dominated, the system must be screened out
    BHS=best_corner_storage(rep,:);
    continue_flag=1;%if this is never changed, BHS is on the inner(worse) side of pareto front
    for phantom_rep=1:num_phantoms
        check=(BHS<=phantoms(phantom_rep,:));
        if sum(check)==num_responses%the BHS already across the pareto front -- we will do quadratic opt problem
            continue_flag=0;
            break
        end
    end
    
    if continue_flag==1%The best corner is dominated, screen it out
        return_index(rep,1)=0;
        continue
    end
    %All pre-check fails, have to do the optimization to check the overlap
    sample_mean=system_info{rep,1};
    sample_covariance_matrix=system_info{rep,2};
    sample_size=sample_size_vec(rep);
    if unique_elipsoid_flag==0
         F_crtical_value=finv(F_splited_confidence,num_responses,sample_size-num_responses);
         ellipsoid_radius=(num_responses*(sample_size-1)*F_crtical_value)/(sample_size*(sample_size-num_responses));
    end
    %Check the overlap of the confidence ellipsoid and the the multiple
    %open 'half-boxes.
    for phantom_rep=1:num_phantoms
        
        A_QP=eye(num_responses);
        b_QP=(phantoms(phantom_rep,:)-sample_mean).';
        options = optimoptions('quadprog','Display','none');
        H=2*sample_covariance_matrix^-1;
        [~, f_val, exitflag] = quadprog(H,[],A_QP,b_QP,[],[],[],[],[],options);
        
        if f_val<=ellipsoid_radius
            return_index(rep,1)=1;%The system shold be reterned and end the loop immediately
            break
        end

        if exitflag~=1 %in case the solver fails. Switch the method used to solve the quadratic programming
        options = optimoptions('quadprog','Display','none','algorithm','active-set','MaxIterations',100000);
       
        if max(b_QP)==inf
            save_index=ones(num_responses,1);
            for check_rep=1:num_responses
                if b_QP(check_rep)==inf
                   save_index(check_rep,1)=0;
                end
            end
            A_QP=A_QP(logical(save_index),:);
            b_QP=b_QP(logical(save_index),:);
        end 
        
        [~, f_val, exitflag] = quadprog(H,[],A_QP,b_QP,[],[],[],[],zeros(num_responses,1),options);

            if f_val<=ellipsoid_radius
                return_index(rep,1)=1;%The system shold be reterned
                break
            end
            
            if exitflag~=1
                return_index(rep,1)=1;
            end

        end

    end
    

end

end

