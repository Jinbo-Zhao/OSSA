%This is function used to screen systems when acceptability is defined as
%pareto Optimality. It takes each system's sampel mean vector and sample
%covariance matrix as the input and returns the index of returned sys. 
%The Common Random Number generator is assumed to be not used, otherwise,
%please update the splited_confidence.
function [return_index] = constrained_RS_confidence_ellipse(system_info,sample_size_vec,alpha,feasible)

%initialization
[num_systems,~]=size(system_info);
[num_responses,~]=size(cell2mat(system_info(1,2)));
worst_corner_storage=zeros(num_systems,num_responses);
best_corner_storage=zeros(num_systems,num_responses);

F_splited_confidence=(1-alpha)^(1/num_systems);

%Calculating the critical value for specifying each confidence region, See
%section 3.3 Confidence Regions and Table 1
if all(sample_size_vec==sample_size_vec(1)*ones(num_systems,1))
    sample_size=sample_size_vec(1);
    F_crtical_value=finv(F_splited_confidence,num_responses,sample_size-num_responses);
    ellipsoid_radius=(num_responses*(sample_size-1)*F_crtical_value)/(sample_size*(sample_size-num_responses));
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
       F_crtical_value=finv(F_splited_confidence,num_responses,sample_size-num_responses);
       ellipsoid_radius=(num_responses*(sample_size-1)*F_crtical_value)/(sample_size*(sample_size-num_responses));
    end

    %Each ellipse can be charicterized by the samllest box containning it
    %We name it Ellipse_Box here
    ellipse_box_radius=sqrt(ellipsoid_radius*sample_variances);
    ellipse_box_right_up_corner=sample_mean+ellipse_box_radius;
    ellipse_box_left_down_corner=sample_mean-ellipse_box_radius;
    worst_corner_storage(rep,:)=ellipse_box_right_up_corner;
    best_corner_storage(rep,:)=ellipse_box_left_down_corner;
end
% Following for loop corrsponds to line3-6
primary_lower_bound=zeros(num_systems,1);
primary_upper_bound=zeros(num_systems,1);
for rep =1:num_systems
    sample_mean=system_info{rep,1}.';
    sample_covariance_matrix=system_info{rep,2};
    sample_size=sample_size_vec(rep);

    if unique_elipsoid_flag==0
       F_crtical_value=finv(F_splited_confidence,num_responses,sample_size-num_responses);
       ellipsoid_radius=(num_responses*(sample_size-1)*F_crtical_value)/(sample_size*(sample_size-num_responses));
    end

    %For getting the lower bound(l_i) of primary response
    %Following line is try to pre-check the overlap problem in line 4
    if sum(best_corner_storage(rep,2:end)>feasible)>0%best corner is NOT in the feasible region
        primary_lower_bound(rep,1)=inf;

    else %We need solve the quadratic programming problem to determine
        A_QP=eye(num_responses);
        A_QP=A_QP(2:end,:);
        b_QP=(feasible-sample_mean(2:end,1));
        options = optimoptions('quadprog','Display','none');
        H=2*sample_covariance_matrix^-1;
        [~, f_val, exitflag] = quadprog(H,[],A_QP,b_QP,[],[],[],[],[],options);
        
        if exitflag~=1 %in case the solver fails. Switch the method used to solve the quadratic programming
            options = optimoptions('quadprog','Display','none','algorithm','active-set','MaxIterations',100000);
            [~, f_val, exitflag] = quadprog(H,[],A_QP,b_QP,[],[],[],[],zeros(num_responses,1),options);
            
            if exitflag~=1%The optimization problem cannot be solved
                f_val=0;%For conservaness, let this system pass the test
            end

        end

        if f_val>ellipsoid_radius
           primary_lower_bound(rep,1)=inf;

        else%have to solve the Quadraticly constrained Linear Opt defined in line 4
           lower_bound = QCLP(sample_covariance_matrix,ellipsoid_radius,sample_mean,feasible,num_responses);
           primary_lower_bound(rep,1)=lower_bound;
        end
    end
    
    %For getting the upper bound(u_i) for primary response
    %Following line is try to pre-check the overlap problem in line 5
    if sum(worst_corner_storage(rep,2:end)>feasible)==0 %worst corner is in the feasible region
        primary_upper_bound(rep,1)=worst_corner_storage(rep,1);
    else
        primary_upper_bound(rep,1)=inf;
    end
        

end
%Following correspond to line 6-9
threshold=min(primary_upper_bound);
return_index=primary_lower_bound<threshold;

temp=primary_lower_bound==inf;
return_index=return_index.*(~temp);
end

function [f_val] = QCLP(sample_covariance_matrix,ellipsoid_radius,sample_mean,feasible,num_responses)
    %This is the case that if we ignore the linear constraints, minimize the
    %primary res, then the optimal solution does not violate the linear
    %constrainrts. We found the solution.
    temp=sample_covariance_matrix(:,1);
    opt_sol=-sqrt(ellipsoid_radius/sample_covariance_matrix(1,1))*temp;%The point on the ellipsoid that minimized the primary response
    opt_sol=opt_sol+sample_mean;

    if opt_sol(2:num_responses,1)<=feasible
        f_val=opt_sol(1,1);
        return
    end

    %The trick failed, using cone programming to solve the problem
    %https://www.mathworks.com/help//optim/ug/convert-quadratic-constraints-to-second-order-cone-constraints.html

    Q=sample_covariance_matrix^-1;
    S=sqrtm(Q);
    b=zeros(num_responses,1);
    d=zeros(num_responses,1);
    gamma=-sqrt(ellipsoid_radius);
    sc = secondordercone(S,b,d,gamma);
    f=[1;zeros(num_responses-1,1)];
    options = optimoptions('coneprog','Display','none','ConstraintTolerance',1e-15,'OptimalityTolerance',1e-15);
    A=eye(num_responses);
    A=A(2:end,:);
    b=(feasible-sample_mean(2:end,1));
    [~,f_val] = coneprog(f,sc,A,b,[],[],[],[],options);
    f_val=f_val+sample_mean(1,1);

end