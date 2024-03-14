%This function is implemented according to the procedure 3 in Wenyu wang's
%Phd dissertation. It is based on pairwise comparison.
%!!This can only be applied to situation where there are 2 responses
%Wang, Wenyu. Sequential Procedures for the” Selection” Problems in Discrete Simulation Optimization. Diss. Purdue University, 2019.
function [return_index] = GSPRT(system_info,sample_size_vec,alpha)

[num_systems,~]=size(system_info);
[num_responses,~]=size(cell2mat(system_info(1,2)));

if all(sample_size_vec==sample_size_vec(1)*ones(num_systems,1))
    sample_size=sample_size_vec(1);
else
    error('Requiring equal sample size')
end

ratio_storage=inf*ones(num_systems,1);
for repi=1:num_systems
    ratio_storage_temp=inf*ones(num_systems,1);
    for repj=1:num_systems
        if repi==repj
            continue
        end
           sample_mean_i=system_info{repi,1}.'; 
           sample_mean_j=system_info{repj,1}.';
           %The sample mean for the pairwise difference
           sample_mean=sample_mean_i-sample_mean_j;

           sample_covariance_matrix_i=system_info{repi,2};
           sample_covariance_matrix_j=system_info{repj,2};
           %The sample coaraince for the pairwise difference
           sample_covariance_matrix=sample_covariance_matrix_i+sample_covariance_matrix_j;
            
            if sum(sample_mean>=0)==num_responses
                solving_part=0;%Don;t have to solve the easy one(down side one:dominated by 0)
            else
                solving_part=1;
            end

            if solving_part==0%Solving the up part of formula 3.20
                H=sample_covariance_matrix^-1;
                
                A1=[-1 0];
                A2=[ 1 0;0 -1];
                b=sample_mean;

                B2=[1;-1];
                B1=-b(1,1);
                B2=B2.*b;
                options = optimoptions('quadprog','Display','none');
                [~, f_val_1, ~] = quadprog(H,[],A1,B1,[],[],[],[],[],options);
                [~, f_val_2, ~] = quadprog(H,[],A2,B2,[],[],[],[],[],options);
                opt_val=min([f_val_1,f_val_2]);

                ratio=-sample_size*opt_val;
                ratio_storage_temp(repj,1)=ratio;

            else%Solving the down side part of formula 3.20
                H=sample_covariance_matrix^-1;
                A=eye(num_responses);
                b=sample_mean;
                options = optimoptions('quadprog','Display','none');
                [~, f_val, ~] = quadprog(H,[],A,b,[],[],[],[],[],options);
                opt_val=f_val;

                ratio=sample_size*opt_val;
                ratio_storage_temp(repj,1)=ratio;
                
            end
        
    end
    ratio_storage(repi,1)=min(ratio_storage_temp);  
end
%The log here is due to our optimization question is corresponding to the
%part in the exp part of Formula 3.20.
return_index=ratio_storage>=log(alpha/num_systems);
end

