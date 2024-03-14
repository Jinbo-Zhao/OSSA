function [quality] = quality_Pareto(returned_subset,Pareto_front)
%This is a function evaluating a set of returned system compared with a
%given pareto front. We calculate each system's closed distance to the
%pareto front. Specifically, What is the closed response that is nor
%dominated by any of the pareto optimal systems.
[~,Pareto_front]=find_pareto_frontier(Pareto_front);
phantoms = FindPhants(Pareto_front);
[num_phantoms,~]=size(phantoms);
[num_systems,num_responses]=size(returned_subset);
distance_record=zeros(num_systems,num_phantoms);

for sys_rep =1:num_systems
    for phantom_rep=1:num_phantoms
        A_QP=eye(num_responses);
        H=2*A_QP;
        true_mean=returned_subset(sys_rep,:);
        b_QP=(phantoms(phantom_rep,:)-true_mean).';
        options = optimoptions('quadprog','Display','none');
        [~, f_val, ~] = quadprog(H,[],A_QP,b_QP,[],[],[],[],[],options);
        distance_record(sys_rep,phantom_rep)=sqrt(f_val);
    end
end
quality=min(distance_record,[],2);
end