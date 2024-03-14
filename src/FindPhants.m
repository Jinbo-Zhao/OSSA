%This Function is used to decompose the region dominated by a set of
%systems. The dominated region will be divided into several hyper unbounded
%rectangle with Phantom systems as the only corner.
%----------------------------------------------------------------------------------
%Applegate, Eric A., et al. "Multi-objective ranking and selection: Optimal
%sampling laws and tractable approximations via SCORE." Journal of Simulation 14.1 (2020): 21-40.

function phantoms = FindPhants(paretos)

maxy=Inf;

% find size of Pareto matrix
[~,numobj] = size(paretos);
phantoms=ones(1,numobj)*0;

v=1:1:numobj;
for b=1:numobj
T = nchoosek(v,b);
for i=1:nchoosek(numobj,b)
        tempPars = paretos(:,T(i,:));
        tempPars = unique(tempPars,'rows');
        % find temp paretos
        [~,tempPars2]=find_pareto_frontier(tempPars);

        %do sweep
        phants = Sweep(tempPars2);

        phan=ones(size(phants,1),numobj)*maxy;
        phan(:,T(i,:)) = phants;
        phantoms=[phantoms; phan];
end
end
phantoms(1,:)=[];
end