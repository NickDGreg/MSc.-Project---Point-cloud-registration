function [gmfit, clusterX] = getAortaDescriptors(aortaverts, centerline, crossSections, crossSection_endIndices, k)
% gmfit is the mixture model struct containing all that statistics of the
% MM. ClusterX is a list indexing each aortavert to a cluster.
% function makes an initial division of the points using a percentage and
% uses GMM to finalise. The initial percentage fitting avoids any bad 
% random initialisations of the GMM such that it can be fully automated.

% need to make into function and call from the simulated and cath stats
% get mean, cov eigs and lower triangular for each MM
% return, does mean a variable size of the stat row which i need
% to incorporate into the createTrainfiles
% check whether making a transformation to the aorta I can make the same
% translation to the aorta and save a lot in processing time.
%

%get mean and cov of aorta quarters to init GMM
numCLpoints = size(centerline,1);
CL_subsets = ones(1,k+1);
for i=1:k-1
    CL_subsets(i+1) = ceil(numCLpoints*(i/k));
end
CL_subsets(k+1) = numCLpoints;

mu = ones(k,3);
CL_fraction = ceil(numCLpoints/k/2);
for i=1:k
    mu(i,:) = centerline(CL_fraction*(k+k-1),:);
end

covar = zeros(3,3,k);
for i=1:k
    pointsAdded = 0;
    crossSubset = zeros(size(aortaverts,1),3);
    for j=CL_subsets(i):CL_subsets(i+1)
        crossSec = crossSections(1:crossSection_endIndices(j),:,j);
        crossSubset(pointsAdded+1:pointsAdded+size(crossSec,1),:) = crossSec;
        pointsAdded = pointsAdded + size(crossSec,1);
    end
    crossSubset = crossSubset(1:pointsAdded,:);
    covar(:,:,i) = cov(crossSubset);
end
resp = ones(1,k)*1/k;

GMM_init = struct('mu',mu,'Sigma',covar,'ComponentProportion',resp);
options = statset('MaxIter',1000); % Increase number of EM iterations
gmfit = fitgmdist(aortaverts,k,'Start',GMM_init,'Options',options);
clusterX = cluster(gmfit,aortaverts);

% groups = zeros(size(aortaverts,1),3*k);
% saveIdx = zeros(1,k);
% colour = 'rgbmyo';
% hold on, grid on, view(-85,11);
% for idx=1:k
%     groupIdxs = clusterX == idx;
%     group = aortaverts(groupIdxs,:);
%     scatter3(group(:,1),group(:,2),group(:,3),1,colour(idx));   
% end    
% plot3(gmfit.mu(:,1),gmfit.mu(:,2),gmfit.mu(:,3),'kx','LineWidth',2,'MarkerSize',5)
% hold off

end
