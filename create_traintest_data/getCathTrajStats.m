function [ cathStats, statNames, statEndIdxs, trajectoriesMat ] = getCathTrajStats( ds_origin_start, ds_GT_cathInfo, numStats, row, rots )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

numCathPoints = size(ds_GT_cathInfo,1);
minSampleSize = 20;
sampleSizeDiff = 10; % n+1 set of samples contain ~3% of total cath points more than n set
%all as percentages of total num points
sampleSteps = inf; %difference between starting positions of samples of the same size
                   %inf means only take samples from the trajectory, start
                   %to finish and do not alter the starting point to
                   %collect more samples
peturbation = 0.00; %random perturbation of start and end of the samples to vary the sizes 

numRows = ceil(numCathPoints/sampleSizeDiff); %determine size of matrix to store the trajectory paths in 
cathStats = zeros(numRows,numStats); 

trajectoriesMat = zeros(numCathPoints+2,8000);
longestTraj = 0;
trajCol = -2;
statRow = 0;

if row==1 && rots==1
    cathStats(1,1:4) = [minSampleSize sampleSteps sampleSizeDiff peturbation ];
    trajectoriesMat(1,1:4) = [ minSampleSize sampleSteps sampleSizeDiff peturbation ];
    statRow = statRow + 1;
end
            
for gap=minSampleSize:sampleSizeDiff:numCathPoints
    peturb = ceil(gap*numCathPoints*peturbation);
    for p=minSampleSize:sampleSteps:numCathPoints-gap 
        %start = ceil(p*numCathPoints);
        start = p;
        if p+gap < numCathPoints %&& start <= lastStartIndex
            statRow = statRow + 1;
            sidx = max( [start+randi([-peturb,peturb],1) 1] );
            eidx = min( [start+gap+randi([-peturb,peturb],1) numCathPoints]);
            %if all percentages
            %eidx = min([start+ceil(gap*numCathPoints)+randi([-peturb,peturb],1) numCathPoints]); 
            n = eidx-sidx;
            cathStats(statRow,3:4) = [sidx eidx];
            [cathStatsRow,statNames,statEndIdxs] = getStatsMatRow(ds_origin_start(sidx:eidx,:),ds_GT_cathInfo(sidx:eidx,:));
            cathStats(statRow,5:statEndIdxs(end)+4) = cathStatsRow;
            trajCol = trajCol + 3;
            trajectoriesMat(2,trajCol) = n;
            trajectoriesMat(3:n+3,trajCol:trajCol+2) = ds_GT_cathInfo(sidx:eidx,1:3);
            if n > longestTraj
                longestTraj = n;
            end
        else
            break
        end
    end
trajectoriesMat = trajectoriesMat(1:longestTraj+2, 1:trajCol+2);
cathStats = cathStats(1:statRow,1:statEndIdxs(end)+4);
end

