function [dataRow, statNames, statEndIdxs] = getStatsMatRow(subsetCathInfo, GT_subset)
    %selects the features to return and returns info about those selected
    %and how to index them, which is used in the python files training the
    %models
    
    numPoints = false;
    posEigVecs = true;
    posEigVals = false;
    posLT = true;
    orientEigVecs = false;
    orientEigVals = false;
    orientLT = false;
    varPos = false;
    varOrient = false;
    meanPos = true;
    meanOrient = true;
    lenDiag = true;
    meanContact = false;
        
    %initialise 
    dataRow = zeros(1,48); 
    statNames = {1,14};
    statNum = 0;
    idx = 0;
    
    %number of points
    if numPoints
        statNum = statNum + 1;
        dataRow(idx+1) = size(subsetCathInfo,1);
        idx = idx + 1;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'n';
    end
    
    if posEigVecs || posEigVals || posLT
        posCov = cov(subsetCathInfo(:,1:3));
        [D,V] = eig(posCov);
    end
    
    %eigenvectors of position cov mat
    if posEigVecs
        statNum = statNum + 1;
        dataRow(idx+1:idx+9) = reshape(D,1,9); 
        idx = idx + 9;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'posEigVecs';
    end
    
    %eigenvalues of position cov mat
    if posEigVals
        statNum = statNum + 1;
        dataRow(idx+1: idx+3) = diag(V); 
        idx = idx + 3;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'posEigVals';
    end
    
    %lower triangular of position cov mat
    if posLT
        statNum = statNum + 1;
        dataRow(idx+1: idx+3) = posCov(logical(tril(ones(3),-1))); 
        idx = idx + 3;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'posLT';
    end
    
    if orientEigVecs || orientEigVals || orientLT
        orientCov = cov(subsetCathInfo(:,4:6));
        [D,V] = eig(orientCov);
    end 
    
    %eigenvectors of orietatation cov mat
    if orientEigVecs
        statNum = statNum + 1;
        dataRow(idx+1:idx+9) = reshape(D,1,9);
        idx = idx + 9;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'orientEigVecs';
    end
    
    %eigenvalues of orietatation cov mat
    if orientEigVals
        statNum = statNum + 1;
        dataRow(idx+1:idx+3) = diag(V);
        idx = idx + 3;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'orientEigVals';
    end
    
    %lower triangular of orietatation cov mat
    if orientLT
        statNum = statNum + 1;
        dataRow(idx+1:idx+3) = orientCov(logical(tril(ones(3),-1)));
        idx = idx + 3;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'orientLT';
    end
        
    %var position
    if varPos
        statNum = statNum + 1;
        dataRow(idx+1:idx+3) =  var(subsetCathInfo(:,1:3)); 
        idx = idx + 3;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'varPos';
    end
    
    %var orientation
    if varOrient
        statNum = statNum + 1;
        dataRow(idx+1:idx+3) =  var(subsetCathInfo(:,4:6)); 
        idx = idx + 3;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'varOrient';
    end

    %mean pos 
    if meanPos
        statNum = statNum + 1;
        dataRow(idx+1:idx+3) = mean(subsetCathInfo(:,1:3)); 
        idx = idx + 3;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'meanPos';
    end
    
    %mean orientation
    if meanOrient
        statNum = statNum + 1;
        dataRow(idx+1:idx+3) = mean(subsetCathInfo(:,4:6)); 
        idx = idx + 3;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'meanOrient';
    end
    
    %length of the diagonal of the minimum enclosing box for the cath points
    if lenDiag
        [~, ~, diagonal] = minEnclosingBox(subsetCathInfo(:,1:3));
        statNum = statNum + 1;
        dataRow(idx+1) = diagonal; 
        idx = idx + 1;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'lenDiag';
    end
    
    %mean contact field of the catheter
    if meanContact
        statNum = statNum + 1;
        dataRow(idx+1) = mean(subsetCathInfo(:,7)); 
        idx = idx + 1;
        statEndIdxs(statNum) = idx;
        statNames{statNum} = 'meanContact';
    end
    
    %Ground truth centering position 
    statNum = statNum + 1;
%     dataRow(idx+1:idx+3) = findCenter(GT_subset(:,1:3));
    dataRow(idx+1:idx+3) = mean(GT_subset(:,1:3));
    idx = idx + 3;
    statEndIdxs(statNum) = idx;
    statNames{statNum} = 'GT_CP';
    
%     statNames = cell(1,numStats);
%     statEndIdxs = zeros(1,numStats);

    dataRow = dataRow(1:idx);
    statNames = statNames(1:statNum);
    statEndIdxs = statEndIdxs(1:statNum);
end