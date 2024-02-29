function [maxM, maxP] = matrixFinder(channelNum, dataNum, pPackageSuccess)
warning('off','all')
tic;

n = channelNum; m = dataNum; pSuccess = pPackageSuccess;

temp = zeros(m, n);
temp(1,1) = 1;

pStack = zeros(m,n);

calculatedLevel = 1;
mMax = temp;
pMax = 0;

while any(any(temp)) ~= 0
    p = zeros(1,n);
    if gfrank(temp) == n
        for stackLevel = calculatedLevel : m
            if length(temp(stackLevel,(temp(stackLevel,:) > 0))) == 1
                p = temp(stackLevel, :) * pSuccess * (1-pSuccess)^(stackLevel - 1) / max(temp(stackLevel,:));
            else
                p = zeros(1,n);
            end
            for i = 1:(stackLevel - 1)
                selectedRows = nchoosek(1:(stackLevel - 1), i);
                for j = 1:length(selectedRows(:,1))
                    selectedRow = selectedRows(j,:);
                    receivedMatrix = temp([stackLevel selectedRow],:);
                    
                    solvedVector = zeros(1, n);
                    for k = 1:n
                        targetVector = zeros(1, n);
                        targetVector(k) = 1;
                        [x, vld] = gflineq(receivedMatrix', targetVector');
                        if vld
                            solvedVector(k) = 1;
                        end
                    end
                    
                    p = p + (solvedVector) * (pSuccess^(i + 1) * (1-pSuccess)^(stackLevel - (i + 1))); 
                end
            end
            if (stackLevel ~= 1)
                p = p + (1-pSuccess) * pStack(stackLevel - 1, :);
            end
            pStack(stackLevel,:) = p;
            calculatedLevel = stackLevel + 1;
        end
        pSum = sum(p);
        if pSum > pMax
            pMax = pSum;
            mMax = temp;
        end
    end
    
    for i = m:-1:1
        if i ~= 1
            if temp(i,:) == temp(i-1,:)
                temp(i,:) = zeros(1,n);
                continue;
            end
        end
        for j = 1:n
            if temp(i,j) == 0
                temp(i,j) = 1;
                calculatedLevel = min(i,calculatedLevel);
                break;
            else
                temp(i,j) = 0;
            end
        end
        break;
    end
end

maxM = mMax;
maxP = pMax;
end