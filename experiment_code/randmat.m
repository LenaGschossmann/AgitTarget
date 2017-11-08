
function block_order = randmat(subject_id)

conditions_wS = [1 2 3 4];
nblocks = 2; %number of blocks per each condition

randlist = zeros(factorial(4), 8);

%create pseudorandom blocklist

iniVec = [1 2 3 4];
randIdx = 1;

%shuffle iniVec
tmpVec = iniVec;
for idx1 = 1:3
    for idx2 = 1:2 
        for idx3 = 1:4
            %first 4 blocks
            randlist(randIdx, 1:4) = tmpVec;
                        
            %second 4 blocks: 1 2 3 4 | 3 4 2 1
            randlist(randIdx, 5:6) = tmpVec(3:4);
            randlist(randIdx, 7) = tmpVec(2);
            randlist(randIdx, 8) = tmpVec(1);
            
            tmpVec(1:3) = tmpVec(2:4);
            tmpVec(4)= setdiff(iniVec, tmpVec);
            randIdx = randIdx+1;
               
        end
        tmpVec(3) = tmpVec(4);
        tmpVec(4) = setdiff(iniVec, tmpVec);
    end
    tmpVec(2:3) = tmpVec(3:4);
    tmpVec(4) = setdiff(iniVec, tmpVec);

end

block_order = randlist(subject_id,:);

end
