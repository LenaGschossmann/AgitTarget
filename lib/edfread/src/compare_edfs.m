function x = compare_edfs(data1, data2,tabs)
% function x = compare_edfs(data1, data2,tabs)
%
% Recursuvely compare two structs and check if all fields contain the
% same data (and if fields are the same and so on).
%
% Input:
%   data1  =  Struct number one to compare
%   data2  =  Struct number two for the comparison.
% Output:
%   x      =  The number of comparisons carried out.
%
%   Niklas Wilming (nwilming@uos.de, June 2009)
%

warning('on', 'MATLAB:concatenation:integerInteraction');
x = 0;
t = [];
for o=1:tabs
    t = [t '   '];
end
if isstruct(data1) && isstruct(data2)
    fields1 = fieldnames(data1);
    fields2 = fieldnames(data2);
    
    if length(fields1) ~= length(fields2)
        error('data 1 and 2 have a different number of data fields');
    end
    for f = 1:length(fields1)
        display([t 'Comparing fields: ' fields1{f} ' / ' fields2{f}])
        if strcmp(fields2{f},'msg')
            x = x + compare_msg_fields(data1, data2);            
        else
            warning('off','MATLAB:concatenation:integerInteraction');
            x = x + compare_edfs([data1.(fields1{f})], [data2.(fields2{f})],tabs+1);
            warning('on', 'MATLAB:concatenation:integerInteraction');
        end
    end
elseif ischar(data1) && ischar(data2)
    if data1~=data2
        x = x + 1;
        error([t 'Difference in data1 and data2 (both are chars)']);
    end
    
elseif isnumeric(data1) && isnumeric(data2)
    if length(data1) ~= length(data2)
        error('t DIFF_DATA_FIELDS have  different length');
    end
    for k = 1:length(data1)
        x = x +1;
        if data1(k) ~= data2(k)
            if ~isnan(data1(k)) && ~isnan(data2(k))
                display([t 'Data fields have different entries: ' num2str(data1(k)) ' / ' num2str(data2(k))])
            end
        end
    end
elseif islogical(data1) && islogical(data2)
    x = x +1;
    if data1 ~= data2
        error('Logicals are not equal')
    end
else
    [data1]
    [data2]
    error('Niklas has forgotten to include a type comparison. Go tell him.')
end

    function  x = compare_msg_fields(m1, m2)
        x = 0;
        l1 = length(m1);
        l2 = length(m2);
        if l1 ~= l2
            display([t 'Message fields have different length']);
        end
        for kk = 1:l1
            x = x+1;
            if ~strcmp(m1(kk).msg,m2(kk).msg)
                display([t 'Messages are different. Index: ' num2str(kk)])
            end
        end
    end

end