function a = concat(a, varargin)
% A = CONCAT(A, STRUCT*)
% concatenates two or more structs on a per field basis
%
% Use e.g. to combine fixmats from different subjects
%
% Example:
%     foo =
%     a: [1 2 3]
%     bar =
%     a: [4 5 6]
%     concat(foo, bar)
%     a: [1 2 3 4 5 6]
    fn = fieldnames(a);
    for i=(1:nargin-1)
        fn = intersect(fn, fieldnames(varargin{i}));
    end
    for i=(1:length(fn))
        for j=(1:nargin-1)
            a.(fn{i}) = [a.(fn{i}), varargin{j}.(fn{i})];
        end;
    end
end