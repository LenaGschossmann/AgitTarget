  function splited = split(fixmat, varargin)
  % splited = split(fixmat, ['field', values]*)
        if (nargin < 1)
            error('splited takes at least one parameter');
        end
        if nargin==1
            splited = fixmat;
            return
        end
        splited = [];
        
        function add_splited(value)
            l = size(splited, 2);
            if l==0
                splited = value;
            else
                splited(l+1:l+size(value,2)) = value;
            end
        end
    
        function ssf = single_split( elems, ffield, frange )
            if isempty(frange)
                frange = unique([elems.(ffield)]);
            end
            for r = (1:length(frange))
                ind = find( [elems.(ffield)] == frange(r) );
                for fn = fieldnames(elems)'
                    ssf(r).(fn{:}) = elems.(fn{:})(:,ind);
                end
            end
        end
   
        for c = fixmat
            add_splited( split( single_split( c, varargin{1}, varargin{2} ), varargin{3:end} ) );
        end
  end