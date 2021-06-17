function figurex(h)
% FIGUREX - figure() without focus stealing for existing figure. 
%
%    If figure window with handle f exists, focus will not be stolen when
%    figurex(f) is called.  Otherwise, focus is stolen (and this is
%    unavoidable).
%
%    See https://stackoverflow.com/q/8488758

fh = get(0,'Children');
for i = 1:length(fh)
    if fh(i).Number == h || strcmp(fh(i).Name,h)
        if fh(i).Number == h
            %fprintf('figurex: Found existing figure %d\n',h);
        end
        if strcmp(fh(i).Name,h)
            %fprintf('figurex: Found existing figure %s\n',h);
        end
        set(0,'CurrentFigure',h);
    else
    	figure(h);
    end
end
