function v = get_option(options, op, opd)

% A utility function to fill a struct

if isfield(options, op)
    if isempty(eval(['options.' op]))
        v = opd;
    else
        v = eval(['options.' op]);
    end
else
    v = opd;
end