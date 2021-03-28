function [defval] = bt_defaultval(config,var,val)
% If the field exists, use its value. If it doesn't, set default based on
% input 'val'.

if isfield(config,var)
    defval = eval(['config.' var]);
else
    defval = val;
end
end
