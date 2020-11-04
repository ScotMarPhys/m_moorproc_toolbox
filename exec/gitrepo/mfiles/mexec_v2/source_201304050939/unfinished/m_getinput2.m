function v = m_getinput(msg,type,opt)

% Call with type = 's' for string or 'd' for double
% If you respond c/r the defaults are ' ' and NaN

% optional third argument opt can take value 'no_ot' to avoid adding to
% MEXEC_A.MARGS_OT


m_common

argsin = MEXEC_A.MARGS_IN;
argsot = MEXEC_A.MARGS_OT;

% If MEXEC_A.MARGS_IN is absent or empty, get input from keyboard
% If MEXEC_A.MARGS_IN has elements, take the input from the first element of the cell array

if isempty(argsin) % MEXEC_A.MARGS_IN is empty, prompt at keyboard
    if strcmp(type,'s')
        v = input(msg,'s');
        if isempty(v); v = ' '; end
    else
        v = input(msg);
        if isempty(v); v = nan; end
    end
else  % MEXEC_A.MARGS_IN has elements, take the input from the first element of the cell array
    v = argsin{1};
    % keyboard input is always converted to character
    % do the same for MEXEC_A.MARGS_IN or varargin input
    if ~ischar(v); v = num2str(v); end
    argsin(1) = [];
    disp(msg);
    disp(v);
end


MEXEC_A.MARGS_IN = argsin;
if exist('opt','var') ~= 1
    argsot = [argsot v];
    MEXEC_A.MARGS_OT = argsot;
elseif ~strcmp(opt,'no_ot')
    argsot = [argsot v];
    MEXEC_A.MARGS_OT = argsot;
end
