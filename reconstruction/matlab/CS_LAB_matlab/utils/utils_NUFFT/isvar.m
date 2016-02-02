 function tf = isvar(name)
%function tf = isvar(name)
%	determine if "name" is a variable in the caller's workspace

if nargin < 1
	help isvar
	error arg
end

tf = logical(1);
evalin('caller', [name ';'], 'tf=logical(0);')
