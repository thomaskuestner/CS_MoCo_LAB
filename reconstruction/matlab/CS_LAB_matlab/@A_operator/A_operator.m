function A_oper = A_operator( A, At )
%A_OPERATOR Linear operator class constructor.
%   A = A_operator(A) creates an A_operator object from the function handle A.
if nargin ~= 2
    error('Syntax: A = A_operator( A, A_transpose )');
elseif isa(A,'A_operator')
    A_oper = A;
elseif isa(A, 'function_handle')
    A_oper.func = A;
    A_oper.transpose = At;
    A_oper = class(A_oper,'A_operator');
else
    error('Input must be a function_handel');
end