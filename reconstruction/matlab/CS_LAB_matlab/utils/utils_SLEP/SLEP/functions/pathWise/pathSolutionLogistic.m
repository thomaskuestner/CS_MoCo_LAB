function [X, C]=pathSolutionLogistic(A, y, z, opts)
%
%% Fuction pathSolution:
%      Solving the pathwise solutions
%
%% Input & Output parameters
%  See the description of the related functions
%
%% Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
%
% You are suggested to first read the Manual.
%
% For any problem, please contact with Jun Liu via j.liu@asu.edu
%
% Last modified 2 August 2009.
%
% Related functions:
%  sll_opts,
%  eppVector, eppMatrix, eplb,
%
%  LogisticR, LogisticC,
%  nnLogisticR, nnLogisticC
%  glLogisticR, mtLogisticR, mclLogisticR
%%

switch(opts.fName)
    case 'LogisticR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X
        C=zeros(1,z_num);            % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=LogisticR(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution x
        C(1,z_ind(1))=c;             % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function LogisticR
            [x, c, funVal]=LogisticR(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution x
            C(1,z_ind(i))=c;         % store the solution c
        end

    case 'LogisticC'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(z);    % sort z in an ascending order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X
        C=zeros(1,z_num);            % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=LogisticC(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution x
        C(1,z_ind(1))=c;             % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function LogisticC
            [x, c, funVal]=LogisticC(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution x
            C(1,z_ind(i))=c;         % store the solution c
        end

    case 'glLogisticR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X
        C=zeros(1,z_num);            % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=glLogisticR(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution x
        C(1, z_ind(1))=c;            % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function glLogisticR
            [x, c, funVal]=glLogisticR(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution x
            C(1, z_ind(i))=c;        % store the solution c
        end

    case 'mtLogisticR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        k=length(opts.ind)-1;        % the number of tasks
        X=zeros(n,k,z_num);          % set the size of output X
        C=zeros(z_num, k);           % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=mtLogisticR(A, y, z_value(1), opts);

        X(:,:,z_ind(1))=x;           % store the solution x
        C(z_ind(1), :)=c;            % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function mtLogisticR
            [x, c, funVal]=mtLogisticR(A, y, z_value(i), opts);

            X(:,:, z_ind(i))=x;      % store the solution x
            C(z_ind(i), :)=c;        % store the solution c
        end

    case 'mcLogisticR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        k=size(y,2);                 % the number of classes (tasks)
        X=zeros(n,k,z_num);          % set the size of output X
        C=zeros(z_num, k);           % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=mcLogisticR(A, y, z_value(1), opts);

        X(:,:,z_ind(1))=x;           % store the solution x
        C(z_ind(1), :)=c;            % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function mcLogisticR
            [x, c, funVal]=mcLogisticR(A, y, z_value(i), opts);

            X(:,:, z_ind(i))=x;      % store the solution x
            C(z_ind(i), :)=c;        % store the solution c
        end

    case 'nnLogisticR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X
        C=zeros(1,z_num);            % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=nnLogisticR(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution x
        C(1,z_ind(1))=c;             % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function LogisticR
            [x, c, funVal]=nnLogisticR(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution x
            C(1,z_ind(i))=c;         % store the solution c
        end

    case 'nnLogisticC'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(z);    % sort z in an ascending order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X
        C=zeros(1,z_num);            % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=nnLogisticC(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution x
        C(1,z_ind(1))=c;             % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function LogisticC
            [x, c, funVal]=nnLogisticC(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution x
            C(1,z_ind(i))=c;         % store the solution c
        end

    case 'mtLogisticC'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(z);    % sort z in an ascending order
        n=size(A,2);                 % the dimensionality of the data
        k=length(opts.ind)-1;        % the number of tasks
        X=zeros(n,k,z_num);          % set the size of output X
        C=zeros(z_num, k);           % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=mtLogisticC(A, y, z_value(1), opts);

        X(:,:,z_ind(1))=x;           % store the solution x
        C(z_ind(1), :)=c;            % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function mtLogisticC
            [x, c, funVal]=mtLogisticC(A, y, z_value(i), opts);

            X(:,:, z_ind(i))=x;      % store the solution x
            C(z_ind(i), :)=c;        % store the solution c
        end

    case 'mcLogisticC'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(z);    % sort z in an ascending order
        n=size(A,2);                 % the dimensionality of the data
        k=size(y,2);                 % the number of classes (tasks)
        X=zeros(n,k,z_num);          % set the size of output X
        C=zeros(z_num, k);           % set the size of output C

        % run the code to compute the first solution
        [x, c, funVal]=mcLogisticC(A, y, z_value(1), opts);

        X(:,:,z_ind(1))=x;           % store the solution x
        C(z_ind(1), :)=c;            % store the solution c

        % set .init for warm start
        opts.init=1;                 % using the defined ones

        for i=2:z_num
            opts.x0=x;               % warm-start of x
            opts.c0=c;               % warm-start of c

            % run the function mcLogisticC
            [x, c, funVal]=mcLogisticC(A, y, z_value(i), opts);

            X(:,:, z_ind(i))=x;      % store the solution x
            C(z_ind(i), :)=c;        % store the solution c
        end

    otherwise
        fprintf('\n The function value specified in opts.fName is not supported!');
end