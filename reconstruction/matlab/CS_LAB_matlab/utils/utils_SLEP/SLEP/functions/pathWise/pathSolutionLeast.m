function X=pathSolutionLeast(A, y, z, opts)
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
%  sll_opts, initFactor,
%  eppVector, eppMatrix, eplb,
%
%  LeastR, LeastC,
%  nnLeastR, nnLeastC
%  glLeastR,  mtLeastR, mcLeastR
%%

switch(opts.fName)
    case 'LeastR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=LeastR(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution

        % set .init for warm start
        opts.init=0;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function LeastR
            [x, funVal]=LeastR(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution
        end

    case 'LeastC'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(z);    % sort z in an ascending order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=LeastC(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution

        % set .init for warm start
        opts.init=0;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function LeastC
            [x, funVal]=LeastC(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution
        end

    case 'glLeastR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=glLeastR(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution

        % set .init for warm start
        opts.init=0;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function glLeastR
            [x, funVal]=glLeastR(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution
        end

    case 'mtLeastR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        k=length(opts.ind)-1;        % the number of tasks
        X=zeros(n,k,z_num);          % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=mtLeastR(A, y, z_value(1), opts);

        X(:,:,z_ind(1))=x;           % store the solution

        % set .init for warm start
        opts.init=0;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function mtLeastR
            [x, funVal]=mtLeastR(A, y, z_value(i), opts);

            X(:,:, z_ind(i))=x;      % store the solution
        end

    case 'mcLeastR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        k=size(y,2);                 % the number of tasks
        X=zeros(n,k,z_num);          % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=mcLeastR(A, y, z_value(1), opts);

        X(:,:,z_ind(1))=x;           % store the solution

        % set .init for warm start
        opts.init=1;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function mcLeastR
            [x, funVal]=mcLeastR(A, y, z_value(i), opts);

            X(:,:, z_ind(i))=x;      % store the solution
        end

    case 'nnLeastR'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(-z);   % sort z in a decresing order
        z_value=-z_value;            % z_value in a decreasing order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=nnLeastR(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution

        % set .init for warm start
        opts.init=0;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function LeastR
            [x, funVal]=nnLeastR(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution
        end

    case 'nnLeastC'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(z);    % sort z in an ascending order
        n=size(A,2);                 % the dimensionality of the data
        X=zeros(n,z_num);            % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=nnLeastC(A, y, z_value(1), opts);

        X(:,z_ind(1))=x;             % store the solution

        % set .init for warm start
        opts.init=0;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function LeastC
            [x, funVal]=nnLeastC(A, y, z_value(i), opts);

            X(:,z_ind(i))=x;         % store the solution
        end

    case 'mtLeastC'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(z);    % sort z in an ascending order
        n=size(A,2);                 % the dimensionality of the data
        k=length(opts.ind)-1;        % the number of tasks
        X=zeros(n,k,z_num);          % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=mtLeastC(A, y, z_value(1), opts);

        X(:,:,z_ind(1))=x;           % store the solution

        % set .init for warm start
        opts.init=0;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function mtLeastC
            [x, funVal]=mtLeastC(A, y, z_value(i), opts);

            X(:,:, z_ind(i))=x;      % store the solution
        end

    case 'mcLeastC'
        z_num=length(z);             % the number of parameters
        [z_value, z_ind]=sort(z);    % sort z in an ascending order
        n=size(A,2);                 % the dimensionality of the data
        k=size(y,2);                 % the number of tasks
        X=zeros(n,k,z_num);          % set the size of output X

        % run the code to compute the first solution
        [x, funVal]=mcLeastC(A, y, z_value(1), opts);

        X(:,:,z_ind(1))=x;           % store the solution

        % set .init for warm start
        opts.init=1;                 % using .initFactor

        for i=2:z_num
            opts.x0=x;               % warm-start

            % run the function mcLeastC
            [x, funVal]=mcLeastC(A, y, z_value(i), opts);

            X(:,:, z_ind(i))=x;      % store the solution
        end

    otherwise
        fprintf('\n The function value specified in opts.fName is not supported!');
end