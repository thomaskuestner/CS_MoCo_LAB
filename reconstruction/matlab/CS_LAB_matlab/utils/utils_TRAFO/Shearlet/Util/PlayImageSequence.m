%% Play an image sequence
%% The input 3-D array must be of type uint8

function PlayImageSequence(varargin)

figure
hold on;

if nargin == 1
    X = varargin{1};
    for n = 1 : size(X, 3)
        imshow(X(:,:,n));
        pause(0.03);
    end
else
    R = varargin{1};
    G = varargin{2};
    B = varargin{3};
    for n = 1 : size(R, 3)
        imshow(cat(3, R(:,:,n), G(:,:,n), B(:,:,n)));
        pause(0.03);
    end
end

hold off;