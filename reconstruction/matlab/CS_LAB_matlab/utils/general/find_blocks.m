function blocks = find_blocks( search_vector )
% find connected blocks and print beginning to end in blocks
% blocks(1,:) contains the beginnings of the blocks
% blocks(2,:) contains the endings of the blocks
% search vector can contain 0 and 1 or any number represantation 
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

% ensure numerical array -> at the moment just 0 and 1 allowed in search
% vector
if(isnumeric(search_vector))
    l_helper = false(1,max(search_vector));
    l_helper(search_vector) = true;
    search_vector = double(l_helper);
else
    search_vector = double(search_vector);
end

blocks = [];
% determine block edges
block_edges = conv(search_vector, [1 -1]);
ind_minus = find(block_edges == -1);
block_edges(ind_minus -1) = 1;
block_edges(ind_minus) = 0;
block_edges = block_edges(1:size(search_vector,2));
ind_block_edges = find(block_edges);

% determine blocks
k = 1;
while(~isempty(ind_block_edges))
    blocks = [blocks, [find(block_edges == 1,1,'first');0]];
    if(any(search_vector(blocks(1,k):end) == 0))
        blocks(2,k) = find(search_vector(blocks(1,k):end) == 0,1,'first')-1 + blocks(1,k)-1;
    else
        % entire last block till end of search_vector
        blocks(2,k) = ind_block_edges(end);
    end
    
% %     if(~all(search_vector(blocks(1,k):blocks(2,k))))
% %         % block is not completely connected
% %         % -> find first occuring zero after start
% %         blocks(2,k) = find(search_vector(blocks(1,k):end) == 0,1,'first')-1;
    ind_block_edges(ismember(ind_block_edges,blocks(:))) = [];
    block_edges(blocks(:,k)) = 0;
    k = k+1;
end


end