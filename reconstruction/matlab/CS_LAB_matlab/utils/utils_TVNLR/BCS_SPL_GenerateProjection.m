% 
% function Phi = BCS_SPL_GenerateProjection(block_size, subrate, filename)
%
%   This function generates the random projection matrix
%   Phi for the given block size and subrate.
%
%   Phi is returned as a M x N matrix, where N = block_size *
%   block_size, and M = round(subrate * N).
%
%   If filename is not specified, the Phi matrix is generated
%   randomly and returned. If filename is specified, but the file does
%   not exist, Phi is generated randomly and written to
%   filename. Finally, if the file does exist, Phi is simply read from
%   the file.
%
%   In all cases, Phi is generated as a random N x N matrix and then
%   truncated to M rows. Phi is stored in filename as a N x N matrix
%   and truncated to M rows upon reading.
%
%   See:
%     S. Mun and J. E. Fowler, "Block Compressed Sensing of Images
%     Using Directional Transforms," submitted to the IEEE
%     International Conference on Image Processing, 2009
%
%   Originally written by SungKwang Mun, Mississippi State University
%

%
% BCS-SPL: Block Compressed Sensing - Smooth Projected Landweber
% Copyright (C) 2009-2011  James E. Fowler
% http://www.ece.mstate.edu/~fowler
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%


function Phi = BCS_SPL_GenerateProjection(block_size, subrate, filename)

N = block_size * block_size;
M = round(subrate * N);

if ((nargin == 3) && exist(filename, 'file'))
    load(filename);
else
  Phi = orth(randn(N, N))';
  %Phi = rand(M,N,'single')- 0.5;
end

if ((nargin == 3) && (~exist(filename, 'file')))
  save(filename, 'Phi');
end

Phi = Phi(1:M, :);
