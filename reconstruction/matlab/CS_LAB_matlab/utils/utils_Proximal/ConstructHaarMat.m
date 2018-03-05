function HaarTransformationMatrix=ConstructHaarMat(WidthOfSquareMatrix)
%--------------------------------------------------------------------------
%Create Haarwavelet transformation matrix H for the matrix vector
%mulplication implimentation of Haar wavelet transformation.
%This function uses the following nice formula to create the Haar
%transformation matrix:
%               H_n=1/sqrt(2)[H_(n/2) kron (1 1)
%                             I_(n/2) kron (1 -1)],
%                              where 'kron' denotes the kronecker product.
%The iteration starts with H_1=[1]. The normalization constant 1/sqrt(2)
%ensure that H_n^T*H_n=I, where I is identity matrix. Haar wavelets are the
%rows of H_n.
%--------------------------------------------------------------------------
%function HaarTransformationMatrix=ConstructHaarWaveletTransformationMatrix(WidthOfSquareMatrix)
% Input:
%       WidthOfSquareMatrix: the width of sqaure Haar wavelet
%                                  transformation matrix. 
% Output:
%       HaarTransformationMatrix: Ceated Haar transformation matrix and
%       it's size is the power of 2,i.e., 2, 4,
%       8,16,32,64,128,256,512,1024,2048,4096,etc.
%--------------------------------------------------------------------------
%Author: Jin Qi
%Email: jqichina@hotmail.com
%Date:  11/4/2011
%--------------------------------------------------------------------------
%Test function for one dimentional signal
%
% n=8;
% S=rand(n,1); % for one dimentional signal
% H=ConstructHaarWaveletTransformationMatrix(n);
% C=H*S; %Haar wavelet transformation
% RS=H'*C;% Reconstruct signal
% fprintf('reconstruction error is %f',norm(S-RS));%plot reconstruction error
%--------------------------------------------------------------------------
%Test for two dimentional signal,ie. image
%
% n=8;
% S=rand(n,n); % for one dimentional signal
% H=ConstructHaarWaveletTransformationMatrix(n);
% C=H*S*H'; %Haar wavelet transformation
% RS=H'*C*H;% Reconstruct signal
% fprintf('reconstruction error is %f',norm(S-RS,'fro'));%plot reconstruction error

n=WidthOfSquareMatrix; % copy the parameter

% check input parameter and make sure it's the power of 2
Level=log2(n);
if 2^Level<n, error('please ensure the value of input parameter is the power of 2');end 

%Initialization
H=[1];
NC=1/sqrt(2);%normalization constant
LP=[1 1]; 
HP=[1 -1];

% iteration from H=[1] 
for i=1:Level
    H=NC*[kron(H,LP);kron(eye(size(H)),HP)];
end
HaarTransformationMatrix=H;
