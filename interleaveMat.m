function interleavedMat = interleaveMat(A, B, varargin)
%inputs: two or three same size matrices
%  Two matrices as input: A is an n-by-m matrix
%                       : B is an n-by-m matrix
%                 output: interleavedMat of A and B of size 2n-by-m
%
%  Three matrices  input: A: n-by-m matrix
%                       : B: n-by-m matrix
%                       : C = varargin{3} n-by-m matrix
%                 output: interleavedMat of A, B and C of size 3n-by-m
%
% e.g for 2 matrices as input:
% A =
%      1     2
%      3     4
% B =
%      5     6
%      7     8
% 
% mat = interleaveMat(A,B)
% mat=
%      1     2
%      5     6
%      3     4
%      7     8

nColumns = size(A,2);

switch nargin
    case 2          
         AB = [A,B]'; 
         interleavedMat = reshape(AB(:),nColumns,[])';
    case 3
         C = varargin{1};
         ABC = [A,B,C]'; 
         interleavedMat = reshape(ABC(:),nColumns,[])';
     case otehrwise
         error('Unexpected inputs for interleaveMat function.')

 end
