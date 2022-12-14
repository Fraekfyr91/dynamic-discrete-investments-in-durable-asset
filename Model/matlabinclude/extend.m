% extend: if number of elements in matrix A is smaller than n, then last element is repeated so that B is n dimensional
% 			matrix. If n is smaller that the number of elements in A, only the first n elements are returned.
% syntax: [B]=extend(A, n);  
% example: extend({1,2},4) returns the 4x1 matrix
%
% 	[1]
% 	[2]
% 	[2]
% 	[2]
%   Alternaively if called with 3 inputs and A is two dimensional matrix, 
%	then both columns and rows are repeated as descrived above
%
%   example: 
%
%	A=[1 2 3
%     4 5 6];
%
%	B=extend(A, 3,5)
%   B =
% 		3x5 matrix
%   	1    2    3    3    3
%   	4    5    6    6    6
%   	4    5    6    6    6

% Bertel Schjerning, University of Copenhagen, April 2019

function [B]=extend(A, n, m);
	if nargin==2 	
		n1=numel(A);
		if n1>n
			B=A(1:n);
			return
		end
		B=cell(n,1);
		B(1:n1)=A;
		if n1<n
			B(n1+1:n) =repmat(A(end), n-n1,1);
		end
	end

	if nargin==3 	
		[n1, m1]= size(A);
		B=zeros(n,m);
		B(1:min(n1,n),1:min(m1,m))=A(1:min(n1,n),1:min(m1,m));
		if m1<m 
			B(:,m1+1:m) = repmat(B(:,m1), 1,m-m1);
		end
		if n1<n
			B(n1+1:n,:) =repmat(B(min(n1,n),:), n-n1,1);
		end
	end

end