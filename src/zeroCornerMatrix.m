function [ output ] = zeroCornerMatrix( input )
%ZEROCORNERMATRIX Set the four corner values of a matrix to zero.
%
%   output = zeroCornerMatrix(input)
%
% INPUTS:
%   input  : [ny x nx] array
%
% OUTPUT:
%   output : input with output(1,1), output(1,end), output(end,1) and
%            output(end,end) set to zero
%
% Brzenski, Jared
% 01/31/2025
% -------------------------------------------------------------------------
    output = input;
    output(1,1)     = 0;
    output(1,end)   = 0;
    output(end,1)   = 0;
    output(end,end) = 0;
end
