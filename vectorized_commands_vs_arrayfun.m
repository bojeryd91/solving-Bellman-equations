clear, close all, clc

X = randn(1000, 1);

tic; for ii = 1:1000; X.^2; end; toc
tic; for ii = 1:1000; arrayfun(@f, X); end; toc
