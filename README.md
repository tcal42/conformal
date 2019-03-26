# conformal

Contributors: Trevor Caldwell, Anne Greenbaum, Kenan Li (University of Washington), and Nick Trefethen (University of Oxford)

This repo contains a Chebfun/Matlab implementation of a conformal mapping method based on the Kerzman-Stein integral equation. The code is intended for mapping regions with smooth boundaries to the unit disk, making use of the accuracy afforded by Chebfun representations. This was originally written by me (Trevor Caldwell) together with Kenan Li and our advisor Anne Greenbaum. Nick Trefethen gave helpful modificiations during Householder XX in 2017, and we also make use of their AAA algorithm for rational approximation in the complex plane (see the Chebfun repository). 

Included in this repo are basic versions of the boundary mapping (kerzstein.m and bcfun_eval.m) and the conformal map derived from the boundary map (conformal.m). Also included are some functions which make use of this conformal code, namely phi_eval.m, which outputs a matrix phi(A), where phi is the conformal map sending the numerical range of A to the unit disk, and optimal_blasckhe.m, which finds the optimal Blaschke product for A by maximizing the 2-norm of B(phi(A)), where B is a Blaschke product of degree length(A) - 1. These codes can be used for studying the behavior of functions of matrices, in addition to computing conformal maps. 
