# Introduction

This project implements the Rolling Guidance Normal Filter method.

# Thirdparty libararies

- Eigen

# Current Problem

The laplacian implementation is
$\Delta v*i = \frac{1}{2A_i} \sum*{j \in N(i)} (\cot \alpha*{ij} + \cot \beta*{ij}) (v_j - vi). $

# Reference

Rolling Guidance Normal Filter for Geometric Processing, Peng-Shuai Wang, Xiao-Ming Fu, Yang Liu, Xin Tong, Shi-Lin Liu and Baining Guo, ACM Transactions on Graphics (SIGGRAPH Asia), 34(6), 2015
Fast High-dimensional Filtering Using the Permutohedral Lattice, Andrew Adams, Jongmin Baek, and Abe Davis, Eurographics 2010
