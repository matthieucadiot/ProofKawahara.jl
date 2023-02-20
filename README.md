# Computer-assisted proofs of solitons, eigencouples and orbital stability in the Kawahara equation.


## Introduction

This Julia code is a complement to the article 
[1] : "Rigorous computation of solutions of semi-linear PDEs on unbounded domains via spectral methods" (M. Cadiot, J-P. Lessard, J-C. Nave) as it provides the necessary rigorous computations that are needed in Section 6.


## The Kawahara equation

The Kawahara equation
 ![equation](http://latex.codecogs.com/gif.latex?%5Clambda_2u%27%27%27%27%20%2B%20%5Clambda_1u%27%27%20%2B%20u%20%2B%20%5Clambda_3u%5E2%20%3D%200)
is known to have solutions on <img src="https://latex.codecogs.com/gif.latex?\mathbb{R}" /> that decay to zero at infinity. These solutions are called solitary waves or soliton (see [1] for a complete description).

The present code provides the rigorous numerics for the proof of solitons of the Kawahara equation using the analysis of [1] (specifically the Section 6). In the beginning of the main code (main_proof_Kawahara), the user can choose the values for N, d, T and c, that are described in [1]. In particular, T and c need to be chosen such that
<img src="https://latex.codecogs.com/gif.latex?0 \leq  T < 0.397" /> 
<img src="https://latex.codecogs.com/gif.latex?c < 1- \frac{a(T)^2}{4b(T)^2}" />
where
<img src="https://latex.codecogs.com/gif.latex?a(T) = \frac{1-3T}{6}" /> 
<img src="https://latex.codecogs.com/gif.latex?b(T) = \frac{19 - 30T - 45T^2}{360}" />.


