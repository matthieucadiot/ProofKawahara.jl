# Computer-assisted proofs of solitons, eigencouples and orbital stability in the Kawahara equation.



Table of contents:


* [Introduction](#introduction)
* [The Kawahara equation](#kawahara)
   * [Proof of solitons](#solitons)
   * [Proof of the first 3 eigencouples](#eigencouples)
   * [Proof of orbital stability](#stability)
* [Utilisation](#utilisation)
* [References](#references)
* [License](#license)



# Introduction

This Julia code is a complement to the article 

[1](#eq1) : "Rigorous computation of solutions of semi-linear PDEs on unbounded domains via spectral methods" (M. Cadiot, J-P. Lessard, J-C. Nave)

as it provides the necessary rigorous computations that are needed in Section 6.


# The Kawahara equation

The Kawahara equation
$$\lambda_2u'''' + \lambda_1u'' + u + \lambda_3u^2 = 0$$
is known to have solutions on $\mathbb{R}$ that decay to zero at infinity. These solutions are called solitary waves or soliton (see [1] for a complete description).

## Proof of solitons

The present code provides the rigorous numerics for the proof of solitons of the Kawahara equation using the analysis of [1] (specifically the Section 6). In the beginning of the main code (main_proof_Kawahara), the user can choose the values for N, d, T and c, that are described in [1]. In particular, T and c need to be chosen such that
 - $0 \leq T < 0.397$ 
 - $c < 1- \frac{a(T)^2}{4b(T)}$    
where
- $a(T) = \frac{1-3T}{6}$
- $b(T) = \frac{19 - 30T - 45T^2}{360}$.   

The code will compute rigorously the needed bounds of Section 6 of [1] and validate of not the computer-assisted proof. If the computer-assisted proof succeeds, the radius for the smallest and biggest ball of contraction is displayed.

## Proof of the first 3 eigencouples

If the proof of the soliton is achieved, the code will then compute approximations for the first 3 eigencouples of the linearization around the proved soliton. Then, the needed bounds for the proof of eigencouples are computed, following the analysis of Section 6 of [1]. In particular, the proof will be valid if T was chosen such that
 $\frac{1}{3} < T < 0.397$. 
  
   

 If the computer-assisted proof succeeds, the radius for the smallest and biggest ball of contraction is displayed for each eigencouple.
 
 
 ## Proof of orbital stability

If the proof of the first 3 eigencouples is achieved, the code will then try to prove Theorem 6.20 in [1]. In particular, we want to prove that the 3 eigenvalues obtained beforehead, are actually the first 3 ones. The algorithm is explained in the proof of Theorem 6.20 and uses Lemma 6.19. 

 If the computer-assisted proof of Theorem 6.20 succeeds, the value for  <img src="https://latex.codecogs.com/gif.latex?\tau" /> (Proposition 6.23) is computed rigorously. In particular, we check that <img src="https://latex.codecogs.com/gif.latex?\tau" /> is striclty negative, and if that is the case, then we obtain that the proved soliton is orbitally stable.
 
 
 # Utilisation
 
 # License
  The theme is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
