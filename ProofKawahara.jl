#Computer assisted proof of solitons and periodic solutions for the Kawahara equation : u + λ1u'' + λ2u'''' + λ3u² = 0
# The following code computes rigorously the needed constants of Section 6 of
# "Rigorous computation of solutions of semi-linear parital differential equations on unbounded domains via spectral methods"  M. Cadiot, J-P. Lessard, J-C. Nave

# In a first time, the code computes the constants of Theorem 5.6 using the analysis of Section 6. From this we can check if the proof of the soliton and the 
# associated periodic function is verified or not (see Remark 5.7 for the periodic part). We essentially prove Theorem 6.13.

# Then we want to apply Theorem 5.11 to prove the first 3 eigencouples of DF(u*) as detialed in Section 6.3. Computing the needed constansts, we prove
# Theorem 6.17 and obtain that the eigencouples are simple.

# At this point, we can prove Theorem 6.20 to obtain that the previously proved eigencouples are actually the first free eigenvalues !
# (using the function no_eigen_outside)

# Finally we compute a value for τ in Proposition 6.23 and prove that τ <0, which yields that the soliton is (orbitally) stable (using Lemma 6.22).

#####################################################################################################################################################################

# Choice of the parameters for the proof of solitons:
#  - 0 < T < 0.397
#  - c < 1- a(T)^2/(4b(T)) where a(T) = (1-3T)/6  and b(T) = (19 - 30T - 45T^2)/360

# Choice of the parameters for the proof the stability of the solitons:
#  - 1/3 < T < 0.397
#  - c < 1- a(T)^2/(4b(T)) where a(T) = (1-3T)/6  and b(T) = (19 - 30T - 45T^2)/360


# Needed packages
using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, FFTW, PrettyTables


#################################### List of the needed functions : go directly to line 355 for the main code ################################################# 

function newton_kawahara(c, T, N, d)

    fourier = CosFourier(N,π/d)
    
    X= 2*d/(2*N+1).*(-N:N)
    u = 2*(c-1)./((cosh.(sqrt(abs(3*(c-1)/(2*(1-3*T))))*X)).^2)
    V = real(FFTW.fftshift(FFTW.fft(FFTW.ifftshift(u))/(2*N+1)))
    V = vec(V[N+1:2*N+1])
    U = Sequence(fourier, vec(real(V)))

    D² = project(Derivative(2), fourier, fourier,Float64)
    D⁴ = project(Derivative(4), fourier, fourier,Float64)
    β = (19 - 30T - 45T^2)/(360*(1-c)) ;  α = (1-3T)/(6*(1-c)) ; γ = 3/(4*(1-c))

    L = I +α*D² + β*D⁴ # petite modification pour minimiser le nombre d'opération
    G = project(Multiplication(γ*U),fourier, fourier,Float64)
    F = L*U +G*U

    n_F = norm(F,2)
    n_max = 20 ; p =0;
    err = 1.0e-16
    while (n_F > err)&&(p <n_max)
        DF = L + 2*G
        U = U - DF\F
        G = project(Multiplication(γ*U),fourier, fourier,Float64)

        F = L*U + G*U
        n_F = norm(F,2) 
        p = p + 1
    end
    if p == n_max
        display("Maximum number of iterations reached")
        return p
    else
        if norm(U) < 1.0e-8
            display("Newton might have converged to the zero solution")
            return 0
        else
            display("Newton converged to a non trivial solution")
            return U
        end
    end
        
    
end



function trace(N,d,id)
    if id==0
        w = interval.((0:N));
        Y = (interval((pi))^2/interval((d))^2)*w.^2;
        X = interval.(((2*(-1).^(0:N)))); X[1] = interval((1));

        f = [X';
            (X.*Y)']
        return f
    else
        w = interval.(big.(0:N));
        Y = (interval(big(pi))^2/interval(big(d))^2)*w.^2;
        X = interval.(big.((2*(-1).^(0:N)))); X[1] = interval(big(1));

        f = [X';
            (X.*Y)']
        return f 
    end
end


function trace_full(N,d)
    w = -N:N;
    p = interval(pi)
    d = interval(d)
    Y2 = (p/d)*interval.(w);
    Y1 = (p^2/d^2)*interval.(w).^2;
    Y3 = (p^3/d^3)*interval.(w).^3;
    X = interval.((-1).^w);

    f = [X';
        (X.*Y1)';
         (X.*Y2)';
         (X.*Y3)']
    return f
end





function proof_eigen_kawahara(U,U0,α,β,γ,N,d,rmin)

    pspace = space(U)
    fourier = space(U0)

    V = component(U,2)
    λ = real(component(U,1)[1])

    D² = project(Derivative(2), fourier, fourier,Complex{Interval{Float64}})
    D⁴ = project(Derivative(4), fourier, fourier,Complex{Interval{Float64}})
    L = 1/(1-λ)*interval.(real.((I + α*D² +β*D⁴)))
    Lᵢ = diagm(ones(2*N+1)./diag(coefficients(L)))
    S = trace_full(N,d); Sᵀ = S' ;

    V = V - Sequence(fourier,Lᵢ*Sᵀ*inv(S*Lᵢ*Sᵀ)*S*coefficients(V))

    αk = α/(1-λ) ; βk = β/(1-λ) ; γk = γ/(1-λ)

    U = Sequence(pspace,vec([λ; coefficients(V)]))
    # display("projection on X⁴₀ created")

    V = 1/norm(V,1)*V
    x = -αk/(2*βk)
    y = sqrt(abs(αk^2-4*βk))/(4*βk)

    a = sqrt( (sqrt(x^2+4*y^2) - x)/2 )
    b = y/a

    
    # Computation of κ defined in Proposition 6.1
    x = sqrt( abs((sqrt(αk^2+12βk) - αk)/(6βk)) )
    Lλ = L - λ/(1-λ)*I;



    Q = project(Multiplication(2*γk*U0),fourier, fourier,Complex{Interval{Float64}})

    T = LinearOperator(pspace,pspace, [[0   -(sqrt(2*d)*coefficients(Lλ*V))'] ; [-sqrt(2*d)*coefficients(Lλ*V) coefficients(Lλ+Q)]] )
    D = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) coefficients(Lλ)]]
    Di = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) Diagonal(ones(2*N+1)./diag(coefficients(Lλ)))]]
    D = LinearOperator(pspace,pspace,D)
    Di = LinearOperator(pspace,pspace,Di)
   

# # # # ################ Z BOUND ######################################################

    A = inv(mid.(T)) ; A = interval.(real.(A)) + im*interval.(imag.(A))
    norm_A = maximum([1 sup(opnorm(LinearOperator(coefficients(D*A*Di)),2))]);

    fourierE = Fourier(2*N, π/d)

    E = Sequence(fourierE ,interval.(zeros(4*N+1)))
    for n = [-2N:-1 1:2N]  
        E[n] = real((-1)^n*( (1-exp(-4*a*d))/(2*a-im*n*π/d) ))*(1-exp(-4*a*d))
    end
    E[0] = 1/(2*a)*(1-exp(-4*a*d)); E = E/(2*d) ;
    
    V0 = 2*γk*U0
    C0 = (a+b)/(4*abs(λ2)*a*b*(a^2+b^2))
    Cd = 4*d + 2*exp(-a*d)/(a*(1-exp(-3*a*d/2))) + 2/(a*(1-exp(-2*a*d)))
    
    norm_dot = maximum([0 sum(coefficients((V0).*(project(V0*E,fourier))))])
    Z11 = sqrt(norm_A^2*C0^2*2*d*norm_dot*( 2/a + exp(-2*a*d)*Cd))

    
    Q_full = project(Multiplication(2*γk*U0),fourier, Fourier(2*N, π/d),Complex{Interval{Float64}})
    Lᵢ = LinearOperator(fourier,fourier,Lᵢ)
    
    Z1 = opnorm(LinearOperator(coefficients(D*(I  - A*T)*Di)) ,2) + abs((1+norm_A)*(norm(2*γk*U0,1))/coefficients(Lλ)[2*N+1,2*N+1]) +  opnorm(LinearOperator(coefficients((Q_full-Q)*Lᵢ)),2)
    
    

    κ = sqrt(3/(8*sqrt(interval(2))*βk^(1/4)))
    r0 = rmin
    r1 = 2*γk*κ*norm(V,1)*norm_A*r0
    r2 = 2*γk*κ*r0*norm_A

    Z = sup(Z11 + Z1 + r2)
    

################# Y BOUND ######################################################
  
    A2 = component(A,2,2)
    Y = sqrt(2*interval(d))*(r1 + sqrt(2*norm(L*A2*(Lλ*V+2*γk*U0*V),2)^2 + 2*norm(2*γk*U0*V - project(2*γk*U0*V,fourier))^2)) + 2*sqrt(2*interval(d))*norm(Lλ*V+2*γk*U0*V,2)*norm(Lλ*V,2)
    
################## Z2 BOUND ######################################################
    C2 = 1/(1-λ)
    Z2 = norm_A*C2

#
################## Verification ###################################################
# # # r_star = 1.0;
# # # interval_proof = interval_of_existence(Y,Interval(0),Z2,r_star,C²Condition())
# # # display(interval_proof)
# # #
x = 1;
if 1- Z>0
    if 1-sup(4*Y*Z2) > 0
      rmin=sup((1-Z - sqrt((1-Z)^2-4*Y*Z2))/(2*Z2));
      rmax=inf((1-Z + sqrt((1-Z)^2-4*Y*Z2))/(2*Z2));
      if rmin<rmax
        display("The computer-assisted proof was successful for the eigenvalue ν with value")
        display(Float64.(mid(λ)))
        # display("Value of the minimal radius for the contraction")
        # display(Float64(rmin))
        # display("Value of the maximal radius for the contraction")
        # display(Float64(rmax))
        
      else
        display("rmin>=rmax")
        x=0
      end
    else
      display("failure: discriminant is negative")
      x=0
    end
  else
      display("failure: 1-z > 0")
      x=0
  end

  return [x==1,rmax,rmin]

end




function no_eigen_outside(r1,r2,r3,U1,U2,U3,U0,Z11,r0,α,β,γ,N)

    fourier = space(U0)


    λ1 = interval(component(real(U1),1)[1])
    λ2 = interval(component(real(U2),1)[1])
    λ3 = interval(component(real(U3),1)[1])

    D² = project(Derivative(2), fourier, fourier,Complex{Interval{Float64}})
    D⁴ = project(Derivative(4), fourier, fourier,Complex{Interval{Float64}})

    L = interval.(real.((I + α*D² +β*D⁴)))
    Q = real.(project(Multiplication(2*γ*U0),fourier, fourier,Complex{Interval{Float64}}))
    Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(L))))

    λ = interval(inf(λ1+r1))
    P = 1
    k=1

    Q_full = project(Multiplication(2*γ*U0),fourier, Fourier(2*N, π/d),Complex{Interval{Float64}})
    nQ = opnorm(LinearOperator(coefficients((Q_full-Q)*Li)),2) 

    n0 = norm(2*γ*U0,1)

    D = interval(sqrt(2))*interval.(Diagonal(vec(ones(2*N+1,1))))
    D[1,1] = interval(1)
    Di = Diagonal(ones(2*N+1)./diag(D))
    D = LinearOperator(fourier,fourier,D)
    Di = LinearOperator(fourier,fourier,Di)


    while  (P==1)&&(inf(λ)<sup(λ2-r2))
        Lλ = 1/(1-λ)*(L - interval(λ)*I)
        Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(Lλ))))

        Df = Lλ + 1/(1-λ)*Q
        A = inv(mid.(Df))
        norm_A = maximum([1 opnorm(LinearOperator(coefficients(D*Lλ*A*Di)),2)]);       # Computation for an upper bound of the norm of A

        Z = opnorm(Lλ*(I  - A*Df)*Li,2)+   abs((1+norm_A)*(n0)/coefficients(Lλ)[2*N+1,2*N+1]) + nQ + Z11*norm_A + r0
        if sup(Z) >= 1
            P=0
            return P==1
        end

        C = maximum([1,norm_A/(1-Z)])
        λ = λ + inf(1/C)

    end

    λ = interval(inf(λ2+r2))
    k=2
    
    while  (P==1)&&(inf(λ)<sup(λ3-r3))
        
        Lλ = 1/(1-λ)*(L - interval(λ)*I)
        Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(Lλ))))

        Df = Lλ + 1/(1-λ)*Q
        A = inv(mid.(Df))
        norm_A = maximum([1 opnorm(LinearOperator(coefficients(D*Lλ*A*Di)),2)]);      # Computation for an upper bound of the norm of A

        Z = opnorm(Lλ*(I  - A*Df)*Li,2)+   abs((1+norm_A)*(n0)/coefficients(Lλ)[2*N+1,2*N+1]) + nQ + Z11*norm_A + r0
        if sup(Z) >= 1
            P=0
            return P==1
        end
        
        
        C = maximum([1,norm_A/(1-Z)])
        λ = λ + inf(1/C)
    end

    λ = interval(inf(λ1-r1))
    nw = norm(2*γ*L*U0,1)
    k=3
    while  (P==1)&&(1-sup(λ)<sup(nw+1/2))
        
        Lλ = 1/(1-λ)*(L - interval(λ)*I)
        Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(Lλ))))

        Df = Lλ + 1/(1-λ)*Q
        A = inv(mid.(Df))
        norm_A = maximum([1 opnorm(LinearOperator(coefficients(D*Lλ*A*Di)),2)]);       # Computation for an upper bound of the norm of A

        Z = opnorm(Lλ*(I  - A*Df)*Li,2)+   abs((1+norm_A)*(n0)/coefficients(Lλ)[2*N+1,2*N+1]) + nQ + Z11*norm_A + r0
        if sup(Z) >= 1
            P=0
            return P==1
        end
        
        

        C = maximum([2,norm_A/(1-Z)])
        λ = λ - inf(1/C)
    end

    return P==1


end



################### PROOF OF SOLITON AND PERIODIC SOLUTIONS : MAIN CODE #################################################################################################################################################

N = 250              # number of Fourier modes : from 0 to N for a cosine series 
d= 50 ; d_big = interval(big(d)) ; d = interval(d)          # size of the domain = half period of the functions
T = 0.35 ; T_big = interval(big(T)) ; T = interval(T)  # Bond number
c = 0.9   ; c_big = interval(big(c)) ; c = interval(c) # velocity of the soliton
fourier = CosFourier(N, π/d)   # definition of the sequence space : cosine series of frequency π/d
fourier0 = CosFourier(N, π/mid(d))
fourier_big = CosFourier(N, π/d_big)
λ2 = (19 - 30T - 45T^2)/(360*(1-c)) ;  λ1 = (1-3T)/(6*(1-c)) ; λ3 = 3/(4*(1-c))  # definition of the parameters of the equation
λ2_big = (19 - 30T_big - 45T_big^2)/(360*(1-c_big)) ;  λ1_big = (1-3T_big)/(6*(1-c_big)) ; λ3_big = 3/(4*(1-c_big))

# Computation of the approximated solution via Newton's method

U0 = newton_kawahara(mid(c), mid(T), N, mid(d)) 
U0_big = Sequence(fourier_big,interval.(big.(coefficients(U0))))
U0 = interval.(U0)

#################################################   Projection on X⁴₀   ##################################################################################
# Projection of U0 in X⁴₀ : U0 needs to represent a function in H⁴₀(Ω)
# We define S as the trace operator (SU = 0 means that U ∈ X⁴₀) and Sᵀ as its adjoint
S = trace(N,d,0); Sᵀ = S' 
S_big = trace(N,d_big,1); Sᵀ_big = S_big' 

# We build the operator L and its inverse Lᵢ. The use of  matrix Lᵢ is described in Theorem 5.1 and Remark 5.2  
D² = project(Derivative(2), fourier, fourier,Interval{Float64})
D⁴ = project(Derivative(4), fourier, fourier,Interval{Float64})
L = I + λ1*D² + λ2*D⁴
Lᵢ = diagm(ones(N+1)./diag(coefficients(L)))
Lᵢ_big = interval.(big.(diagm(ones(N+1)./diag(coefficients(mid.(L))))))

# Finally we can build the projection of U0 on X⁴₀ that we denote U0 again : U0 = U0 - DᵢSᵀinv(SDᵢSᵀ)SU0

U0_big = U0_big - Sequence(fourier_big,vec(Lᵢ_big*Sᵀ_big*inv(S_big*Lᵢ_big*Sᵀ_big)*S_big*coefficients(U0_big)))
Lᵢ = LinearOperator(fourier,fourier,Lᵢ)

U0 = Interval.(Float64.(inf.(U0_big),RoundDown),Float64.(sup.(U0_big),RoundUp) )  # go back to simple precision
U0 = Sequence(fourier,coefficients(U0))

##################################################   Computation of C  ##################################################################################
# Computation of the constant C(d). To gain precision on the computation of C, one could use multiprecision in the computations

# We first compute the roots of l(ξ) given in the beginnig of Section 6.2
w₁ = -λ1/(2*λ2)
w₂ = sqrt(abs(λ1^2-4*λ2))/(4*λ2)
a = sqrt( (sqrt(w₁^2+4*w₂^2) - w₁)/2 )
b = w₂/a


#Computation of max 1/l(ξ)
max1 = maximum([1, 1/(1 + λ1*(-λ1/(2λ2)) + λ2*(-λ1/(2λ2))^2)])



# # We define an operator D that help us to switch between the cosine and exponential series
# # (as the theoretical analysis is done in exponential series)
# # For a linear operator B between cosine fourier series, D*B*inv(D) gives the equivalent operator
# # on exponential series for the 0:N modes (the negative modes can be found using the even symmetry)
# # In particular, if B is diagonal, then D*B*inv(D) = B
D = interval(sqrt(2))*interval.(Diagonal(vec(ones(N+1,1))))
D[1,1] = interval(1)
Di = Diagonal(ones(N+1)./diag(D))
D = LinearOperator(fourier,fourier,D)
Di = LinearOperator(fourier,fourier,Di)

#######################################  Computation of the C bound ######################################################

# We build the operator Q(U0) = DF(U0) - L and the operator DF(U0)
DG = project(Multiplication(2*λ3*U0),fourier, fourier,Interval{Float64})      # multiplication operator by 2λ3U0
DF = L + DG       # Jacobien of F at U0

# Then we construct the operator A and compute its norm
A = interval.(inv(mid.(DF))) ;     # A is an approximated inverse for DF(U0)
norm_A = maximum([1 opnorm(LinearOperator(coefficients(D*L*A*Di)),2)])       # Computation for an upper bound of the norm of A



#Computation of the bound Z11

fourierE = CosFourier(2*N, π/d)

E = Sequence(fourierE ,interval.(zeros(2*N+1)))
for n = 1:2N  
    E[n] = real((-1)^n*( (1-exp(-4*a*d))/(2*a-im*n*π/d) ))*(1-exp(-4*a*d))
end
E[0] = 1/(2*a)*(1-exp(-4*a*d)); E = E/(2*d) ;

V0 = 2λ3*U0
C0 = (a+b)/(4*abs(λ2)*a*b*(a^2+b^2))
Cd = 4*d + 2*exp(-a*d)/(a*(1-exp(-3*a*d/2))) + 2/(a*(1-exp(-2*a*d)))

norm_dot = maximum([0 sum(coefficients((D*V0).*(D*project(V0*E,fourier))))])
Z11 = sqrt(norm_A^2*C0^2*2*d*norm_dot*( 2/a + exp(-2*a*d)*Cd))





DG_full = project(Multiplication(2*λ3*U0),fourier, CosFourier(2*N, π/d),Interval{Float64}) 
# Computation of Z defined in Lemma ...
Z1 = opnorm(D*L*(I  - A*DF)*Lᵢ*Di,2) +   abs((1+norm_A)*(norm(2*λ3*U0,1))/coefficients(L)[N+1,N+1])  + opnorm(LinearOperator(coefficients(D*(DG_full-DG)*Lᵢ*Di)),2)

Z = sup(Z1 + Z11)


# ################ Y BOUND ######################################################
# Computation of the Y bound, defined in Theorem 5.6
U2 = U0*U0
Y0 = sqrt(2*d)*sqrt( norm( L*A*project(L*U0+λ3*U2,fourier) ,2)^2 + norm(U2-project(U2,fourier),2)^2 )
# display("value of Y0")
# display(Y0)
#
################################ Z2 BOUND ######################################################

# Computation of the constants κ and C(r) 
κ = maximum([1,1/(1+λ1*(sqrt(abs(λ1)/λ2))^2 + λ2*(sqrt(abs(λ1)/λ2))^4)])*maximum([1,sqrt(2/(2+λ1/λ2))])*sqrt(3/(8*sqrt(interval(2))*λ2^(1/4)))

# Make sure that κ is big enough for the proof of the periodic solution 
κ = maximum([κ 1/sqrt(2d)*sqrt(1+d/π)*κ])

Cr = 2*λ3*κ

# Computation of the bound Z2 
Z2 = norm_A*Cr
#
# # # ############### Verification ###################################################
# In this section we verify that the hypotheses of the radii-polynomial approach are satisfied. If it is the case 
# we display the smallest and the largest radius for which we obtain a contraction in Hˡ
if 1- Z>0
  if 1-sup(4*Y0*Z2) > 0
    rmin=sup((1-Z - sqrt((1-Z)^2-4*Y0*Z2))/(2*Z2))
    rmax=inf((1-Z + sqrt((1-Z)^2-4*Y0*Z2))/(2*Z2))
    if rmin<rmax
      display("The computer-assisted proof for the soliton and the periodic solution was successful !")
      display("Value of the minimal radius for the contraction")
     display(rmin)
      # display("Value of the maximal radius for the contraction")
      # display(rmax)
      
    else
      display("rmin>=rmax")
      return(0)
    end
  else
    display("failure: discriminant is negative")
    return(0)
  end
else
    display("failure: 1-z > 0")
    return(0)
end




####################### PROOF OF EIGENCOUPLES ######################################################################################################################################

# We build the full exponential series from the cosine series
U0_f = coefficients(U0)
U0_f = [reverse(U0_f[2:N+1]); U0_f[1] ; U0_f[2:N+1]]
U0_f_big = coefficients(U0_big)
U0_f_big = [reverse(U0_f_big[2:N+1]); U0_f_big[1] ; U0_f_big[2:N+1]]

fourier = Fourier(N,π/d)
U0_f = Sequence(fourier,U0_f)
U0_f_big = Sequence(fourier,U0_f_big)

# We first build the operators L and Q(U0)
D² = project(Derivative(2), fourier, fourier,Complex{Interval{Float64}})
D⁴ = project(Derivative(4), fourier, fourier,Complex{Interval{Float64}})
Lf = I + λ1*D² + λ2*D⁴         # operator L
Qf = project(Multiplication(2*λ3*U0_f),fourier, fourier,Complex{Interval{Float64}})      # multiplication operator by 2λ3U0

# We make sure that Df is Hermitian (true theoretically)
Df = Hermitian(coefficients(mid.(Lf + Qf)))
E = eigen(Df,1:3)   # computation of the first 3 eigencouples 

pspace = ParameterSpace()×fourier

# Construction of the eigencouples
ν1 = E.values[1] ; ν2 = E.values[2] ; ν3 = E.values[3]
ν2 = 0 # Theoretically we expect the second eigenvalue to be exactly zero with eigenvector (u^*)' 
V1 = E.vectors[:,1] ; V2 = E.vectors[:,2] ; V3 = E.vectors[:,3] 
U1 = Sequence(pspace,vec([ν1;V1])) ; U2 = Sequence(pspace,vec([ν2;V2])) ; U3 = Sequence(pspace,vec([ν3;V3]))
U1 = interval.(real(U1)) + im*interval.(imag(U1))  
U2 = interval.(real(U2)) + im*interval.(imag(U2))
U3 = interval.(real(U3)) + im*interval.(imag(U3))


# Proof of the eigencouples !
P1 = proof_eigen_kawahara(U1,U0_f,λ1,λ2,λ3,N,d,rmin)
P2 = proof_eigen_kawahara(U2,U0_f,λ1,λ2,λ3,N,d,rmin)
P3 = proof_eigen_kawahara(U3,U0_f,λ1,λ2,λ3,N,d,rmin)


########################### PROOF THAT THE PREVIOUS EIGENCOUPLES ARE THE FIRST 3 ONES #####################################################################################################################################

# Computation of a and b in the case ν = ν3
λ1k = λ1/(1-ν3) ; λ2k = λ2/(1-ν3) ; λ3k = λ3/(1-ν3)
w₁ = -λ1k/(2*λ2k)
w₂ = sqrt(abs(λ1k^2-4*λ2k))/(4*λ2k)
a = sqrt( (sqrt(w₁^2+4*w₂^2) - w₁)/2 )
b = w₂/a

# Computation of Z11 in the case in the case ν = λ3. This quantity will be usefull in the next proof (no eigenvalue in (ν_i, ν_{i+1})) and no eigen smaller that ν1)
fourierE = Fourier(2*N, π/d)

E = Sequence(fourierE ,interval.(zeros(4*N+1)))
    for n = [-2N:-1 1:2N]  
        E[n] = real((-1)^n*( (1-exp(-4*a*d))/(2*a-im*n*π/d) ))*(1-exp(-4*a*d))
    end
    E[0] = 1/(2*a)*(1-exp(-4*a*d)); E = E/(2*d) ;
    
    V0 = 2*λ3k*U0_f
    C0 = (a+b)/(4*abs(λ2k)*a*b*(a^2+b^2))
    Cd = 4*d + 2*exp(-a*d)/(a*(1-exp(-3*a*d/2))) + 2/(a*(1-exp(-2*a*d)))
    
    norm_dot = maximum([0 sum(coefficients((V0).*(project(V0*E,fourier))))])
    Z11 = sqrt(C0^2*2*d*norm_dot*( 2/a + exp(-2*a*d)*Cd))# return to simple precision

# Proof of Theorem 6.20
P = no_eigen_outside(P1[2],P2[2],P3[2],U1,U2,U3,U0_f,Z11,rmin,λ1,λ2,λ3,N)

if P==1
  display("There is no eigenvalue in between ν1 and ν2 and between ν2 and ν3")
  display("Moreover, there is no eigenvalue smaller than ν1")
else
  display("The inverse becomes singular for some value of ν")
end

############################### PROOF OF ORBITAL STABILITY ####################################################################################################

# Proof of stability of the soliton (works only if 1/3 <= T <= 0.397...) using Proposition 6.23.
if (P1[1]==P2[1]==P3[1]==P)&&(ν1<0)&&(ν3>0)
    V0 = A*U0
    C1 = norm_A/(1-Z)*1/(1-2*abs(λ3)*norm_A/(1-Z)*rmin)
    
    V =  V0 -  Sequence(CosFourier(N, π/d),coefficients(Lᵢ)^4*Sᵀ*inv(S*coefficients(Lᵢ)^4*Sᵀ)*S*coefficients(V0))

    ϵ = C1*sqrt(interval(2*d))*norm(D*(U0-(L+DG)*V),2) + 2*λ3*κ*sqrt(interval(2*d))*C1*rmin*norm(D*L*V,2)

    τ = 2*d*dot(coefficients(D*U0),coefficients(D*V)) + ϵ*sqrt(interval(2*d))*norm(U0,2) + 2*C1*(sqrt(interval(2*d))*norm(U0,2)+rmin)*rmin 


    if τ<0
      #display("P1, P2 and P3 are satisfied !")
      display("The soliton is stable and a value for τ is")
      display(sup(τ))
    else
      display("τ is positive, the test is inconclusive")
    end

  else 
    display("Assumption P1, P2 or P3 is not sastisfied")
end


display("Values in the computer-assisted proof of the soliton")

header1 = (["||DF(u0)^{-1}||_{2,l}", "Z1", "Zu", "Y0", "Minimal radius of contraction", "Maximal value of contraction"])
data = [sup(C1) sup(Z1) sup(Z11) sup(Y0) sup(rmin) inf(rmax) ]

pretty_table(data;
             header = header1
             )


display("Values in the computer-assisted proof of the eigencouples")

header = (["Eigenvalue", "Minimal radius of contraction", "Maximal value of contraction"])
data = [mid(ν1)  sup(P1[3]) inf(P1[2]);
        mid(ν2)  sup(P2[3]) inf(P2[2]);
        mid(ν3)  sup(P3[3]) inf(P3[2])]
             
pretty_table(data;
            header = header )
