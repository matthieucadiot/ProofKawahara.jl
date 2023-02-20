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


#################################### List of the needed functions : go directly to line 504 for the main code ################################################# 

function newton_kawahara(c, T, N, d)

    fourier = CosFourier(N,π/d)
    
    X= 2*d/(2*N+1).*(-N:N)
    u = 2*(c-1)./((cosh.(sqrt(abs(3*(c-1)/(2*(1-3*T))))*X)).^2)
    V = real(FFTW.fftshift(FFTW.fft(FFTW.ifftshift(u))/(2*N+1)))
    V = vec(V[N+1:2*N+1])
    U = Sequence(fourier, vec(real(V)))

    D² = project(Derivative(2), fourier, fourier,Float64)
    D⁴ = project(Derivative(4), fourier, fourier,Float64)
    λ2 = (19 - 30T - 45T^2)/(360*(1-c)) ;  λ1 = (1-3T)/(6*(1-c)) ; λ3 = 3/(4*(1-c))

    L = I +λ1*D² + λ2*D⁴ # petite modification pour minimiser le nombre d'opération
    G = project(Multiplication(λ3*U),fourier, fourier,Float64)
    F = L*U +G*U

    n_F = norm(F,2)
    n_max = 20 ; p =0;
    err = 1.0e-16
    while (n_F > err)&&(p <n_max)
        DF = L + 2*G
        U = U - DF\F
        G = project(Multiplication(λ3*U),fourier, fourier,Float64)

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
        w = interval.((0:N))
        Y = (interval((pi))^2/interval((d))^2)*w.^2
        X = interval.(((2*(-1).^(0:N)))); X[1] = interval((1))

        f = [X';
            (X.*Y)']
        return f
    else
        w = interval.(big.(0:N))
        Y = (interval(big(pi))^2/interval(big(d))^2)*w.^2
        X = interval.(big.((2*(-1).^(0:N)))); X[1] = interval(big(1))

        f = [X';
            (X.*Y)']
        return f 
    end
end


function trace_full(N,d)
    w = -N:N
    p = interval(pi)
    d = interval(d)
    Y2 = (p/d)*interval.(w);
    Y1 = (p^2/d^2)*interval.(w).^2
    Y3 = (p^3/d^3)*interval.(w).^3
    X = interval.((-1).^w)

    f = [X';
        (X.*Y1)';
         (X.*Y2)';
         (X.*Y3)']
    return f
end


function compute_boundary(V,N,d,a)
    S = 0.

    for n = 1:N
    
        for m = 1:N
               b = 4*a*(1/(((π)/d*(n-m))^2+4*a^2) + 1/((π/d*(n+m))^2+4*a^2))
               u= V[n]
               v= V[m]
               S = S+ ((-1)^(n-m))*u*v*b
         end
               b = 4*a/((π/d*(n))^2+4*a^2) 
               u= V[n]
               v= V[0]
               S = S+(((-1)^(n)))*u*v*b
    end

    for m = 1:N
           b = 4*a/((π/d*(m))^2+4*a^2)  
           u= V[0]
           v= V[m]
           S = S+ (((-1)^(m)))*u*v*b
    end
    b =  1/(2*a) 
    u= V[0]
    v= V[0]
    S = S+ u*v*b

    f = sqrt(abs(S))
end




function compute_boundary_1(V,N,d,a)
    S = 0.


    for n = 1:N
        for m = 1:N
               b = 4*a*(1/((π/d*(n-m))^2+4*a^2) - 1/((π/d*(n+m))^2+4*a^2)) 
               u= (n)*π/d*V[n]
               v= (m)*π/d*V[m]
               S = S+ (((-1)^(n-m)))*u*v*b
        end

    end

    f = sqrt(abs(S))
end


function compute_boundary_total(V,N,d,a,b,λ2,λ1,max1,max2)

# Computation of the constants of Lemma 6.6
        Cv0 = 2*compute_boundary(V,N,d,a)             #computation of C(v₀)e^{-ad} 
        Cv02 = 2*compute_boundary(V*V,2*N,d,a)        #computation of C(v₀²)e^{-ad}  
        Cv0d =  2*compute_boundary_1(V,N,d,a)         #computation of C(v₀')e^{-ad}  

        δ₀  = (a+b)/(4*λ2*a*b*(a^2+b^2))
        δ₁ = (a^3 + a^2*b + b^2*a + b^3)/(4*λ2*a*b*(a^2+b^2))

        
        C1 = 1/(4*λ2*a*b)*(max1*Cv0 + (2*d+1/a)*δ₀*Cv0 + (2*d+1/a)*δ₀*Cv02)
        C3 = δ₁*(max2*Cv0 + Cv0d + (2*d+1/a)/(4*λ2*a*b)*Cv0 + (2*d+1/a)*Cv02/(4*λ2*a*b))

# Computation of b1ᴺ
        N = 100;
        b1ᴺ = λ1^2
        for n = 1:100
            tn = n*π/d
            b1ᴺ = b1ᴺ + 2*(λ1 + λ2*tn^2)^2/(1+(1+λ1*tn^2 + λ2*tn^4)^2)
        end


# Computation of the constants of lemma 6.3 
        
        b1 = (1+max1)/d*( b1ᴺ + 2d^2/(3*π^4*N^3) )^(1/2)
        b3 = (1+max1)*λ2/d*( 1/2 + 3*π*d/(4*sqrt(interval(2))*λ2^(1/4))*maximum([1,2/(2+λ1/sqrt(λ2))]) )^(1/2)

        Z1 = max(1,norm(V,1))*( b1*C1 + b3*C3)

    return abs(Z1)
end





function compute_boundary_eig(V,N,d,a,b,λ2,λ1,max1)


    d0 = abs((a+b)/(4*λ2*a*b*(a^2+b^2)))
    d1 = abs((a^3 + a^2*b + a*b^2 + b^3)/(4*λ2*a*b*(a^2+b^2)))
        
    W = Sequence(CosFourier(N, π/d),V[0:N])
    W2 = Sequence(CosFourier(2*N, π/d),(V*V)[0:2*N])

    Cv0 = 2*compute_boundary(W,N,d,a)
    Cv02 =  2*compute_boundary(W2,2*N,d,a)
    Cv0d =  2*compute_boundary_1(W,N,d,a)


    C0 = d0*( Cv0 + (2*d+1/a+1)*d0*Cv0 + (2*d+1/a+1)*d0*Cv02 )
    C1 = 1/(4*λ2^2*a*b)*( Cv0 + (2*d+1/a+1)*d0*Cv0 + (2*d+1/a+1)*d0*Cv02 )
    C2 = d1*( Cv0 + (2*d+1/a+1)*d0*Cv0 + (2*d+1/a+1)*d0*Cv02 )
    C3 = d1*( Cv0*max1 + Cv0d + (2*d+1/a+1)/(4*λ2*a*b)*Cv0 + (2*d+1/a+1)*Cv02/(4*λ2*a*b) )

    # Computation of bᴺ0, bᴺ1 and bᴺ2

    N = 100;
    b0ᴺ = interval(0.)
    b1ᴺ = λ1^2
    b2ᴺ = interval(0.)
    for n = 1:100
        tn = n*π/d
        b0ᴺ = b0ᴺ + 2*tn^2*(λ1 + λ2*tn^2)^2/(1+(1+λ1*tn^2 + λ2*tn^4)^2)
        b1ᴺ = b1ᴺ + 2*(λ1 + λ2*tn^2)^2/(1+(1+λ1*tn^2 + λ2*tn^4)^2)
        b2ᴺ = b2ᴺ + 2*(λ2*tn)^2/(1+(1+λ1*tn^2 + λ2*tn^4)^2)
    end


    b0 = 1/d*( b0ᴺ + d^2/(π^2*N) )^(1/2)
    b1 = 1/d*( b1ᴺ + 2*d^2/(3*π^4*N^3) )^(1/2)
    b2 = 1/d*( b2ᴺ +  d^6/(5*π^6*N^5))^(1/2)
    b3 = λ2/d*( 1/2 + 3*π*d/(4*sqrt(interval(2))*λ2^(1/4)) )^(1/2)

    Z1 = 10/3*max(1,norm(V,1))*(b0*C0 + b1*C1 + b2*C2 + b3*C3)

    return abs(Z1)
end



function proof_eigen_kawahara(U,U0,λ1,λ2,λ3,N,d,rmin)

    pspace = space(U)
    fourier = space(U0)

    V = component(U,2) 
    ν = real(component(U,1)[1])

    D² = project(Derivative(2), fourier, fourier,Complex{Interval{Float64}})
    D⁴ = project(Derivative(4), fourier, fourier,Complex{Interval{Float64}})
    L = 1/(1-ν)*interval.(real.((I + λ1*D² +λ2*D⁴)))
    Lᵢ = diagm(ones(2*N+1)./diag(coefficients(L)))
    S = trace_full(N,d); Sᵀ = S' 

    V = V - Sequence(fourier,Lᵢ*Sᵀ*inv(S*Lᵢ*Sᵀ)*S*coefficients(V))

    λ1k = λ1/(1-ν) ; λ2k = λ2/(1-ν) ; λ3k = λ3/(1-ν)

    U = Sequence(pspace,vec([ν; coefficients(V)]))
    # display("projection on X⁴₀ created")

    # Computation of a and b
    x = -λ1k/(2*λ2k)
    y = sqrt(abs(λ1k^2-4*λ2k))/(4*λ2k)

    a = sqrt( (sqrt(x^2+4*y^2) - x)/2 )
    b = y/a

    

    # Computation of maximum of |ξ|/l(ξ) 
    x = sqrt( abs((sqrt(λ1k^2+12λ2k) - λ1k)/(6λ2k)) )
    max2 = 1/(2λ1k*x + 4λ2k*x^3)

    Z1 = compute_boundary_eig(2*λ3k*U0,N,interval(big(d)),a,b,λ2k,λ1k,max2)

    λ3k = Interval.(Float64.(inf.(λ3k),RoundDown),Float64.(sup.(λ3k),RoundUp) )
    ν = Interval.(Float64.(inf.(ν),RoundDown),Float64.(sup.(ν),RoundUp) )
    U0 = Interval.(Float64.(inf.(U0),RoundDown),Float64.(sup.(U0),RoundUp) )
    L = Interval.(Float64.(inf.(real(L)),RoundDown),Float64.(sup.(real(L)),RoundUp) )
    Lᵢ = Interval.(Float64.(inf.(real(Lᵢ)),RoundDown),Float64.(sup.(real(Lᵢ)),RoundUp) )
    V = Interval.(Float64.(inf.(real(V)),RoundDown),Float64.(sup.(real(V)),RoundUp) ) + im*Interval.(Float64.(inf.(imag(V)),RoundDown),Float64.(sup.(imag(V)),RoundUp) )
    Z1 = Interval.(Float64.(inf.(Z1),RoundDown),Float64.(sup.(Z1),RoundUp) )
    

    Lν = L - ν/(1-ν)*I

     # We make sure that ||V||ₗ ≤ 1 in order to have ||v||ₗ ≤ 1
    V = inf(1/norm(Lν*V,2))*V

    Q = project(Multiplication(2*λ3k*U0),fourier, fourier,Complex{Interval{Float64}})

    T = LinearOperator(pspace,pspace, [[0   -(coefficients(Lν*V))'] ; [-coefficients(V) coefficients(Lν+Q)]] )
    D = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) coefficients(Lν)]]
    Di = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) Diagonal(ones(2*N+1)./diag(coefficients(Lν)))]]
    D = LinearOperator(pspace,pspace,D)
    Di = LinearOperator(pspace,pspace,Di)
   

# # # # ################ Z BOUND ######################################################

    A = inv(mid.(T)) ; A = interval.(real.(A)) + im*interval.(imag.(A))
    norm_A = maximum([1 sup(opnorm(LinearOperator(coefficients(D*A)),2))]);
    
    Q_full = project(Multiplication(2*λ3k*U0),fourier, Fourier(2*N, π/d),Complex{Interval{Float64}})
    Lᵢ = LinearOperator(fourier,fourier,Lᵢ)
    
    Z = opnorm(LinearOperator(coefficients(D*(I  - A*T)*Di)) ,2) + abs((1+norm_A)*(norm(2*λ3k*U0,1))/coefficients(Lν)[2*N+1,2*N+1]) +  opnorm(LinearOperator(coefficients((Q_full-Q)*Lᵢ)),2)
    

    C1 = maximum( [2 ; norm_A/(1-Z)] )
    

    κ = sqrt(3*π/(4*sqrt(interval(2))*λ2k^(1/4)))
    r0 = rmin
    r1 = 2*λ3k*κ*norm(V,1)*r0
    r2 = 2*λ3k*κ*r0

    C = interval(sqrt(3))*C1/(1-Z1*C1)
    C = C/(1-r2*C)
     

################# Y BOUND ######################################################

    Y = sqrt(2*interval(d))*C*(r1 + norm((Lν*V+2*λ3k*U0*V),2))
    
################## Z2 BOUND ######################################################

    Cr = 1/(1-ν)
    Z2 = C*Cr
    
################## Verification ###################################################
# # # r_star = 1.0;
# # # interval_proof = interval_of_existence(Y,Interval(0),Z2,r_star,C²Condition())
# # # display(interval_proof)
# # #
x = 1;
if 1- Z>0
    if 1-sup(4*Y*Z2) > 0
      rmin = Float64.(sup.((1 - sqrt(1-4*Y*Z2))/(2*Z2)),RoundUp)
      rmax = Float64.(inf.((1 + sqrt(1-4*Y*Z2))/(2*Z2)),RoundDown)
      if rmin<rmax
        display("The computer-assisted proof was successful for the eigenvalue ν with value")
        display(Float64.(mid(ν)))
        #display("Value of the minimal radius for the contraction")
        #display(Float64(rmin))
        #display("Value of the maximal radius for the contraction")
        #display(Float64(rmax))
        
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




function no_eigen_outside(r1,r2,r3,U1,U2,U3,U0,λ1,λ2,λ3,N,Z1)

    fourier = space(U0)

    ν1 = component(real(U1),1)[1];     ν2 = component(real(U2),1)[1];    ν3 = component(real(U3),1)[1]
    ν1 = Interval.(Float64.(inf.(ν1),RoundDown),Float64.(sup.(ν1),RoundUp) ) 
    ν2 = Interval.(Float64.(inf.(ν2),RoundDown),Float64.(sup.(ν2),RoundUp) ) 
    ν3 = Interval.(Float64.(inf.(ν3),RoundDown),Float64.(sup.(ν3),RoundUp) ) 


    D² = project(Derivative(2), fourier, fourier,Complex{Interval{Float64}})
    D⁴ = project(Derivative(4), fourier, fourier,Complex{Interval{Float64}})

    L = interval.(real.((I + λ1*D² +λ2*D⁴)))
    Q = real.(project(Multiplication(2*λ3*U0),fourier, fourier,Complex{Interval{Float64}}))
    Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(L))))

    ν = interval(inf(ν1+r1))
    P = 1
    k=1

    # Computation of the non-changing upper bounds
    # Computation of \|(π^{2N} - π^N)DG(U₀)L^{-1}\|
    Q_full = project(Multiplication(2*λ3*U0),fourier, Fourier(2*N, π/d),Complex{Interval{Float64}})
    nQ = opnorm(LinearOperator(coefficients((Q_full-Q)*Li)),2) 

    # Computation of \|2λ3U0\|_1
    n0 = norm(2*λ3*U0,1)
    
    while  (P==1)&&(inf(ν)<sup(λ2-r2))
        Lν = 1/(1-ν)*(L - interval(ν)*I)
        Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(Lν))))

        Df = Lν + 1/(1-ν)*Q
        A = inv(mid.(Df))
        norm_A = maximum([1 opnorm(LinearOperator(coefficients(Lν*A)),2)])       # Computation for an upper bound of the norm of A

        Z = opnorm(Lν*(I  - A*Df)*Li,2)+   1/(1-ν)*abs((1+norm_A)*(n0)/coefficients(Lν)[2*N+1,2*N+1]) + nQ
        if sup(Z) >= 1
            P=0
            return P==1
        end
        
        C1 = maximum([2,norm_A/(1-Z)])

        if sup(Z1*C1) >= 1
            P=0
            return P==1
        end

        C = C1/(1 - C1*Z1)
        ν = ν + inf(1/C)
    end

    ν = interval(inf(ν2+r2))
    k=2
    
    while  (P==1)&&(inf(ν)<sup(ν3-r3))
        
        Lν = 1/(1-ν)*(L - interval(ν)*I)
        Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(Lν))))

        Df = Lν + 1/(1-ν)*Q
        A = inv(mid.(Df))
        norm_A = maximum([1 opnorm(LinearOperator(coefficients(Lν*A)),2)])       # Computation for an upper bound of the norm of A

        Z = opnorm(Lν*(I  - A*Df)*Li,2)+   1/(1-ν)*abs((1+norm_A)*(n0)/coefficients(Lν)[2*N+1,2*N+1]) + nQ
        if sup(Z) >= 1
            P=0
            return P==1
        end

        C1 = maximum([2,norm_A/(1-Z)])
        
        if sup(C1*Z1) >= 1
            P=0
            return P==1
        end

        C = C1/(1 - C1*Z1)
        ν = ν + inf(1/C)
    end

    ν = interval(inf(ν1-r1))
    nw = norm(2*λ3*L*U0,1)
    k=3
    while  (P==1)&&(1-sup(ν)<sup(nw+1/2))
        
        Lν = 1/(1-ν)*(L - interval(ν)*I)
        Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(Lν))))

        Df = Lν + 1/(1-ν)*Q
        A = inv(mid.(Df))
        norm_A = maximum([1 opnorm(LinearOperator(coefficients(Lν*A)),2)])      # Computation for an upper bound of the norm of A

        Z = opnorm(Lν*(I  - A*Df)*Li,2)+   1/(1-ν)*abs((1+norm_A)*(n0)/coefficients(Lν)[2*N+1,2*N+1]) + nQ
        if sup(Z) >= 1
            P=0
            return P==1
        end
        
        C1 = maximum([2,norm_A/(1-Z)])

        if sup(Z1*C1) >= 1
            P=0
            return P==1
        end

        C = C1/(1 - C1*Z1)
        ν = ν - inf(1/C)
    end

    return P==1


end



################### PROOF OF SOLITON AND PERIODIC SOLUTIONS : MAIN CODE #################################################################################################################################################

N = 250              # number of Fourier modes : from 0 to N for a cosine series 
d= 50          # size of the domain = half period of the functions
T = 0.35 ; T_big = interval(big(T)) ; T = interval(T)  # Bond number
c = 0.9   ; c_big = interval(big(c)) ; c = interval(c) # velocity of the soliton
fourier = CosFourier(N, π/d)   # definition of the sequence space : cosine series of frequency π/d

λ2 = (19 - 30T - 45T^2)/(360*(1-c)) ;  λ1 = (1-3T)/(6*(1-c)) ; λ3 = 3/(4*(1-c))  # definition of the parameters of the equation
λ2_big = (19 - 30T_big - 45T_big^2)/(360*(1-c_big)) ;  λ1_big = (1-3T_big)/(6*(1-c_big)) ; λ3_big = 3/(4*(1-c_big))

# Computation of the approximated solution via Newton's method

U0 = newton_kawahara(mid(c), mid(T), N, d) 
U0_big = interval.(big.(U0))
U0 = interval.(U0)

# fig = lines(LinRange(0, 2pi, 100), t -> real(solution(t)[0]))
# save(fig, “my_plot.png”)


#################################################   Projection on X⁴₀   ##################################################################################
# Projection of U0 in X⁴₀ : U0 needs to represent a function in H⁴₀(Ω)
# We define S as the trace operator (SU = 0 means that U ∈ X⁴₀) and Sᵀ as its adjoint
S = trace(N,d,0); Sᵀ = S' 
S_big = trace(N,d,1); Sᵀ_big = S_big' 

# We build the operator L and its inverse Lᵢ. The use of  matrix Lᵢ is described in Theorem 5.1 and Remark 5.2  
D² = project(Derivative(2), fourier, fourier,Interval{Float64})
D⁴ = project(Derivative(4), fourier, fourier,Interval{Float64})
L = I + λ1*D² + λ2*D⁴
Lᵢ = diagm(ones(N+1)./diag(coefficients(L)))
Lᵢ_big = interval.(big.(diagm(ones(N+1)./diag(coefficients(mid.(L))))))

# Finally we can build the projection of U0 on X⁴₀ that we denote U0 again : U0 = U0 - DᵢSᵀinv(SDᵢSᵀ)SU0

U0_big = U0_big -  Sequence(fourier,vec(Lᵢ_big*Sᵀ_big*inv(S_big*Lᵢ_big*Sᵀ_big)*S_big*coefficients(U0_big)))
Lᵢ = LinearOperator(fourier,fourier,Lᵢ)

U0 = Interval.(Float64.(inf.(U0_big),RoundDown),Float64.(sup.(U0_big),RoundUp) )  # go back to simple precision

##################################################   Computation of C  ##################################################################################
# Computation of the constant C(d). To gain precision on the computation of C, one could use multiprecision in the computations

# We first compute the roots of l(ξ) given in the beginnig of Section 6.2
w₁ = -λ1/(2*λ2)
w₂ = sqrt(abs(λ1^2-4*λ2))/(4*λ2)
a = sqrt( (sqrt(w₁^2+4*w₂^2) - w₁)/2 )
b = w₂/a

# Computation of max ξ/(l(ξ)) 
x = sqrt( abs((sqrt(λ1^2+12λ2) - λ1)/(6λ2)) )
max2 = 1/(2λ1*x + 4λ2*x^3)

#Computation of max 1/l(ξ)
max1 = maximum([1, 1/(1 + λ1*(-λ1/(2λ2)) + λ2*(-λ1/(2λ2))^2)])

# We compute the value of Z₁
Z1 = compute_boundary_total(2*λ3*U0,N,interval(d),a,b,λ2,λ1,max1,max2)


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
Q = project(Multiplication(2*λ3*U0),fourier, fourier,Interval{Float64})      # multiplication operator by 2λ3U0
DF = L + Q        # Jacobien of F at U0

# Then we construct the operator A and compute its norm
A = interval.(inv(mid.(DF))) ;     # A is an approximated inverse for DF(U0)
norm_A = maximum([1 opnorm(LinearOperator(coefficients(D*L*A*Di)),2)])       # Computation for an upper bound of the norm of A

Q_full = project(Multiplication(2*λ3*U0),fourier, CosFourier(2*N, π/d),Interval{Float64}) 
# Computation of Z defined in Lemma ...
Z = opnorm(D*L*(I  - A*DF)*Lᵢ*Di,2) +   abs((1+norm_A)*(norm(2*λ3*U0,1))/coefficients(L)[N+1,N+1])  + opnorm(LinearOperator(coefficients(D*(Q_full-Q)*Lᵢ*Di)),2)
#display(Z)

# Computation of C1 
C1 = maximum( [2 ; norm_A/(1-Z)] )

# Finally we can compute an upper bound for the norm of the inverse of L + DG(u0) that we denote by C
C = C1/(1-Z1*C1)
#display("value of C")
#display(sup(C))
#
# ################ Y BOUND ######################################################
# Computation of the Y bound, defined in Theorem 5.6
Y = sqrt(2*interval(d))*C*norm( L*U0+λ3*U0*U0 ,2)
#display("value of Y")
#display(sup(Y))
#
################################ Z2 BOUND ######################################################

# Computation of the constants κ and C(r) 
κ = maximum([1,1/(1+λ1*(sqrt(abs(λ1)/λ2))^2 + λ2*(sqrt(abs(λ1)/λ2))^4)])*maximum([1,sqrt(2/(2+λ1/λ2))])*sqrt(3*π/(4*sqrt(interval(2))*λ2^(1/4)))
Cr = 2*λ3*κ

# Computation of the bound Z2 
Z2 = C*Cr
#
# # # ############### Verification ###################################################
# In this section we verify that the hypotheses of the radii-polynomial approach are satisfied. If it is the case 
# we display the smallest and the largest radius for which we obtain a contraction in Hˡ
if 1- Z>0
  if 1-sup(4*Y*Z2) > 0
    rmin=sup((1 - sqrt(1-4*Y*Z2))/(2*Z2))
    rmax=inf((1 + sqrt(1-4*Y*Z2))/(2*Z2))
    if rmin<rmax
      display("The computer-assisted proof for the soliton and the periodic solution was successful !")
      #display("Value of the minimal radius for the contraction")
      #display(rmin)
      #display("Value of the maximal radius for the contraction")
      #display(rmax)
      
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
V1 = E.vectors[:,1] ; V2 = E.vectors[:,2] ; V3 = E.vectors[:,3] 
U1 = Sequence(pspace,vec([big(ν1);V1])) ; U2 = Sequence(pspace,vec([big(ν2);V2])) ; U3 = Sequence(pspace,vec([big(ν3);V3]))
U1 = interval.(real(U1)) + im*interval.(imag(U1))  
U2 = interval.(real(U2)) + im*interval.(imag(U2))
U3 = interval.(real(U3)) + im*interval.(imag(U3))


# Proof of the eigencouples !
P1 = proof_eigen_kawahara(U1,U0_f,λ1,λ2,λ3,N,d,rmin)
P2 = proof_eigen_kawahara(U2,U0_f,λ1,λ2,λ3,N,d,rmin)

# Sometimes, the proof might fail in standard precision. We can use multiprecision instead
P3 = proof_eigen_kawahara(big.(U3),U0_f_big,λ1_big,λ2_big,λ3_big,N,d,rmin)


########################### PROOF THAT THE PREVIOUS EIGENCOUPLES ARE THE FIRST 3 ONES #####################################################################################################################################

# Computation of maximum |ξ|/l(ξ) in the case ν = λ3
λ1k = λ1_big/(1-λ3) ; λ2k = λ2_big/(1-λ3) ; λ3k = λ3_big/(1-λ3)
x = sqrt( abs((sqrt(λ1k^2+12λ2k) - λ1k)/(6λ2k)) )
max_1 = 1/(2λ1k*x + 4λ2k*x^3)

# Computation of a and b in the case ν = λ3
w₁ = -λ1k/(2*λ2k)
w₂ = sqrt(abs(λ1k^2-4*λ2k))/(4*λ2k)
a = sqrt( (sqrt(w₁^2+4*w₂^2) - w₁)/2 )
b = w₂/a

# Computation of Z1 in the case in the case ν = λ3. This quantity will be usefull in the next proof (no eigenvalue in (ν_i, ν_{i+1})) and no eigen smaller that ν1)
Z1_λ3 = compute_boundary_eig(2*λ3k*U0_f_big,N,interval(big(d)),a,b,λ2k,λ1k,max_1)
Z1_λ3 = Interval.(Float64.(inf.(Z1_λ3),RoundDown),Float64.(sup.(Z1_λ3),RoundUp) ) # return to simple precision

# Proof of Theorem 6.20
P = no_eigen_outside(P1[2],P2[2],P3[2],U1,U2,U3,U0_f,λ1,λ2,λ3,N,Z1_λ3)

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
    
    V =  V0 -  Sequence(CosFourier(N, π/d),coefficients(Lᵢ)^4*Sᵀ*inv(S*coefficients(Lᵢ)^4*Sᵀ)*S*coefficients(V0))

    ϵ = C/(1-2λ3*κ*C*rmin)*sqrt(interval(2*d))*norm(D*(U0-(L+Q)*V),2) + 2*λ3*κ*sqrt(interval(2*d))*C/(1-2λ3*κ*C*rmin)*rmin*norm(D*L*V,2)

    τ = 2*d*dot(coefficients(D*U0),coefficients(D*V)) + ϵ*sqrt(interval(2*d))*norm(U0,2) + 2*C/(1-2λ3*κ*C*rmin)*(sqrt(interval(2*d))*norm(U0,2)+rmin)*rmin 


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

header = (["C", "Z1", "Z", "Y", "Minimal radius of contraction", "Maximal value of contraction"])
data = [sup(C) sup(Z1) sup(Z) sup(Y) sup(rmin) inf(rmax) ]

pretty_table(data;
             header = header
             )


display("Values in the computer-assisted proof of the eigencouples")

header = (["Eigenvalue", "Minimal radius of contraction", "Maximal value of contraction"])
data = [mid(ν1)  sup(P1[3]) inf(P1[2]);
        mid(ν2)  sup(P2[3]) inf(P2[2]);
        mid(ν3)  sup(P3[3]) inf(P3[2])]
             
pretty_table(data;
            header = header )
