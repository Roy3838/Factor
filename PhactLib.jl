"""
    customGCD(a::T, r::T, N::T, k::T=1) where T <: Union{Int32, Int64, BigInt}

Returns the point (a,r) when it is not a trivial case
"""
function customGCD(a::Union{Int64, Int128,  BigInt}, r::Union{Int64, Int128,  BigInt}, N::Union{Int64, Int128,  BigInt}, k::Union{Int64, Int128,  BigInt}=1)
    # Define/Preallocate
    #f = zero(T)
    #g1 = zero(T)
    #g2 = zero(T)
    #gf = zero(T)
    #g = zero(T)
    #gi = zero(T)

    f = powermod(a, k*r, N)
    g1 = gcd(f+1,N)
    g2 = gcd(f-1,N)
    gf = max(g1,g2)
    g = gf*(1<gf<N) # this is the returned value!
    #g = 1*(1<gf<N) # this is the returned value!
    return g
end

"""
    FNFfact(N::T, criticN::T) where T <: Union{Int32, Int64, BigInt}

Return a non-trivial factor if succeed
"""
function FNFfact(N::Union{Int64, Int128,  BigInt}, b::Union{Int64, Int128,  BigInt})

#    a0 = div(N,2)
    a0 = one(typeof(N)) + one(typeof(N))
    #r0 = sN
    a2 = N
    #r2 = div(r0,2)
    #k = 1
    factorP = zero(typeof(N))
    counter = zero(typeof(N))

    # ie the while loop more efficient compared to the for-loop?
    while factorP==0 && counter < N
        counter = counter+1
        a = a0 + counter*BigInt(exp10(b))
        factorP = customGCD(a, a, N)
    end
    return factorP, counter
end

"""
    FactAnalysis(N::T) where T <: Union{Int32, Int64, BigInt}   **** FALTA ARREGLAR!

Finds the first argument of the function f(x) that throws a non-trivial factor.
"""
function FactAnalysis(N::T) where T <: Union{Int32, Int64, BigInt}

maxN=size(N,1);
y=zeros(T,maxN);

    for x=1:maxN,
        limN=N[x];
        for i=1:limN;
            G=customGCD((i+isqrt(limN)),(i+isqrt(limN)),limN)
    #         G=customGCD(i,i,limN)
    #         println(G);
            if G!=0
                y[x]=i;
    #             println(y[x])
                break;
            end
        end
    end
return y;
end
