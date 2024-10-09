function mdotmax(P0,T0,A,gamma,R)
    return ( P0*A*sqrt(gamma)*((gamma+1)/2)^(-1*(gamma+1)/(2*(gamma-1))) )/sqrt(R*T0)
end

function mdotfrommach(M,P0,T0,A,gamma,R)
    return ( P0*A*sqrt(gamma)*M*(1+((gamma-1)/2)*M^2)^(-1*(gamma+1)/(2*(gamma-1))) )/sqrt(R*T0)
end

function pstar(P0,gamma)
    return P0 * ((gamma+1)/2)^(-1*gamma/(gamma -  1))
end

function machfrompressure(P0,P,gamma)
    return sqrt( (2/(gamma-1)) * ((P0/P)^((gamma-1)/gamma) - 1) )
end

function pressurefrommach(M,P0,gamma)
    return P0*(1+((gamma-1)/2)*M^2)^(-1*(gamma)/(gamma-1))
end

function astarfrommach(M,A,gamma)
    return A*M* ((2/(gamma+1))*(1+ ((gamma-1)/2)*M^2))^(-1*(gamma+1)/(2*(gamma-1)))
end

function machfromarea(arearatio,gamma,supersonic::Bool)
    f(M) = ( (1/M)* ((2/(gamma+1))*(1+ ((gamma-1)/2)*M^2))^((gamma+1)/(2*(gamma-1))) ) - arearatio
    #=fx = ZeroProblem(f, (0,1))
    if supersonic
        fx = ZeroProblem(f, (1,1e3))
    end
    return solve(fx, Bisection(); xatol=1/16)=#
    if supersonic
        return fzero(f,1,1e3)
    else
        #println(f(0))
        #println(f(1))
        return fzero(f,0,1)
    end
end
