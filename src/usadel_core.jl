#Usadel core
using SparseArrays
using LinearAlgebra
#Params struct
struct params
    E::Float64
    σn::Float64
    σs::Float64
    Γin::Float64 
    D::Vector{Float64}
    Δ::Vector{ComplexF64}
    dx::Float64
    N::Int
    Ln::Float64
    Ls::Float64    
    i0L::Int
    i0R::Int
    nodes::Dict{Int,Symbol}
end 

#setup
function setup_simulation_NS(Ln,Ls,dx,E,Γin,σn,σs,Dn,Ds)
    #just a bit cleaner
    nn = round(Int, Ln/dx)
    ns = round(Int, Ls/dx)
    N  = nn + ns
    i0L=nn
    i0R=nn+1
    #structures
    node_map=Dict{Int,Symbol}()
    D_tab=zeros(Float64,N)
    Δ_tab=zeros(ComplexF64,N)
    #inicialize structures
    for n in 1:N
        if n==1
            node_map[n]=:Nbulk
            D_tab[n]=Dn
            Δ_tab[n]=0.0
        elseif n<nn
            node_map[n]=:N
            D_tab[n]=Dn
            Δ_tab[n]=0.0
        elseif n==nn
            node_map[n]=:NS
            D_tab[n]=Dn
            Δ_tab[n]=0.0
        elseif n==nn+1
            node_map[n]=:SN
            D_tab[n]=Ds
            Δ_tab[n]=1.0
        elseif n>nn+1 && n<N
            node_map[n]=:S
            D_tab[n]=Ds
            Δ_tab[n]=1.0
        else
            node_map[n]=:bulk
            D_tab[n]=Ds
            Δ_tab[n]=1.0
        end
    end
    return params(E,σn,σs,Γin,D_tab,Δ_tab,dx,N,Ln,Ls,i0L,i0R,node_map)
end

function setup_simulation_SNS(Ls1,Ln,Ls2,dx,E,Γin,σn,σs,Dn,Ds)
    #just a bit cleaner
    ns1 = round(Int, Ls1/dx)
    nn  = round(Int, Ln/dx)
    ns2 = round(Int, Ls2/dx)
    N   = ns1 + nn + ns2
    i0L=ns1
    i0R=ns1+nn
    #structures
    node_map=Dict{Int,Symbol}()
    D_tab=zeros(Float64,N)
    Δ_tab=zeros(ComplexF64,N)
    #inicialize structures
    for n in 1:N
        if n==1
            node_map[n]=:Sbulk
            D_tab[n]=Ds
            Δ_tab[n]=1.0
        elseif n<ns1
            node_map[n]=:S
            D_tab[n]=Ds
            Δ_tab[n]=1.0
        elseif n==ns1
            node_map[n]=:SN
            D_tab[n]=Ds
            Δ_tab[n]=1.0
        elseif n==ns1+1
            node_map[n]=:NS
            D_tab[n]=Dn
            Δ_tab[n]=0.0
        elseif n>ns1+1 && n<ns1+nn+1
            node_map[n]=:N
            D_tab[n]=Dn
            Δ_tab[n]=0.0
        elseif n==ns1+nn+1
            node_map[n]=:NS
            D_tab[n]=Dn
            Δ_tab[n]=0.0 
        elseif n==ns1+nn+2 
            node_map[n]=:SN 
            D_tab[n]=Ds 
            Δ_tab[n]=1.0 
        elseif n>ns1+nn+2 && n<N 
            node_map[n]=:S 
            D_tab[n]=Ds 
            Δ_tab[n]=1.0 
        else 
            node_map[n]=:Sbulk 
            D_tab[n]=Ds 
            Δ_tab[n]=1.0 
        end 
     end 
     return params(E,σn,σs,Γin,D_tab,Δ_tab,dx,N,Ls1,Ls2,i0L,i0R,node_map)
end 

#index & DOS helpers
node_index(x::Real,p::params)=clamp((x>=0 ? p.i0R : p.i0L)+round(Int,x/p.dx),1,p.N)
dos_at_node(theta::AbstractVector{ComplexF64},idx::Int)=abs(real(cos(theta[idx])))

#initial guess
function get_theta_0(p::params)
    theta_0=zeros(ComplexF64,p.N)
    ω=-im*p.E+p.Γin
    atan_val=atan(p.Δ[end]/ω)
    for i in 1:p.N
        node=p.nodes[i]
        if node==:vacc || node==:N
            theta_0[i]=0.0
        elseif node==:S || node==:bulk
            theta_0[i]=atan_val
        elseif node==:NS || node==:SN
            theta_0[i]=atan_val/2
        elseif node==:Nbulk || node==:N
            theta_0[i]=0.0
        end
    end
    return theta_0
end

#build equation system
function build_eq_sys(theta,p::params,matsubara::Bool=false,ωn::Float64=0.0)
    #structure
    I,J,V=Int[],Int[],ComplexF64[]
    r=zeros(ComplexF64,p.N)
    h=p.dx
    h2=h^2
    ω=matsubara ? complex(ωn) : (-im*p.E+p.Γin)
    #helper functions
    function der_r(theta_i,i)
        return ω*cos(theta_i)+p.Δ[i]*sin(theta_i)
    end
    function get_rh(theta_i,i)
        return ω*sin(theta_i)-p.Δ[i]*cos(theta_i)
    end 
    for i in 1:p.N
        type=p.nodes[i]
        if type==:vacc
            #left boundary condition
            push!(I,i); push!(J,i); push!(V,1)
            push!(I,i); push!(J,i+1); push!(V,-1)
            r[i]=theta[i]-theta[i+1]
        elseif type==:N || type==:S
            #interior 
            push!(I,i); push!(J,i-1); push!(V,p.D[i]/(2*h2))
            push!(I,i); push!(J,i+1); push!(V,p.D[i]/(2*h2))
            d_r=der_r(theta[i],i)
            push!(I,i); push!(J,i); push!(V,(-p.D[i]/h2)-d_r)
            #r[i]=0
            d2th=(theta[i-1]-2*theta[i]+theta[i+1])/h2
            rh=get_rh(theta[i],i)
            r[i]=(p.D[i]/2)*d2th-rh
        elseif type==:SN
            # interface: sigma_N*(θ_NS - θ_Nprev)/h - sigma_S*(θ_Snext - θ_SN)/h = 0
            push!(I,i); push!(J,i-2); push!(V,-p.σn/h)
            push!(I,i); push!(J,i-1); push!(V, p.σn/h)
            push!(I,i); push!(J,i);   push!(V, p.σs/h)
            push!(I,i); push!(J,i+1); push!(V,-p.σs/h)
            r[i] = (p.σn*(theta[i-1] - theta[i-2]) - p.σs*(theta[i+1] - theta[i]))/h
        elseif type==:NS
            #interface
            push!(I,i); push!(J,i); push!(V,cos(theta[i]))
            push!(I,i); push!(J,i+1); push!(V,-cos(theta[i+1]))
            r[i]=sin(theta[i])-sin(theta[i+1])
        elseif type==:bulk
            #right boundary condition
            push!(I,i); push!(J,i); push!(V,1)
            thet=atan(p.Δ[i]/(-im*p.E+p.Γin)) #arctan(delta/omega)
            r[i]=theta[i]-thet
        elseif type==:Nbulk
            # left bulk normal reservoir: theta = 0
            push!(I,i); push!(J,i); push!(V,1)
            r[i] = theta[i]
        end
    end
    return sparse(I,J,V,p.N,p.N),r
end

#newton solver
function newton_basic(theta_0,p::params,max_iters::Int=50,tol::Real=1e-10,lambda::Real=0.5,matsubara::Bool=false,ωn::Float64=0.0)
    theta=copy(theta_0)
    for k in 1:max_iters
        J,r=build_eq_sys(theta,p,matsubara,ωn)
        if maximum(abs.(r))<=tol
            return theta, true, k
        end
        η=1e-8
        #I should have been identity
        J_reg = J + η*I
        dtheta=J_reg\(-r)
        theta.+=lambda.*dtheta
        if maximum(abs.(dtheta))<=tol #maybe dtheta*lambda??
            return theta, true, k
        end
    end
    return theta, false, max_iters
end

#DOS
function compute_DOS(energies::Vector{Float64},p::params,x::Real,maxIters::Int=50,tol::Real=1e-10,lambda::Real=0.5)
    idx=node_index(x,p)
    dos=zeros(Float64,length(energies))
    #slight change in order
    p0 = params(energies[1], p.σn, p.σs, p.Γin, p.D, p.Δ, p.dx, p.N, p.Ln, p.Ls, p.i0L, p.i0R, p.nodes)
    theta = get_theta_0(p0)
    for (k,E) in pairs(energies)
        p_E = params(E,p.σn,p.σs,p.Γin,p.D,p.Δ,p.dx,p.N,p.Ln,p.Ls,p.i0L,p.i0R,p.nodes)
        theta, converged, iters = newton_basic(theta, p_E, maxIters, tol, lambda)
        if !converged
            @warn "Did not converge for E=$E after $iters iterations"
        end
        dos[k] = dos_at_node(theta, idx)
    end
    return dos,idx
end

function self_consistent_delta(p::params,T::Float64,Tc::Float64,maxIters::Int=100,nMats::Int=50,tol::Real=1e-10,alpha::Float64=0.1)
    omegas=[(2n+1)*π*T for n in 0:nMats-1]
    Δ=copy(p.Δ)
    theta=get_theta_0(p)
    function update_delta(Δ,F,omegas)
        N=length(Δ)
        new_Δ=zeros(ComplexF64,N)
        for i in 1:N
            new_Δ[i]=(2π*T*sum(F[k,i] for k in 1:length(omegas)))
        end
        return new_Δ
    end
    for iter in 1:maxIters
        F=zeros(ComplexF64,nMats,p.N)
        for (k,ωn) in pairs(omegas)
            p_k=params(0.0,p.σn,p.σs,p.Γin,p.D,Δ,p.dx,p.N,p.Ln,p.Ls,p.i0L,p.i0R,p.nodes)
            theta,_,_=newton_basic(theta,p_k,maxIters,tol,lambda=0.5,matsubara=true,ωn=ωn)
            F[k,:].=sin.(theta)
        end
        new_Δ=update_delta(Δ,F,omegas)
        #enforce physics
        new_Δ[1:p.i0L].=0.0
        new_Δ[p.N]=p.Δ[p.N]
        mixed_Δ=(1-alpha).*Δ.+alpha.*real(new_Δ)
        mixed_Δ.=clamp(mixed_Δ,-5.0,5.0)

        diff=maximum(abs.(mixed_Δ.-Δ))
        Δ.=mixed_Δ

        if diff<tol
            println("Converged after $iter iterations with max Δ change of $diff")
            return Δ,theta
        end
    end
    @warn "Did not converge after $maxIters iterations, final max Δ change: $diff"
    return Δ,theta
end