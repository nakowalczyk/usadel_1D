#Usadel core
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
    i0L::Int
    i0R::Int
    Ln::Float64
    Ls::Float64
    nodes::Dict{Int,Symbol}
end 

#setup
function setup_simulation(Ln,Ls,dx,E,Γin,σn,σs)
    nn=Int(Ln/dx)
    i0L=nn
    ns=Int(Ls/dx)
    i0R=nn+1
    N=Int((Ln+Ls)/dx)
    #structures
    node_map=Dict{Int,Symbol}()
    D_tab=zeros(Float64,N)
    Δ_tab=zeros(ComplexF64,N)
    #inicialize structures
    for n in 1:N
        if n==1
            node_map[n]=:vacc
            D_tab[n]=0.0152
            Δ_tab[n]=0.0
        elseif n<nn
            node_map[n]=:N
            D_tab[n]=0.0152
            Δ_tab[n]=0.0
        elseif n==nn
            node_map[n]=:NS
            D_tab[n]=0.0152
            Δ_tab[n]=0.0
        elseif n==nn+1
            node_map[n]=:SN
            D_tab[n]=0.00636
            Δ_tab[n]=1.0
        elseif n>nn+1 && n<N
            node_map[n]=:S
            D_tab[n]=0.00636
            Δ_tab[n]=1.0
        else
            node_map[n]=:bulk
            D_tab[n]=0.00636
            Δ_tab[n]=1.0
        end
    end
    return params(E,σn,σs,Γin,D_tab,Δ_tab,dx,N,Ln,Ls,i0L,i0R,node_map)
end

#index & DOS helpers
node_index(x::Real,p::params)=clamp((x>=0 ? p.i0R : p.i0L)+round(Int,x/p.dx),1,p.N)
dos_at_node(theta::AbstractVector{ComplexF64},idx::Int)=abs(real(cos(theta[idx])))

#initial guess
function get_theta_0(p::params)
    theta_0=zeros(ComplexF64,p.N)
    ω=-im*p.E+p.Γin
    atan_val=atan(p.Δ[1]/ω)
    for i in 1:p.N
        node=p.nodes[i]
        if node==:vacc || node==:N
            theta_0[i]=0.0
        elseif node==:S || node==:bulk
            theta_0[i]=atan_val
        elseif node==:NS || node==:SN
            theta_0[i]=atan_val/2
        end
    end
    return theta_0
end

#build equation system
function build_eq_sys(theta,p::params)
    #structure
    I,J,V=Int[],Int[],ComplexF64[]
    r=zeros(ComplexF64,p.N)
    h=p.dx
    h2=h^2
    #helper functions
    function der_r(theta_i,i)
        ω=(-im*p.E+p.Γin)
        return ω*cos(theta_i)+p.Δ[i]*sin(theta_i)
    end
    function get_rh(theta_i,i)
        ω=(-im*p.E+p.Γin)
        return ω*sin(theta_i)-p.Δ[i]*cos(theta_i)
    end 
    for i in 1:p.N
        if type==:vacc
            #left boundary condition
            push!(I,i); push!(J,i); push!(V,1)
            push!(I,i); push!(J,i+1); push!(V,-1). 
            r[i]=theta[i]-theta[i+1]
        elseif type==:N || type==:S
            #interior 
            push!(I,i); push!(J,i-1); push!(V,p.D[i]/2*h2)
            push!(I,i); push!(J,i+1); push!(V,p.D[i]/2*h2)
            d_r=der_r(theta[i],i)
            push!(I,i); push!(J,i); push!(V,(p.D[i]/2*h2)-d_r)
            #r[i]=0
            d2th=(theta[i-1]-2*theta[i]+theta[i+1])/h2
            rh=get_rh(theta[i],i)
            r[i]=(p.D[i]/2)*d2th-rh
        elseif type==:NS
            #interface
            push!(I,i); push!(J,i); push!(V,p.σn/h)
            push!(I,i); push!(J,i-1); push!(V,p.σs/h)
            push!(I,i); push!(J,i+1); push!(V,-p.σn/h)
            push!(I,i); push!(J,i-2); push!(V,-p.σs/h)
            r[i]=(p.σs*theta[i-1]-p.σs*theta[i-2]-p.σn*theta[i+1]+p.σn*theta[i])/h
        elseif type==:SN
            #interface
            push!(I,i); push!(J,i); push!(V,-cos(theta[i]))
            push!(I,i); push!(J,i+1); push!(V,cos(theta[i+1]))
            r[i]=sin(theta[i]-sin(theta[i+1]))
        elseif type==:bulk
            #right boundary condition
            push!(I,i); push!(J,i); push!(V,1)
            thet=atan(p.Δ[i]/(-im*p.E+p.Γin)) #arctan(delta/omega)
            r[i]=theta[i]-thet
        end
    end
    return sparse(I,J,V,p.N,p.N),r
end

#newton solver
function newton_basic(theta_0,p::params,max_iters::Int=50,tol::Real=1e-10,lambda::Real=0.5)
    theta=copy(theta_0)
    for k in 1:max_iters
        J,r=build_eq_sys(theta,p)
        if maximum(abs.(r))<=tol
            return theta
        end
        η=1e-8
        J_reg=J+η*I
        dtheta=J_reg\(-r)
        theta.+=lambda.*dtheta
        if maximum(abs.(dtheta))<=tol #maybe dtheta*lambda??
            return theta
        end
    end
    return theta
end

#DOS
function compute_DOS(energies::Vector{Float64},p::params,x::Real,maxIters::Int=50,tol::Real=1e-10,lambda::Real=0.5)
    n=p.N
    idx=node_index(x,p)
    theta=get_theta_0(p)
    dos=zeros(Float64,length(energies))
    for (k,E) in pairs(energies)
        p_E=params(E,p.σn,p.σs,p.Γin,p.D,p.Δ,p.dx,p.N,p.Ln,p.Ls,p.nodes)
        theta,converged,iters=newton_basic(theta,p_E,maxIters,tol,lambda)
        dos[k]=dos_at_node(theta,idx)
    end
    return dos,idx
end
