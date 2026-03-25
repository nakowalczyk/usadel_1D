#Usadel core
#Definicja parametrów
struct params
    E::Float64
    σn::Float64
    σs::Float64
    Γin::Float64 
    D::Vector{Float64}
    Δ::Vector{ComplexF64}
    dx::Float64
    N::Int
    nodes::Dict{Int,Symbol}
end 

#Parametry symulacji
Ln, Ls =200,100
dx=1.0

function setup_simulation(Ln,Ls,dx,E,Γin)
    nn=Int(Ln/dx)
    ns=Int(Ls/dx)
    N=Int((Ln+Ls)/dx)
    #struktury
    node_map=Dict{Int,Symbol}()
    D_tab=zeros(N)
    Δ_tab=zeros(ComplexF64,N)
    #inicjalizacja słownika
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
    return params(E,Γin,D_tab,Δ_tab,dx,N,node_map)
end


#budowa macierzy rzadkiej
function build_eq_sys(theta,p::params)
    #struktura
    I,J,V=Int[],Int[],ComplexF64[]
    r=zeros(ComplexF64,p.N)
    h=p.dx
    h2=h^2
    #pochodne
    function der_r(theta_i,i)
        term_E=(-im*p.E+p.Γin)
        return term_E*cos(theta_i)+p.Δ[i]*sin(theta_i)
    end
    function get_rh(theta_i,i)
        term_E=(-im*p.E+p.Γin)
        return term_E*sin(theta_i)-p.Δ[i]*cos(theta_i)
    end 

    for i in 1:p.N
        if type==:vacc
            #warunek brzegowy 
            push!(I,i); push!(J,i); push!(V,1)
            push!(I,i); push!(J,i-1); push!(V,-1)
            r[i]=theta[i]-theta[i-1]
        elseif type==:N || type==:S
            #interior 
            push!(I,i); push!(J,i-1); push!(V,D[i]/2*h2)
            push!(I,i); push!(J,i+1); push!(V,D[i]/2*h2)
            d_r=der_r(theta[i],i)
            push!(I,i); push!(J,i); push!(V,(D[i]/2*h2)-d_r)
            #r[i]=0
            d2th=(theta[i-1]-2*theta[i]+theta[i+1])/h2
            rh=get_rh(theta[i],i)
            r[i]=(p.D[i]/2)*d2th-rh
        elseif type==:NS
            #interface
            push!(I,i); push!(J,i); push!(V,σn/h)
            push!(I,i); push!(J,i-1); push!(V,σs/h)
            push!(I,i); push!(J,i+1); push!(V,-σn/h)
            push!(I,i); push!(J,i-2); push!(V,-σs/h)
            r[i]=(σs*theta[i-1]-σs*theta[i-2]-σn*theta[i+1]+σn*theta[i])/h
        elseif type==:SN
            #interface
            push!(I,i); push!(J,i); push!(V,-cos(theta[i]))
            push!(I,i); push!(J,i+1); push!(V,cos(theta[i+1]))
            r[i]=sin(theta[i]-sin(theta[i+1]))
        elseif type==:bulk
            #warunek brzegowy
            push!(I,i); push!(J,i); push!(V,1)
            thet=atan(p.Δ[i]/(-im*p.E+p.Γin)) #arctan(delta/omega)
            r[i]=theta[i]-thet
        end
    end
    return sprase(I,J,V,p.N,p.N),r
end
