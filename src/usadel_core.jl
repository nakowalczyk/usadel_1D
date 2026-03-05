#Usadel core

#Definicja parametrów
struct params
    E::Float64
    Γin::Float64 
    D::Vector{Float64}
    Δ::Vector{Float64}
    dx::Float64
    N::Int
end 

#Parametry symulacji
Ln, Ls =200,100
dx=1.0
nn=Int(Ln/dx)
ns=Int(Ls/dx)
N=Int((Ln+Ls)/dx)

nodes=Dict{Int,Dict{Symbol,Any}}()
#Inicjalizacja słownika
for n in 1:N
    nodes[n]=Dict()
    if n==1
        nodes[n][:type]=:vacc
    elseif n<nn
        nodes[n][:type]=:N
    elseif n==nn
        nodes[n][:type]=:NS
    elseif n==nn+1
        nodes[n][:type]=:SN
    elseif n>nn+1 && n<N
        nodes[n][:type]=:S
    else
        nodes[n][:type]=:vacc
    end
end