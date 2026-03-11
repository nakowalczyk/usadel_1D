#Usadel core
#inicjalizacja słownika i tablic mozna zrobić w jednej pętli, ale dla przejrzystości zostawiłam to narazie w dwóch osobnych
#Definicja parametrów
struct params
    E::Float64
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
nn=Int(Ln/dx)
ns=Int(Ls/dx)
N=Int((Ln+Ls)/dx)

#Inicjalizacja słownika
node_map=Dict{Int,Symbol}()
for n in 1:N
    if n==1
        node_map[n]=:vacc
    elseif n<nn
        node_map[n]=:N
    elseif n==nn
        node_map[n]=:NS
    elseif n==nn+1
        node_map[n]=:SN
    elseif n>nn+1 && n<N
        node_map[n]=:S
    else
        node_map[n]=:bulk
    end
end

#Inicjalizacja tablic
D_tab=zeros(N)
Δ_tab=zeros(ComplexF64,N)

for n in 1:N
    type=node_map[n]
    if type==:vacc || type==:N || type==:NS
        D_tab[n]=0.0152
        Δ_tab[n]=0.0
    elseif type==:SN || type==:S || type==:bulk
        D_tab[n]=0.00636
        Δ_tab[n]=1.0
    end
end

#sprawdzenie
println(params(0.8,12-3,D_tab,Δ_tab,dx,N,node_map))
println(node_map[256])
println(D_tab[256],Δ_tab[256])