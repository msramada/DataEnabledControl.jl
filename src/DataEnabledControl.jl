
"""
DataEnabledControl.jl is a Julia package implementing the DeePC algorithm in the paper 
"Coulson, J., Lygeros, J., & Dörfler, F. (2019, June). Data-enabled predictive control: 
In the shallows of the DeePC. In 2019 18th European Control Conference (ECC) (pp. 307-312). IEEE."
"""
module DataEnabledControl

# ==== imports ====
using Convex, SCS

# ==== exports ====
export run_DeePC, DeePC_struct, ConstraintMatrices

"""
    ConstraintMatrices
Stores the constraint matrices representing polyhedron constraints of the form:
`Ay * y <= by` and `Au * u <= bu`.
"""
mutable struct ConstraintMatrices
    Ay::Union{Matrix{Float64}, Nothing}
    by::Union{Vector{Float64}, Nothing}
    Au::Union{Matrix{Float64}, Nothing}
    bu::Union{Vector{Float64}, Nothing}
end
"""
    DeePC_struct
Used to formulate the DeePC optimization problem
Q:                  the weighting matrix of the output signal
R:                  the weighting matrix of the input signal
PastDataMatrix:     the left tall matrix [Uₚ; Yₚ] in the paper
FutureInputMatrix:  Uf according to the paper
FutureOutputMatrix: Yf according to the paper
Constrained:        true if the system is input or output constrained
constraintMatrices: params of polyhedral constraints
"""
mutable struct DeePC_struct
    # Data and parameters required to run DeePC
    Q::Matrix{Float64}
    R::Matrix{Float64}
    N::Int
    PastDataMatrix::Matrix{Float64} # [Uₚ; Yₚ]
    FutureInputMatrix::Matrix{Float64}
    FutureOutputMatrix::Matrix{Float64}
    Constrained::Bool
    constraintMatrices::ConstraintMatrices

    function DeePC_struct(Q::Matrix{Float64}, R::Matrix{Float64}, 
                            N::Int, U::Matrix{Float64}, Y::Matrix{Float64},
                            Tᵢₙᵢ::Int, Constrained::Bool)
        constraintMatrices = ConstraintMatrices(nothing, nothing, nothing, nothing)
        PastDataMatrix, FutureInputMatrix, FutureOutputMatrix = build_data_matrix(U, Y, Tᵢₙᵢ, N)
            m, T = size(U)
            if T< (m+1) * (Tᵢₙᵢ + N + 1) - 1
                error("Length of training data (T) needs to be bigger to at least detect a 1st order module.")
            end
        new(Q, R, N, PastDataMatrix, FutureInputMatrix, FutureOutputMatrix, Constrained, constraintMatrices)        
    end
end

"""
    build_data_matrix(U::Matrix{Float64}, Y::Matrix{Float64}, Tᵢₙᵢ::Int, N::Int)
Takes the data collected `U,Y` through a persistently excited experiment and returns
the matrices `[Up; Yp], Uf, Yf` used in the construction of the DeePC algorithm.
"""
function build_data_matrix(U::Matrix{Float64}, Y::Matrix{Float64}, Tᵢₙᵢ::Int, N::Int)
    Hu = hankellize(U, Tᵢₙᵢ+N)
    Hy = hankellize(Y, Tᵢₙᵢ+N)

    m = size(U)[1]
    p = size(Y)[1]

    Up = Hu[1:Tᵢₙᵢ*m,:];      Uf = Hu[Tᵢₙᵢ*m+1:end,:];  
    Yp = Hy[1:Tᵢₙᵢ*p,:];      Yf = Hy[Tᵢₙᵢ*p+1:end,:];

    return [Up; Yp], Uf, Yf
end

"""
    hankellize(U::Matrix{Float64}, L::Int) 
Constructs the hankel matrix `H`` of a matrix `U`, with `L` number of shifts. `H` here corresponds to `H_L(U)` in the DeePC paper.
"""
function hankellize(U::Matrix{Float64}, L::Int)
    m, T = size(U)

    if T<L
        error("T is less than L, cannot construct a Hankel matrix.")
    end
    H = zeros(m*L, T-L+1)
    for row=1:L
        H[(row-1)*m+1:row*m,:] = U[:,row:T-L+row]
    end
    return H
end

"""
    run_DeePC(prob::DeePC_struct, Uᵢₙᵢ::Vector{Float64}, Yᵢₙᵢ::Vector{Float64}; slack_var = false, λy = 1.0, λg = 0.0)
Constructs and solves the DeePC optimization problem, with a slack variable or without (noisy/nonlinear or not).
It returns `uₜ` and `yₜ`, the input and output sequence solutions of the DeePC over the N-horizon.
The regularization weights (the λ's) are set to zero and the slack variable is deactivated in the default settings.
"""
function run_DeePC(prob::DeePC_struct, Uᵢₙᵢ::Vector{Float64}, Yᵢₙᵢ::Vector{Float64}; 
    slack_var = false, λy = 1.0, λg = 0.0)


if slack_var == false
    λy = 1.0; λg = 0.0; # Some weights with no effect on optimization
end



p = size(prob.R)[1]
m = size(prob.Q)[1]

uₜ = Variable(p, prob.N)
yₜ = Variable(m, prob.N)

g = Variable(size(prob.PastDataMatrix)[2])
σ = Variable(size(Yᵢₙᵢ)[1])

# Define objective function:
cost =  λg * norm(g,1) + λy * norm(σ,1) # this part will have an effect if slack variable is true
for i in 1:prob.N
    cost += quadform(uₜ[:,i], prob.R; assume_psd=true) + quadform(yₜ[:,i], prob.Q; assume_psd=true)
end

# Define optimization problem
problem = minimize(cost)

# Define dynamic constraints (data-driven representation of the dynamics):
problem.constraints += [uₜ[:,i] == prob.FutureInputMatrix[(i-1)*p+1:i*p,:] * g for i in 1:prob.N]
problem.constraints += [yₜ[:,i] == prob.FutureOutputMatrix[(i-1)*m+1:i*m,:] * g for i in 1:prob.N]

# Include the slack variable for some flexibility when facing noise and nonlinearities:
if slack_var == false
    problem.constraints += [prob.PastDataMatrix * g == [Uᵢₙᵢ; Yᵢₙᵢ]]
elseif slack_var == true
    problem.constraints += [prob.PastDataMatrix * g == [Uᵢₙᵢ; Yᵢₙᵢ + σ]]
end

# Include the input and output constraints (if they exist)
if prob.Constrained == true
    Mats = prob.constraintMatrices
    if !(Mats.Au isa Nothing || Mats.bu isa Nothing)
        problem.constraints += [Mat.Au * uₜ[:,i] <= Mats.bu for i in 1:prob.N]
    end
    if !(Mats.Ay isa Nothing || Mats.by isa Nothing)
        problem.constraints += [Mats.Ay * yₜ[:,i] <= Mats.by for i in 1:prob.N]
    end
end

solve!(problem, SCS.Optimizer; silent_solver = true)

return uₜ.value, yₜ.value
end

end
