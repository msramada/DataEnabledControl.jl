using DataEnabledControl
using Test

@testset "DataEnabledControl.jl" begin
    p = m = 2
    Q = R = [1.0 0;0 1.0]
    Tᵢₙᵢ = 5
    N = 10
    T = 60
    U = zeros(p, T); Y = zeros(m, T)
    Constrained = true
    prob = DeePC_struct(Q, R, N, U, Y, Tᵢₙᵢ,Constrained)
    @test prob.constraintMatrices.Ay isa Nothing
    @test size(prob.PastDataMatrix)[1] == Tᵢₙᵢ*m + Tᵢₙᵢ*p
    @test DataEnabledControl.hankellize([1.0 2.0;], 2) == [1.0; 2.0 ;;]
    Uv = vec(zeros(m, Tᵢₙᵢ))
    Yv = vec(zeros(p, Tᵢₙᵢ))
    uₜ, yₜ = run_DeePC(prob, Uv, Yv, slack_var = false, λy=100.0, λg=1.0)
    @test uₜ == zeros(m, N)
    @test yₜ == zeros(p, N)
    
end
