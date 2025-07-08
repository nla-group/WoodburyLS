
function WoodburyLS(A, b, U, V, x0, AtAsolver)
    # Follow-up solve using the Woodbury identity
    r = size(U, 2)
    n = size(A,2)

    # X = [V  A' * U]
    X = Matrix{Float64}(undef, n, 2r)
    @views X[:, 1:r] .= V
    @views mul!(X[:, r+1:2r], A', U)

    # Yt = [X[:, r+1:2r]' + (U' * U) * V'; V']
    Yt = Matrix{Float64}(undef, 2r, n)
    @views Yt[1:r, :] .= X[:, r+1:2r]'
    @views mul!(Yt[1:r,:], (U' * U), V', 1.0, 1.0)
    @views Yt[r+1:end, :] .= V'
    
    Z = AtAsolver(X)

    # x = x0 + Z[:, 1:r] * (U' * b)
    x = copy(x0)
    @views mul!(x, Z[:, 1:r], U' * b, 1.0, 1.0)

    M = I + Yt * Z

    # x = x - Z * (M \ (Yt * x))
    mul!(x, Z, M \ (Yt * x), -1.0, 1.0)
    return x
end

function WoodburyLS_setup(A, b)
    # Initial solve: return x0 = A \ b and an in-place solver for (A'A)x = b
    F = qr(A)
    x = F \ b
    RtR = Cholesky(UpperTriangular(F.R))
    AtAsolver = X -> ldiv!(RtR, X)
    return x, AtAsolver
end
