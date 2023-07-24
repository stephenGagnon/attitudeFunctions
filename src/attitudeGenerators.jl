"""
    Generates a set of random attitudes

    inputs:
    N - number of attitudes to generate
    T - type of output
        valid types:
        quaternion
        MRP
        GRP
        DCM
    customTypes:
        if true, function will return an Nx1 array of custom attitude parameter
            types as defined in this package
        if false, function will return an array of size MxN where M is the dimension
            of the parameterization (4 for quaternion, 3x3 for DCM etc.) and N is
            the requested number of attitudes
    a,f:
        parameters that specify the exact GRP tranformation. Default to a=f=1
        which corresponds to an MRP
"""
function randomAtt(N :: Int64, T=MRP, a = 1, f = 1; vectorize = false, customTypes = false, randomType = :LHS)

    if randomType == :LHS
        val = lhs(3,N)
        val[2:3,:] = val[2:3,:].*2*pi;
    elseif randomType == :uniform
        val = rand(3,N)
        val[2:3,:] = val[2:3,:].*2*pi;
    end

    if N != 1
        # q = zeros(4,N)
        # q[1,:] = sqrt.(val[1,:]).*cos.(val[2,:])
        # q[2,:] = sqrt.(val[1,:]).*sin.(val[2,:])
        # q[3,:] = sqrt.(1 .- val[1,:]).*sin.(val[3,:])
        # q[4,:] = sqrt.(1 .- val[1,:]).*cos.(val[3,:])
        q = [randomQGen.(eachcol(val))...]

        if vectorize & customTypes
            error("Both vectorize and custom types cannot be selected")
        end

        if vectorize
            q = hcat(q...)
        end

        if customTypes
            q = quaternion.(q)
        end
    else
        # q = Array{Float64,1}(undef,4)
        # q[1] = sqrt(val[1])*cos(val[2])
        # q[2] = sqrt(val[1])*sin(val[2])
        # q[3] = sqrt(1 - val[1])*sin(val[3])
        # q[4] = sqrt(1 - val[1])*cos(val[3])
        q = randomQGen(val)
    end

    if T == quaternion
        return q
    elseif T == MRP
        return q2p(q)
    elseif T == GRP
        return q2p(q,a,f)
    elseif T == DCM
        return q2A(q)
    else
        throw(error("Please provide a valid attitude parameterization"))
    end
end

# parallelized version
function randomAttPar(N :: Int64, T=MRP)

    x = rand(3,N)

    Threads.@threads for i = 1:N
        x[1,i] = x[1,i]*1/N + (i-1)/N
        x[2,i] = (x[2,i]*1/N + (i-1)/N)*2*pi
        x[3,i] = (x[3,i]*1/N + (i-1)/N)*2*pi
    end

    for i = 1:3
        x[i,:]  = x[i,randperm(N)]
    end
    # transpose(hcat([row[randperm(d)] for row in eachrow(x)]...))

    if T == quaternion
        out = Matrix{Float64}(undef,4,N)
    elseif T == MRP
        out = Matrix{Float64}(undef,3,N)
    end

    Threads.@threads for i = 1:N
        if T == quaternion
            out[1,i] = sqrt(x[1,i])*cos(x[2,i])
            out[2,i] = sqrt(x[1,i])*sin(x[2,i])
            out[3,i] = sqrt(1 - x[1,i])*sin(x[3,i])
            out[4,i] = sqrt(1- x[1,i])*cos(x[3,i])
        elseif T == MRP
            q = Vector{Float64}(undef,4)
            q[1] = sqrt(x[1,i])*cos(x[2,i])
            q[2] = sqrt(x[1,i])*sin(x[2,i])
            q[3] = sqrt(1 - x[1,i])*sin(x[3,i])
            q[4] = sqrt(1 - x[1,i])*cos(x[3,i])
            out[:,i] = q2p(q)
        end
    end

    return out
end

function randomBoundedAngularVelocity(N :: Int64, upperBound :: Number, vectorize = false)
    w = lhs(3,N)
    # wn_max = 0;
    # for i = 1:N
    #     wn = norm(w[:,i])
    #     if wn > wn_max
    #         wn_max = wn
    #     end
    # end
    w = (w .- .5).*(2*upperBound)

    if !vectorize
        if N != 1
            w = [w[:,i] for i in 1:size(w,2)]
        else
            w = w[:,1]
        end
    end

    return w
end

function randomAttState(N :: Int64, upperBound :: Number, T=MRP, a = 1, f = 1; randomType = :LHS)

    atts = randomAtt(N , T, a, f, randomType = randomType)
    ws = randomBoundedAngularVelocity(N, upperBound)
    if N > 1
        Out = Array{Array{Float64,1},1}(undef, N)
        for i = 1:N
            Out[i] = vcat(atts[i], ws[i])
        end
    elseif N == 1
        Out = vcat(atts, ws)
    end
    return Out
end

"""
    latin hypercube sampling:
    generates d samples with N dimensions
    samples have N elements with values in the range 0-1
    there is guarenteed to be one sample with a value in the
    ranges 0-1/d, 1/d-2/d, ... , (1-d)/d-1 for each dimension
"""
function lhs(N :: Int64, d :: Int64)

    x = rand(N,d)

    x = x.*1/d .+ reshape((collect(0:d-1))./d,1,d)

    return transpose(hcat([row[randperm(d)] for row in eachrow(x)]...))
end

function randomQGen(v)

    q = Array{typeof(v[1]),1}(undef,4)
    q[1] = sqrt(v[1])*cos(v[2])
    q[2] = sqrt(v[1])*sin(v[2])
    q[3] = sqrt(1 - v[1])*sin(v[3])
    q[4] = sqrt(1 - v[1])*cos(v[3])

    return q
end

"""
    Function to randomly perturb attitudes by a small amount.
    inputs:
    atts/att -- an attitude or 1D array of attitudes (see the types file for details on supported attitude reresentations (anyAttitude, arrayofAtts))
    p -- Controls the size of perturbation.
    outputs:
    out -- an attitude or 1D array of attitudes. There is one output attitude for each input, which is slightly perturbed (a small rotation is applied)
"""
function attitudeRoughening(atts :: arrayofAtts, p = .001)

    out = similar(atts)
    for i = eachindex(atts)
        out[i] = attitudeRoughening(atts[i],p)
    end
    return out
end

function attitudeRoughening(att :: quaternion, p = .001)
    qp = randomQGen(rand(3) .* p .* [1;2*pi;2*pi])
    return qprod(att,qp)
end

function attitudeRoughening(att :: MRP, p = .001)
    qp = randomQGen(rand(3) .* p .* [1;2*pi;2*pi])
    return q2p(qprod(p2q(att),qp))
end

function attitudeRoughening(att :: DCM, p = .001)
    qp = randomQGen(rand(3) .* p .* [1;2*pi;2*pi])
    return q2A(qprod(A2q(att),qp))
end

function attitudeRoughening(att :: Mat, p = .001)
    qp = randomQGen(rand(3) .* p .* [1;2*pi;2*pi])
    return q2A(qprod(A2q(att),qp))
end

function attitudeRoughening(att :: Vec, p = .001)
    qp = randomQGen(rand(3) .* p .* [1;2*pi;2*pi])
    if length(att) == 3
        return q2p(qprod(p2q(att),qp))
    elseif length(att) == 4
        return qprod(att,qp)
    end
end