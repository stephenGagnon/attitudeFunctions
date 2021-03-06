"""
    Computes the error between two quaternions by finding the rotation between them
    and returning 2x the vector part of the quaternion. This corresponds to the roll
    pitch yaw errors for small error angles.

        !!!! Assumes inputs are quaternions !!!!
"""
function attitudeErrors(qE :: Vec, q :: Vec)

    if !(sign(qE[4]) == sign(q[4]) == -1) | !(sign(qE[4]) == sign(q[4]) == 1)
            qE = -qE
    end
    return 2*qprod(qE,qinv(q))[1:3]
end

function attitudeErrors(qE :: Union{Mat,Vec}, q :: Union{Mat,Vec})

    s1 = size(qE)
    s2 = size(q)

    if s1[1] < 4 | s1[1] > 4
        qE = qE'
    end

    if s2[1] < 4 | s2[1] > 4
        q = q'
    end

    if (length(s1) == 1) & (length(s2) > 1)
        dalpha = zeros(3,s2[2])
        for i = 1:size(q,2)
            dalpha[:,i] = attitudeErrors(qE,q[:,i])
        end
    elseif (length(s1) > 1) & (length(s2) == 1)
        dalpha = zeros(3,s1[2])
        for i = 1:size(qE,2)
            dalpha[:,i] = attitudeErrors(qE[:,i],q[:])
        end
    elseif (length(s1) == length(s2) > 1)
        dalpha = zeros(3,s2[2])
        for i = 1:size(q,2)
            dalpha[:,i] = attitudeErrors(qE[:,i],q[:,i])
        end
    else
        error("quaternion arrays must be the same size or a single quaternion")
    end

    return dalpha
end

function attitudeErrors(qE :: Union{Vecs,Vec}, q :: Union{Vecs,Vec})

    if (typeof(qE) <: Vec) & (typeof(q) <: Vecs)
        dalpha = Array{typeof(qE),1}(undef,length(q))
        for i = 1:length(q)
            dalpha[i] = attitudeErrors(qE,q[i])
        end
    elseif (typeof(q) <: Vec) & (typeof(qE) <: Vecs)
        dalpha = Array{typeof(q),1}(undef,length(qE))
        for i = 1:length(qE)
            dalpha[i] = attitudeErrors(qE[i],q)
        end
    elseif (typeof(q) <: Vecs) & (typeof(qE) <: Vecs)
        dalpha = Array{typeof(q[1]),1}(undef,length(q))
        for i = 1:length(qE)
            dalpha[i] = attitudeErrors(qE[i],q[i])
        end
    else
        error("quaternion arrays must be the same size or a single quaternion")
    end

    return dalpha
end

function attitudeErrors(qE :: quaternion, q :: quaternion)

    if !(sign(qE.s) == sign(q.s) == -1) | !(sign(qE.s) == sign(q.s) == 1)
            qE = -qE
    end
    return 2*qprod(qE,qinv(q))[1:3]
end

function attitudeErrors(qE :: Union{Array{quaternion, 1}, quaternion},
                        q :: Union{Array{quaternion, 1}, quaternion})

    if typeof(qE) == quaternion
        dalpha = Array{quaternion,1}(undef,length(q))
        for i = 1:size(q,2)
            dalpha[i] = attitudeErrors(qE,q[i])
        end
    elseif typeof(qE) == quaternion
        dalpha = Array{quaternion,1}(undef,length(qE))
        for i = 1:size(qE,2)
            dalpha[i] = attitudeErrors(qE[i],q)
        end
    else
        dalpha = Array{quaternion,1}(undef,length(q))
        for i = 1:size(q)
            dalpha[i] = attitudeErrors(qE[i],q[i])
        end
    end

    return dalpha
end

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

function randomBoundedAngularVelocity(N :: Int64, upperBound :: Float64, vectorize = false)
    w = lhs(3,N)
    # wn_max = 0;
    # for i = 1:N
    #     wn = norm(w[:,i])
    #     if wn > wn_max
    #         wn_max = wn
    #     end
    # end
    w = w.*upperBound

    if !vectorize
        if N != 1
            w = [w[:,i] for i in 1:size(w,2)]
        else
            w = w[:,1]
        end
    end

    return w
end

function randomAttState(N :: Int64, upperBound :: Float64, T=MRP, a = 1, f = 1; randomType = :LHS)

    atts = randomAtt(N , T, a, f, randomType = randomType)
    ws = randomBoundedAngularVelocity(N, upperBound)
    Out = Array{Array{Float64,1},1}(undef,N)
    for i = 1:N
        Out[i] = vcat(atts[i],ws[i])
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

function attitude2Array(x :: Union{Array{MRP,1},Array{GRP,1},Array{quaternion,1}})

    if (typeof(x[1]) == MRP) | (typeof(x[1]) == GRP)
        out = Array{Float64,2}(undef,3,length(x))
    elseif typeof(x[1]) == quaternion
        out = Array{Float64,2}(undef,4,length(x))
    end

    for i = 1:length(x)
        if (typeof(x[1]) == MRP) | (typeof(x[1]) == GRP)
            out[:,i] = x[i].p
        elseif typeof(x[1]) == quaternion
            out[1:3,i] = x[i].v
            out[4,i] = x[i].s
        end
    end
    return out
end

function attitude2Vecs(x :: Union{Array{MRP,1},Array{GRP,1},Array{quaternion,1}})

    out = Array{Array{Float64,1}}(undef,length(x))

    for i = 1:length(x)
        if (typeof(x[1]) == MRP) | (typeof(x[1]) == GRP)
            out[i] = x[i].p
        elseif typeof(x[1]) == quaternion
            out[i][1:3] = x[i].v
            out[i][4] = x[i].s
        end
    end
    return out
end

function toBodyFrame(att :: anyAttitude, usun :: Vec, uobs :: MatOrVecs, a = 1, f = 1)

    if (typeof(att) <: Vec) & (length(att) == 3)
        rotFunc = ((A,v) -> p2A(A,a,f)*v)
    elseif ((typeof(att) <: Vec) & (length(att) == 4)) | (typeof(att) == quaternion)
        rotFunc = qRotate
    elseif (typeof(att) <: Mat) & (size(att) == (3,3))
        rotFunc = ((A,v) -> A*v)
    elseif typeof(att) <: Union{DCM,MRP,GRP}
        rotFunc = ((A,v) -> any2A(A).A*v)
    else
        error("Please provide a valid attitude. Attitudes must be represented
        as a single 3x1 or 4x1 float array, a 3x3 float array, or any of the
        custom attitude types defined in the attitueFunctions package.")
    end
    return _toBodyFrame(att,usun,uobs,rotFunc)
end

function _toBodyFrame(att :: anyAttitude{T}, usun :: Vec, uobs :: ArrayOfVecs, rotFunc :: Function) where {T <: Real}
    #
    usunb = rotFunc(att,usun)

    uobsb = Array{Array{T,1},1}(undef,length(uobs))

    for i = 1:length(uobs)
        uobsb[i] = rotFunc(att,uobs[i])
    end
    # uobsb = map(x -> rotFunc(att,x), uobs)

    return usunb :: Vec, uobsb :: ArrayOfVecs
end

function _toBodyFrame!(att :: anyAttitude, usun :: Vec, uobs :: ArrayOfVecs, rotFunc :: Function, usunb, uobsb)
    rotFunc(att, usun, usunb)

    for i = 1:length(uobs)
        rotFunc(att, uobs[i], uobsb[i])
    end
end

function _toBodyFrame(att :: anyAttitude{T}, usun :: Vec, uobs :: Mat, rotFunc :: Function) where {T <: Real}

    usunb = rotFunc(att,view(usun,:))

    uobsb = Array{T,2}(undef,size(uobs))

    for i = 1:size(uobs,2)
        uobsb[:,i] = rotFunc(att,view(uobs,:,i))
    end

    return usunb :: Vec, uobsb :: ArrayOfVecs
end

function _toInertialFrame(att :: anyAttitude{T}, un :: Mat, uu :: Mat, uv :: Mat, rotFunc :: Function, parameterization) where {T <: Real}

    unb = Array{T,2}(undef,size(un))
    uub = Array{T,2}(undef,size(uu))
    uvb = Array{T,2}(undef,size(uv))

    if parameterization == MRP
        atti = -att
    elseif parameterization == quaternion
        atti = qinv(att)
    elseif parameterization == DCM
        atti = att'
    end

    for i = 1:size(un,2)
        unb[:,i] = rotFunc(atti,view(un,:,i))
        uub[:,i] = rotFunc(atti,view(uu,:,i))
        uvb[:,i] = rotFunc(atti,view(uv,:,i))
    end

    return unb,uub,uvb
end

function _toInertialFrame(att :: anyAttitude{T}, un :: ArrayOfVecs, uu :: ArrayOfVecs, uv :: ArrayOfVecs, rotFunc :: Function, parameterization) where {T <: Real}

    unb = Array{Array{T,1},1}(undef,length(un))
    uub = Array{Array{T,1},1}(undef,length(uu))
    uvb = Array{Array{T,1},1}(undef,length(uv))

    if parameterization == MRP
        atti = -att
    elseif parameterization == quaternion
        atti = qinv(att)
    elseif parameterization == DCM
        atti = att'
    elseif parameterization == att2D
        atti = -att
    else
        error("please provide valid attitude parameterization")
    end

    for i = 1:length(un)
        unb[i] = rotFunc(atti,un[i])
        uub[i] = rotFunc(atti,uu[i])
        uvb[i] = rotFunc(atti,uv[i])
    end

    return unb,uub,uvb
end

function sMRP(p :: Vec)
    return -p./(dot(p,p))
end

function vecAlignAttGen(v1 :: Vec{T}, v2, tht) where {T <: Real}

    if norm(v1 .+ v2) != 0
        e = (v1 .+ v2)./norm(v1 .+ v2)
    else
        e = cross(v1,v2)./norm(cross(v1,v2))
    end

    atts = Array{T,1}(undef,length(tht))
    for i = 1:length(tht)
        atts[i] = [e*cos(tht[i]/2) + cross(v2,e).*sin(tht[i]/2); -dot(e,v2)*sin(tht[i]/2)]
    end
    return atts
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
    for i = 1:length(atts)
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

"""
    Function to imitate MATLABs identity matrix functionalliy
"""
function eye(dim, T = Float64)
    return Matrix{T}(I,dim,dim)
end

"""
    Function to get the parameterizaiton of an attitude
    Assumes that any array with multiple attitudes will have attitudes stored as columns
"""
function getAttParam(Att, fullState = false)
    if typeof(Att) <: AbstractVector
        l = length(Att)
        if (l == 3) & (!fullState)
            return GRP
        elseif (l == 4) & (!fullState)
            return quaternion
        elseif (l == 6) & (fullState)
            return GRP
        elseif (l == 7) & (fullState)
            return quaternion
        else
            error("Invalide Attitude Parameterization")
        end
    elseif typeof(Att) <: AbstractArray
        s = size(Att)
        if s == (3,3)
            return DCM
        elseif (s[1] == 3) & (!fullState)
            return GRP
        elseif (s[1] == 4) & (!fullState)
            return quaternion
        elseif (s[1] == 6) & (fullState)
            return GRP
        elseif (s[1] == 7) & (fullState)
            return quaternion
        else
            error("Invalide Attitude Parameterization")
        end
    elseif typeof(Att) <: AbstractVector{V} where {V <: AbstractVector}
        l = length(Att[])
        if (l == 3) & (!fullState)
            return GRP
        elseif (l == 4) & (!fullState)
            return quaternion
        elseif (l == 6) & (fullState)
            return GRP
        elseif (l == 7) & (fullState)
            return quaternion
        else
            error("Invalide Attitude Parameterization")
        end
    elseif typeof(Att) <: Union{MRP,GRP,quaternion,DCM}
        return typeof(Att)
    else
        error("Invalide Attitude Parameterization")
    end
end
