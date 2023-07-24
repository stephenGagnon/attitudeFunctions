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

function frameConversion(attitude :: anyAttitude{T}, vectors :: Array{Vec,1}, rotFunc :: Function) where {T <: Real}
    transformedVectors = similar(vectors)
    for i = eachindex(vectors)
        transformedVectors[i] = rotFunc(attitude,vectors[i])
    end
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

function _toBodyFrame(att :: anyAttitude{T}, usun :: Vec, uobs :: Vector{Vector{Float64}}, rotFunc :: Function) where {T <: Real}
    #
    usunb = rotFunc(att,usun)

    uobsb = Array{Array{T,1},1}(undef,length(uobs))

    for i = eachindex(uobs)
        uobsb[i] = rotFunc(att,uobs[i])
    end
    # uobsb = map(x -> rotFunc(att,x), uobs)

    return usunb :: Vec, uobsb :: ArrayOfVecs
end

function _toBodyFrame(att :: anyAttitude, usun :: Vec, uobs :: ArrayOfVecs{T}, rotFunc :: Function) where {T <: Vec}
    #
    usunb = rotFunc(att,usun)

    uobsb = Array{Array{eltype(T),1},1}(undef,length(uobs))

    for i = eachindex(uobs)
        uobsb[i] = rotFunc(att,uobs[i])
    end
    # uobsb = map(x -> rotFunc(att,x), uobs)

    return usunb :: Vec, uobsb :: ArrayOfVecs
end

function _toBodyFrame!(att :: anyAttitude, usun :: Vec, uobs :: ArrayOfVecs, rotFunc :: Function, usunb, uobsb)
    rotFunc(att, usun, usunb)

    for i = eachindex(uobs)
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
    end
    for uvi in uv
           uvb[i] = rotFunc(atti,uvi)
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

function rotate_2D(tht,v)
    return [cos(tht) sin(tht);-sin(tht) cos(tht)]*v
end
