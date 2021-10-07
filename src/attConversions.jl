# Contains a set of functions to convert between various attitude parameterizations
# currently, MRP, GRP, quaternion, and DCM are supported

"""
    Converts quaterions to dirrection cosine matricies

    accepted inputs:
    4x1 or 4xn float array where collumns correspond to quaternions
    single or 1xn array of quaternion type

    outputs:
    if input is a float array, output is a 3x3 or 3x3xn float array of direction cosine matrices
    if input is a quaternion array, output is a DCM array of the same size
        only supports 1d arrays
"""
function q2A(q :: Vec)

    A = Array{Float64,2}(undef,3,3)
    A[1,1] = (q[1]^2 - q[2]^2 - q[3]^2 + q[4]^2)
    A[1,2] = (2*(q[1]*q[2] + q[3]*q[4]))
    A[1,3] = (2*(q[1]*q[3] - q[2]*q[4]))
    A[2,1] = (2*(q[1]*q[2] - q[3]*q[4]))
    A[2,2] = (-q[1]^2 + q[2]^2 - q[3]^2 + q[4]^2)
    A[2,3] = (2*(q[2]*q[3] + q[1]*q[4]))
    A[3,1] = (2*(q[1]*q[3] + q[2]*q[4]))
    A[3,2] = (2*(q[2]*q[3] - q[1]*q[4]))
    A[3,3] = (-q[1]^2 - q[2]^2 + q[3]^2 + q[4]^2)

    return A
end

function q2A(q :: Mat)

    A = Array{Float64,3}(undef,3,3,size(q,2))

    for i = 1:size(q,2)
        A[:,:,i] = q2A(q[:,i])
    end
    return A
end

function q2A(q :: Vecs)

    A = Array{Array{Float64,2},1}(undef,length(q))

    for i = 1:length(q)
        A[i] = q2A(q[i])
    end
    return A
end

function q2A(q :: quaternion)
    A = Array{Float64,2}(undef,3,3)
    A[1,1] = (q.v[1]^2 - q.v[2]^2 - q.v[3]^2 + q.s^2)
    A[1,2] = (2*(q.v[1]*q.v[2] + q.v[3]*q.s))
    A[1,3] = (2*(q.v[1]*q.v[3] - q.v[2]*q.s))
    A[2,1] = (2*(q.v[1]*q.v[2] - q.v[3]*q.s))
    A[2,2] = (-q.v[1]^2 + q.v[2]^2 - q.v[3]^2 + q.s^2)
    A[2,3] = (2*(q.v[2]*q.v[3] + q.v[1]*q.s))
    A[3,1] = (2*(q.v[1]*q.v[3] + q.v[2]*q.s))
    A[3,2] = (2*(q.v[2]*q.v[3] - q.v[1]*q.s))
    A[3,3] = (-q.v[1]^2 - q.v[2]^2 + q.v[3]^2 + q.s^2)
    return DCM(A)
end

function q2A(q :: Array{quaternion,1})
    A = Array{DCM,1}(undef,length(q))

    for i = 1:size(q,2)
        A[i] = q2A(q[i])
    end
    return A
end

"""
    Function to convert rodrigues parameters to quaternions

    accepted inputs:
    float array of rodrigues parameters where each column is a 3 element MRP or GRP
        (a=f=1 is assumed unless provided)
    MRP type or 1D array of MRP types
    GRP type or 1D array of GRP types

    outputs:
    if a float array is provided, return a float array of rodrigues parameters
        where the ith column corresponds to the ith column of the input array
    if an MRP or GRP type array is provided, returns a quaternion array of the same size
        only supports 1d arrays
"""
function p2q(p :: Vec, a=1, f=1)

    q = Array{Float64,1}(undef,4)
    pd = dot(p,p)
    q[4] = (-a*pd + f*sqrt(f^2 + (1-a^2)*pd))/(f^2 + pd)
    # q[1:3] = (a + q[4]).*p./f
    q[1] = (a + q[4])*p[1]/f
    q[2] = (a + q[4])*p[2]/f
    q[3] = (a + q[4])*p[3]/f
    return q
end

function p2q(p :: Mat, a=1, f=1)

    q = zeros(4,size(p,2))
    for i = 1:size(p,2)
        q[:,i] = p2q(p[:,i],a,f)
    end
    return q
end

function p2q(p :: MRP)

    qv = Array{Float64,1}(undef,3)
    pd = p.p'*p.p
    qs = (-pd + 1)/(1 + pd)
    qv = (1 + qs).*(p.p)
    return quaternion(qv,qs)
end

function p2q(p :: Array{MRP,1})

    q = Array{quaternion,1}(undef,length(p))
    for i = 1:size(p)
        q[i] = p2q(p[i])
    end
    return q
end

function p2q(p :: GRP)

    qv = Array{Float64,1}(undef,3)
    pd = p.p'*p.p
    qs = (-p.a*pd + p.f*sqrt(p.f^2 + (1-p.a^2)*pd))/(p.f^2 + pd)
    qv = (p.a + qs).*(p.p)./p.f
    return quaternion(qv,qs)
end

function p2q(p :: Array{GRP,1})

    q = Array{quaternion,1}(undef,length(p))
    for i = 1:size(p)
        q[i] = p2q(p[i])
    end
    return q
end

"""
    Function to convert quaternions to rodrigues parameters

    accpeted inputs:
    4x1 or 4xn float array
    quaternion type or 1xn array of quaternion type

    outputs:
    if float arrays are provided, return a float array of rodrigues parameters
        where the ith column corresponds to the ith column of the quaternion array
        a=f=1 is assumed unless provided
    if a quaternion type is provided with a | f != 1, returns a GRP type
    if a quaternion type is provided without a and f values, returns an MRP type
    if a quaternion type with a=f=1 is provided, returns an MRP type
"""
function q2p(q :: Vec, a = 1, f = 1)
    return f*q[1:3]./(a + q[4])
end

function q2p(q :: Mat, a = 1, f = 1)

    p = zeros(3,size(q,2))
    for i = 1:size(q,2)
        p[:,i] = q2p(q[:,i],a,f)
    end
    return p
end

function q2p(q :: Vecs, a = 1, f = 1)

    p = Array{Array{typeof(q[1][1]),1},1}(undef,length(q))
    for i = 1:length(q)
        p[i] = q2p(q[i],a,f)
    end
    return p
end

function q2p(q :: quaternion)
    return MRP((q.v)./(1 + q.s))
end

function q2p(q :: Array{quaternion,1})

    p = Array{MRP,1}(undef,length(q))
    for i = 1:length(q)
        p[i] = q2p(q[i])
    end
    return p
end

function q2p(q :: quaternion, a, f)

    if a | f != 1
        return GRP(f.*(q.v)./(a + q.s),a,f)
    else
        return q2p(q)
    end
end

function q2p(q :: Array{quaternion,1}, a, f)

    if a | f !=1
        p = Array{GRP,1}(undef,length(q))
        for i = 1:size(q,2)
            p[i] = q2p(q[i],a,f)
        end
        return p
    else
        return q2p(q)
    end
end


"""
    convert direction cosine matrices to quaternions

    accpeted inputs:
    3x3 or 3x3xn float array
    DCM type or 1xn array of DCM type

    outputs:
    if float arrays are provided, return a float array of quaternions
        where the ith column corresponds to the 3x3xith element of the DCM array
    if a DCM type array is provided returns a quaternion type array of the same size
"""
function A2q(A :: Mat)

    q = Array{Float64,1}(undef,4)
    q[4] = .5*sqrt(1 + tr(A))
    q[1] = .25*(A[2,3]-A[3,2])/q[4]
    q[2] = .25*(A[3,1]-A[1,3])/q[4]
    q[3] = .25*(A[1,2]-A[2,1])/q[4]
    return q
end

function A2q(A :: T where {Num <: Number, T <: AbstractArray{Num,3}})

    q = Array{Float64,2}(undef,4,size(A,3))

    for i = 1:size(A,3)
        q[:,i] = A2q(A[:,:,i])
    end
    return q
end

function A2q(A :: DCM)

    qs = .5*sqrt(1 + tr(A.A))
    qv = Array{Float64,1}(undef,3)
    qv[1] = .25*(A.A[2,3]-A.A[3,2])/qs
    qv[2] = .25*(A.A[3,1]-A.A[1,3])/qs
    qv[3] = .25*(A.A[1,2]-A.A[2,1])/qs
    return quaternion(qv,qs)
end

function A2q(A :: Array{DCM,1})

    q = Array{quaternion,1}(undef,length(A))

    for i = 1:size(A)
        q[i] = A2q(A[i])
    end
    return q
end


"""
    convert rodrigues parameters to direction cosine matricies

    accepted inputs:
    float array of rodrigues parameters where each column is a 3 element MRP or GRP
        (a=f=1 is assumed unless provided)
    MRP type or 1D array of MRP types
    GRP type or 1D array of GRP types

    outputs:
    if a float array is provided, return a float array of rodrigues parameters
        where the ith column corresponds to the 3x3xith element of the DCM array
    if an MRP or GRP type array is provided, returns a DCM array of the same size
"""
function p2A(p :: Vec, a = 1.0, f = 1.0)
    pd = dot(p,p)
    q4 = (-a*pd + f*sqrt(f^2 + (1-a^2)*pd))/(f^2 + pd)
    A = Array{Float64,2}(undef,3,3)
    A[1,1] = (((a + q4)*p[1]/f)^2 - ((a + q4)*p[2]/f)^2 - ((a + q4)*p[3]/f)^2 + q4^2)
    A[1,2] = (2*(((a + q4)*p[1]/f)*((a + q4)*p[2]/f) + ((a + q4)*p[3]/f)*q4))
    A[1,3] = (2*(((a + q4)*p[1]/f)*((a + q4)*p[3]/f) - ((a + q4)*p[2]/f)*q4))
    A[2,1] = (2*(((a + q4)*p[1]/f)*((a + q4)*p[2]/f) - ((a + q4)*p[3]/f)*q4))
    A[2,2] = (-((a + q4)*p[1]/f)^2 + ((a + q4)*p[2]/f)^2 - ((a + q4)*p[3]/f)^2 + q4^2)
    A[2,3] = (2*(((a + q4)*p[2]/f)*((a + q4)*p[3]/f) + ((a + q4)*p[1]/f)*q4))
    A[3,1] = (2*(((a + q4)*p[1]/f)*((a + q4)*p[3]/f) + ((a + q4)*p[2]/f)*q4))
    A[3,2] = (2*(((a + q4)*p[2]/f)*((a + q4)*p[3]/f) - ((a + q4)*p[1]/f)*q4))
    A[3,3] = (-((a + q4)*p[1]/f)^2 - ((a + q4)*p[2]/f)^2 + ((a + q4)*p[3]/f)^2 + q4^2)
    return A
    # return q2A(p2q(p,a,f))
end

function p2A(p :: Mat, a = 1, f = 1)

    A = zeros(3,3,size(p,2))
    for i = 1:size(p,2)
        A[:,:,i] = p2A(view(p,:,i),a,f)
    end
    return A
end

function p2A(p :: Union{MRP,GRP}, a=1, f=1)
    return q2A(p2q(p))
end

function p2A(p :: Array{Union{MRP,GRP},1}, a=1, f=1)

    A = Array{DCM,1}(undef,length(p))
    for i = 1:size(p)
        A[i] = q2A(p2q(p[i]))
    end
    return A
end

"""
    convert direction cosine matrix to Rodrigues parameters

    accpeted inputs:
    3x3 or 3x3xn float array
    DCM type or 1xn array of DCM type

    outputs:
    if float arrays are provided, return a float array of rodrigues parameters
        where the ith column corresponds to the 3x3xith element of the DCM array
        a=f=1 is assumed unless provided
    if a DCM type is provided with a | f != 1, returns a GRP type
    if a DCM is provided without a and f values, returns an MRP type
    if a DCM with a=f=1 is provided, returns an MRP type
"""
function A2p(A :: Mat, a = 1, f = 1)
    return q2p(A2q(A),a,f)
end

function A2p(A :: T where {Num <: Number, T <: AbstractArray{Num,3}},
     a = 1, f = 1)

    q = Array{Float64,2}(undef,4,size(A,3))

    for i = 1:size(A,3)
        q[:,i] = A2q(A[:,:,i])
    end
    return q
end

function A2p(A :: DCM)

    return q2p(A2q(A))
end

function A2p(A :: Array{DCM,1})

    p = Array{MRP,1}(undef,length(A))

    for i = 1:size(A,3)
        p[i] = A2p(A[i])
    end
    return p
end

function A2p(A :: DCM, a, f)
    if a | f != 1
        return q2p(A2q(A),a,f)
    else
        return q2p(A2q(A))
    end
end

function A2p(A :: Array{DCM,1}, a, f)

    p = Array{GRP,1}(undef,length(A))

    for i = 1:size(A,3)
        p[i] = A2p(A[i])
    end
    return p
end

"""
    functions to convert any attitude type to a DCM
"""
function any2A(att :: quaternion)
    return q2A(att)
end

function any2A(att :: Array{quaternion,1})
    A = Array{DCM,1}(undef,length(att))
    for i = 1:length(att)
        A[i] = q2A(att[i])
    end
    return A
end

function any2A(att :: Union{MRP,GRP})
    return p2A(att)
end

function any2A(att :: Union{Array{GRP,1},Array{MRP,1}})

    A = Array{DCM,1}(undef,length(att))
    for i = 1:length(att)
        A[i] = p2A(att[i])
    end
    return A
end

function any2A(att :: DCM)
    return att
end

function any2A(att :: Mat, attType = DCM, a = 1.0 , f = 1.0)
    if attType == DCM
        return att
    elseif (attType == MRP) | (attType == GRP)
        return p2A(att,a,f)
    elseif attType == quaternion
        return q2A(att)
    end
end

function any2A(att :: Vec, attType = MRP, a = 1.0 , f = 1.0)
    if (attType == MRP) | (attType == GRP)
        return p2A(att,a,f)
    elseif attType == quaternion
        return q2A(att)
    end
end
