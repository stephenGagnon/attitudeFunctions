module attitudeFunctions

using LinearAlgebra
using Random
#using Infiltrator

export q2A, p2q, q2p, A2q, p2A, A2p, qprod, qinv, attitudeErrors, randomAtt,
    quaternion, GRP, MRP, DCM


"""
    Custom type for quaternions with 2 fields:
    v - the vector part
    s - the scalar part
"""
struct quaternion
    v :: Array{Float64,1}#vector part
    s :: Float64 #scalar part
end

function quaternion(a :: Array{Float64,2},b :: Float64)
    if (size(a,1) == 3 & size(a,2) == 1) | (size(a,2) == 3 & size(a,1) == 1)
        return quaternion(a[:],b)
    else
        throw(error("Vector of quaternion part must have length 3"))
    end
end

"""
    Custom type for generalized Rodrigues parameters with 3 fields:
    p - the 3 element vector specifying the attitude
    a,f - the parameters specifying the exact GRP transformations
"""
struct GRP
    # GRP values
    p :: Array{Float64,1}
    # a=f=1 gives the standard modified rodrigues parameters
    a :: Float64
    f :: Float64
end

"""
    Custom type for modified Rodrigues parameters with one field:
    p - the 3 element vector specifying the attitude
"""
struct MRP
    #Modified Rodrigues Parameters
    # MRP values
    p :: Array{Float64,1}
end

"""
    Custom type for direction cosine matrices with one field:
    A - the DCM represented as a 2D array
"""
struct DCM
    A :: Array{Float64,2} #full attitude matrix
end

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
function q2A(q :: Array{Float64,1})

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

function q2A(q :: Array{Float64,2})

    A = Array{Float64,3}(undef,3,3,size(q,2))

    for i = 1:size(q,2)
        A[:,:,i] = q2A(q[:,i])
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
    A = Array{DCM,1}(undef,size(q))

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
function p2q(p :: Array{Float64,1}, a, f)

    q = Array{Float64,1}(undef,4)
    pd = dot(p,p)
    q[4] = (-a*pd + f*sqrt(f^2 + (1-a^2)*pd))/(f^2 + pd)
    q[1:3] = (a + q[4]).*p./f
    return q
end

function p2q(p :: Array{Float64,2}, a, f)

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

    q = Array{quaternion,1}(undef,size(p))
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

    q = Array{quaternion,1}(undef,size(p))
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
function q2p(q :: Array{Float64,1}, a = 1, f = 1)
    return f*q[1:3]./(a + q[4])
end

function q2p(q :: Array{Float64,2}, a = 1, f = 1)

    p = zeros(3,size(q,2))
    for i = 1:size(q,2)
        p[:,i] = q2p(q[:,i],a,f)
    end
    return p
end

function q2p(q :: quaternion)
    return MRP((q.v)./(1 + q.s))
end

function q2p(q :: Array{quaternion,1})

    p = Array{MRP,1}(undef,size(q))
    for i = 1:size(q,2)
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
        p = Array{GRP,1}(undef,size(q))
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
function A2q(A :: Array{Float64,2})

    q = Array{Float64,1}(undef,4)
    q[4] = .5*sqrt(1 + tr(A))
    q[1] = .25*(A[2,3]-A[3,2])/q[4]
    q[2] = .25*(A[3,1]-A[1,3])/q[4]
    q[3] = .25*(A[1,2]-A[2,1])/q[4]
    return q
end

function A2q(A :: Array{Float64,3})

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

    q = Array{quaternion,1}(undef,size(A))

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
function p2A(p :: Array{Float64,1}, a = 1, f = 1)

    return q2A(p2q(p,a,f))
end

function p2A(p :: Array{Float64,2}, a = 1, f = 1)

    A = zeros(3,3,size(p,2))
    for i = 1:size(p,2)
        A[:,:,i] = q2A(p2q(p[:,i],a,f))
    end
    return A
end

function p2A(p :: Union{MRP,GRP})
    return q2A(p2q(p))
end

function p2A(p :: Array{Union{MRP,GRP},1})

    A = Array{DCM,1}(undef,size(p))
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
function A2p(A :: Array{Float64,2}, a = 1, f = 1)
    return q2p(A2q(A),a,f)
end

function A2p(A :: Array{Float64,3}, a = 1, f = 1)

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

    p = Array{MRP,1}(undef,size(A))

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

    p = Array{GRP,1}(undef,size(A))

    for i = 1:size(A,3)
        p[i] = A2p(A[i])
    end
    return p
end

"""
    Computes the quaternion product of two quaternions

    accepted inputs:
    two 4x1 or 4xn float arrays where each column is a quaternion
        inputs must be either the same size, or a single quaternion and an array
        the single quaternion must be in a 1D (not 2D) float array

    two quaternion arrays, these must be the same size or a single quaternion and
        an array

    outputs:
    if inputs are float arrays, output is a float array where columns correspond
        to the quaternion product of the inputs
    if the inputs are quaternion type arrays, then the outputs are of quaternion type
"""
function qprod(q1 :: Array{Float64,1}, q2 :: Array{Float64,1})
    qp = zeros(4,1)
    qp[1] =  q1[4]*q2[1] + q1[3]*q2[2] - q1[2]*q2[3] + q1[1]*q2[4]
    qp[2] = -q1[3]*q2[1] + q1[4]*q2[2] + q1[1]*q2[3] + q1[2]*q2[4]
    qp[3] =  q1[2]*q2[1] - q1[1]*q2[2] + q1[4]*q2[3] + q1[3]*q2[4]
    qp[4] = -q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3] + q1[4]*q2[4]
    return qp
end

function qprod(q1 :: Array{Float64,2}, q2 :: Array{Float64,2})

    if size(q1,1) < 4 | size(q1,1) > 4
        q1 = q1'
    end
    if size(q2,1) < 4 | size(q2,1) > 4
        q1 = q1'
    end

    if size(q1,2) == 1
        qp = zeros(4,size(q2,2))
        for i = 1:size(q2,2)
            qp[:,i] = qprod(q1[:],q2[:,i])
        end
        return qp
    elseif size(q2,2) == 1
        qp = zeros(4,size(q1,2))
        for i = 1:size(q1,2)
            qp[:,i] = qprod(q1[:,i],q2[:])
        end
        return qp
    elseif size(q1,2) == size(q2,2)
        qp = zeros(4,size(q2,2))
        for i = 1:size(q2,2)
            qp[:,i] = qprod(q1[:,i],q2[:,i])
        end
        return qp
    else
        throw(ArgumentError("quaternion input matricies must be the same size or one dimensional"))
    end
end

function qprod(q1 :: quaternion, q2 :: quaternion)
    # return quaternion((cross(q1.v,q2.v) + (q1.s .* q2.v) + (q2.s .* q1.v)),(q1.s*q2.s - dot(q1.v,q2.v)))
    qv = Vector{Float64}(undef,3)
    qv[1] =  q1.s*q2.v[1] + q1.v[3]*q2.v[2] - q1.v[2]*q2.v[3] + q1.v[1]*q2.s
    qv[2] = -q1.v[3]*q2.v[1] + q1.s*q2.v[2] + q1.v[1]*q2.v[3] + q1.v[2]*q2.s
    qv[3] =  q1.v[2]*q2.v[1] - q1.v[1]*q2.v[2] + q1.s*q2.v[3] + q1.v[3]*q2.s
    qs = -q1.v[1]*q2.v[1] - q1.v[2]*q2.v[2] - q1.v[3]*q2.v[3] + q1.s*q2.s
    return quaternion(qv,qs)
end

function qprod(q1 :: Union{Array{Array{quaternion,1}, 1},quaternion},
               q2 :: Union{Array{Array{quaternion,1}, 1},quaternion})
    if typeof(q1) == quaternion
        qp = Array{quaternion,1}(undef,size(q2))
        for i = 1:size(q2)
            qp[i] = qprod(q1,q2[i])
        end
        return qp
    elseif typeof(q2) == quaternion
        qp = Array{quaternion,1}(undef,size(q1))
        for i = 1:size(q1,2)
            qp[i] = qprod(q1[i],q2)
        end
        return qp
    elseif size(q1) == size(q2)
        qp = Array{quaternion,1}(undef,size(q2))
        for i = 1:size(q2)
            qp[i] = qprod(q1[i],q2[i])
        end
        return qp
    else
        throw(ArgumentError("quaternion input vectors must be the same size or singular"))
    end
end

"""
    computes the inverse of a quaternion

    accepted inputs:
    a 4x1 or 4xn float array where each collumn is a quaternion
    a singleton or 1xn array of quaternion type

    outputs:
    if input is a float array, output is a float array where collumns correspond
        to the quaternion inverse of the input array collumn
    if the input is a quaternion type array, then the output is alsoa  quaternion
        type array
"""
function qinv(q :: Array{Float64,1})

    qi = zeros(4,1)
    qi[1:3] = -q[1:3]
    qi[4] = q[4]
    return qi[:]
end

function qinv(q :: Array{Float64,2})
    if size(q,1) < 4 | size(q,1) > 4
        q = q'
    end
    qi = zeros(size(q))
    for i = 1:size(q,2)
        qi[:,i] = qinv(q[:,i])
    end
    return qi
end

function qinv(q :: quaternion)

    return quaterion(-q.v,q.s)
end

function qinv(q :: Array{quaternion,1})
    qi = Array{quaternion,1}(undef,size(q))
    for i = 1:size(q)
        qi[i] = qinv(q[i])
    end
    return qi
end

"""
    Computes the error between two quaternions by finding the rotation between them
    and returning 2x the vector part of the quaternion. This corresponds to the roll
    pitch yaw errors for small error angles.
"""
function attitudeErrors(qE :: Array{Float64,1},q :: Array{Float64,1})

    if !(sign(qE[4]) == sign(q[4]) == -1) | !(sign(qE[4]) == sign(q[4]) == 1)
            qE = -qE
    end
    return 2*qprod(qE,qinv(q))
end

function attitudeErrors(q1 :: Union{Array{Array{Float64,1}, 1},Array{Float64,1}},
                        q2 :: Union{Array{Array{Float64,1}, 1},Array{Float64,1}})

    if size(q,1) < 4 | size(q,1) > 4
        q = q'
    end

    if size(qE,1) < 4 | size(qE,1) > 4
        qE = qE'
    end

    dalpha = zeros(3,size(q,2))

    if size(qE,2) == 1
        for i = 1:size(q,2)
            dalpha[:,i] = attitudeErrors(qE[:],q[:,i])
        end
    elseif size(qE,2) == 1
        for i = 1:size(qE,2)
            dalpha[:,i] = attitudeErrors(qE[:,i],q[:])
        end
    else
        for i = 1:size(q,2)
            dalpha[:,i] = attitudeErrors(qE[:,i],q[:,i])
        end
    end

    return dalpha
end

function attitudeErrors(qE :: quaternion, q :: quaternion)

    if !(sign(qE.s) == sign(q.s) == -1) | !(sign(qE.s) == sign(q.s) == 1)
            qE = -qE
    end
    return 2*qprod(qE,qinv(q))
end

function attitudeErrors(q1 :: Union{Array{quaternion, 1}, quaternion},
                        q2 :: Union{Array{quaternion, 1}, quaternion})

    if typeof(qE) == quaternion
        dalpha = Array{quaternion,1}(undef,size(q))
        for i = 1:size(q,2)
            dalpha[i] = attitudeErrors(qE,q[i])
        end
    elseif typeof(qE) == quaternion
        dalpha = Array{quaternion,1}(undef,size(qE))
        for i = 1:size(qE,2)
            dalpha[i] = attitudeErrors(qE[i],q)
        end
    else
        dalpha = Array{quaternion,1}(undef,size(q))
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
function randomAtt(N :: Int64, T; customTypes = false, a = 1, f = 1)

    val = lhs(3,N)

    val[2:3,:] = val[2:3,:].*2*pi;

    q = zeros(4,N)

    q[1,:] = sqrt.(val[1,:]).*cos.(val[2,:])
    q[2,:] = sqrt.(val[1,:]).*sin.(val[2,:])
    q[3,:] = sqrt.(1 .- val[1,:]).*sin.(val[3,:])
    q[4,:] = sqrt.(1 .- val[1,:]).*cos.(val[3,:])

    if customTypes
        q = [quaternion(col[1:3],col[4]) for col in eachcol(q)]
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

"""
    latin hypercube sampling:
    generates N samples with d dimensions
    samples have d elements with values in the range 0-1
    there is guarenteed to be one sample with the value in each dimesnion in the
    ranges 0-1/d, 1/d-2/d, ... , (1-d)/d-1
"""
function lhs(N :: Int64, d :: Int64)

    x = rand(N,d)

    x = x.*1/d .+ reshape((collect(0:d-1))./d,1,d)
    x[x.>1] .= 1.0

    return transpose(hcat([row[randperm(d)] for row in eachrow(x)]...))
end

end
