# this file contains functions for the manipulation and combination of quaternions

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
function qprod(q1 :: Vec{T}, q2 :: Vec) where {T <: Real}
    # qp = zeros(4,)
    qp = Array{T,1}(undef,4)
    qp[1] =  q1[4]*q2[1] + q1[3]*q2[2] - q1[2]*q2[3] + q1[1]*q2[4]
    qp[2] = -q1[3]*q2[1] + q1[4]*q2[2] + q1[1]*q2[3] + q1[2]*q2[4]
    qp[3] =  q1[2]*q2[1] - q1[1]*q2[2] + q1[4]*q2[3] + q1[3]*q2[4]
    qp[4] = -q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3] + q1[4]*q2[4]
    return qp
end

function qprod(q1 :: Union{Mat,Vec}, q2 :: Union{Mat,Vec})

    if size(q1,1) < 4 | size(q1,1) > 4
        q1 = q1'
    end
    if size(q2,1) < 4 | size(q2,1) > 4
        q1 = q1'
    end

    if (length(size(q1)) == 1) | size(q1,2) == 1
        qp = zeros(4,size(q2,2))
        for i = 1:size(q2,2)
            qp[:,i] = qprod(q1[:],q2[:,i])
        end
        return qp
    elseif (length(size(q2)) == 1) | size(q2,2) == 1
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
    qv = Vector{getDataType(q1)}(undef,3)
    qv[1] =  q1.s*q2.v[1] + q1.v[3]*q2.v[2] - q1.v[2]*q2.v[3] + q1.v[1]*q2.s
    qv[2] = -q1.v[3]*q2.v[1] + q1.s*q2.v[2] + q1.v[1]*q2.v[3] + q1.v[2]*q2.s
    qv[3] =  q1.v[2]*q2.v[1] - q1.v[1]*q2.v[2] + q1.s*q2.v[3] + q1.v[3]*q2.s
    qs = -q1.v[1]*q2.v[1] - q1.v[2]*q2.v[2] - q1.v[3]*q2.v[3] + q1.s*q2.s
    return quaternion(qv,qs)
end

function qprod(q1 :: Union{Array{Array{quaternion,1}, 1},quaternion},
               q2 :: Union{Array{Array{quaternion,1}, 1},quaternion})
    if typeof(q1) == quaternion
        qp = Array{quaternion,1}(undef,length(q2))
        for i = 1:size(q2)
            qp[i] = qprod(q1,q2[i])
        end
        return qp
    elseif typeof(q2) == quaternion
        qp = Array{quaternion,1}(undef,length(q1))
        for i = 1:size(q1,2)
            qp[i] = qprod(q1[i],q2)
        end
        return qp
    elseif size(q1) == size(q2)
        qp = Array{quaternion,1}(undef,length(q2))
        for i = 1:size(q2)
            qp[i] = qprod(q1[i],q2[i])
        end
        return qp
    else
        throw(ArgumentError("quaternion input vectors must be the same size or singular"))
    end
end

"""
    Rotate a vector by a quaternion

    inputs:
    q - quaternion
    v - vector

    outputs:
    vp - rotated vector
"""
function qRotate(q :: Vec{T}, v :: Vec) where {T <: Real}
    vp = Array{T,1}(undef,3)
    qRotate!(q,v,vp)
    return vp
end

function qRotate(q :: Vec, v :: Vec, vp :: Vec)
    vp[1] = -(-q[1]*v[1] - q[2]*v[2] - q[3]*v[3])*q[1] -
            (q[2]*v[1] - q[1]*v[2] + q[4]*v[3])*q[2] +
            (-q[3]*v[1] + q[4]*v[2] + q[1]*v[3])*q[3] +
            (q[4]*v[1] + q[3]*v[2] - q[2]*v[3])*q[4]
    vp[2] =  (q[2]*v[1] - q[1]*v[2] + q[4]*v[3])*q[1] -
            (-q[1]*v[1] - q[2]*v[2] - q[3]*v[3])*q[2] -
            (q[4]*v[1] + q[3]*v[2] - q[2]*v[3])*q[3] +
             (-q[3]*v[1] + q[4]*v[2] + q[1]*v[3])*q[4]
    vp[3] = -(-q[3]*v[1] + q[4]*v[2] + q[1]*v[3])*q[1] +
            (q[4]*v[1] + q[3]*v[2] - q[2]*v[3])*q[2] -
            (-q[1]*v[1] - q[2]*v[2] - q[3]*v[3])*q[3] +
            (q[2]*v[1] - q[1]*v[2] + q[4]*v[3])*q[4]
    return vp
end

function qRotate!(q, v, vp)
    vp[1] = -(-q[1]*v[1] - q[2]*v[2] - q[3]*v[3])*q[1] -
            (q[2]*v[1] - q[1]*v[2] + q[4]*v[3])*q[2] +
            (-q[3]*v[1] + q[4]*v[2] + q[1]*v[3])*q[3] +
            (q[4]*v[1] + q[3]*v[2] - q[2]*v[3])*q[4]
    vp[2] =  (q[2]*v[1] - q[1]*v[2] + q[4]*v[3])*q[1] -
            (-q[1]*v[1] - q[2]*v[2] - q[3]*v[3])*q[2] -
            (q[4]*v[1] + q[3]*v[2] - q[2]*v[3])*q[3] +
             (-q[3]*v[1] + q[4]*v[2] + q[1]*v[3])*q[4]
    vp[3] = -(-q[3]*v[1] + q[4]*v[2] + q[1]*v[3])*q[1] +
            (q[4]*v[1] + q[3]*v[2] - q[2]*v[3])*q[2] -
            (-q[1]*v[1] - q[2]*v[2] - q[3]*v[3])*q[3] +
            (q[2]*v[1] - q[1]*v[2] + q[4]*v[3])*q[4]
end

function qRotate(q :: quaternion, v :: Vec)
    q1 = zeros(getDataType(q), 4)
    q1[1] =  q.s*v[1] + q.v[3]*v[2] - q.v[2]*v[3]
    q1[2] = -q.v[3]*v[1] + q.s*v[2] + q.v[1]*v[3]
    q1[3] =  q.v[2]*v[1] - q.v[1]*v[2] + q.s*v[3]
    q1[4] = -q.v[1]*v[1] - q.v[2]*v[2] - q.v[3]*v[3]

    vp = zeros(getDataType(q), 3)
    vp[1] = -q1[4]*q.v[1] - q1[3]*q.v[2] + q1[2]*q.v[3] + q1[1]*q.s
    vp[2] =  q1[3]*q.v[1] - q1[4]*q.v[2] - q1[1]*q.v[3] + q1[2]*q.s
    vp[3] = -q1[2]*q.v[1] + q1[1]*q.v[2] - q1[4]*q.v[3] + q1[3]*q.s

    return vp
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
function qinv(q :: Vec{T}) where {T <: Real}

    qi = Array{T,1}(undef,4)
    qi[1:3] = -q[1:3]
    qi[4] = q[4]
    return qi[:]
end

function qinv(q :: Mat{T}) where {T <: Real}
    if size(q,1) < 4 | size(q,1) > 4
        q = q'
    end
    qi = Array{T,1}(undef,4)
    for i = 1:size(q,2)
        qi[:,i] = qinv(q[:,i])
    end
    return qi
end

function qinv(q :: quaternion)

    return quaterion(-q.v,q.s)
end

function qinv(q :: Array{quaternion,1})
    qi = Array{quaternion,1}(undef,length(q))
    for i = 1:size(q)
        qi[i] = qinv(q[i])
    end
    return qi
end

"""
    funciton for computing a distance metric between two quaternions
    calling with a single array of quaternions (represented as vectors) will return a matrix where the i,jth element is the distance between the ith and jth quaternions in the array
    calling with two arrays of quaternions will compute the distance between each possible pair of quaternions from the two sets
"""
function quaternionDistance(q :: ArrayOfVecs{Vec{T}}) where {T <: Real}
    dist = Array{T,2}(undef,length(q),length(q))
    for i = 1:length(q)
        for j = i:length(q)
            dist[i,j] = 1 - abs(q[i]'*q[j])
            dist[j,i] = dist[i,j]
        end
    end
    return dist
end

function quaternionDistance(q1 :: ArrayOfVecs{Vec{T}}, q2 :: ArrayOfVecs) where {T <: Real}
    dist = Array{T,2}(undef,length(q1),length(q2))
    for i = 1:length(q1)
        for j = i:length(q2)
            dist[i,j] = 1 - abs(q[i]'*q[j])
            dist[j,i] = dist[i,j]
        end
    end
    return dist
end

"""
    function for computing an average quaternion from a set of quaternions
    Markely showed that the average quaternion can be found by computing a matrix as the weighted sum of the outer product of each quaternion with itself, and then taking the eigenvector of that matrix corresponding to the maximum eigenvalue.
"""
function quaternionAverage(q::Array{V}, w) where {T<:Real,V<:Vec{T}}
    M = zeros(4, 4)
    for i = 1:lastindex(q)
        M += w[i] * (q[i] * q[i]')
    end
    # if any(isnan.(M)) | any(M .== Inf)
    #     print("test\n")
    #     @infiltrate
    # end
    E = eigen(M)
    return E.vectors[:, argmax(E.values)]
end

function quaternionAverage(q::Array{V}) where {T<:Real,V<:Vec{T}}
    M = zeros(4, 4)
    for i = 1:lastindex(q)
        M += (q[i] * q[i]')
    end

    E = eigen(M)
    return E.vectors[:, argmax(E.values)]
end