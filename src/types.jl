const Vec{T<:Number} = AbstractArray{T,1}
const Mat{T<:Number} = AbstractArray{T,2}
const Vecs{T<:Number} = Array{V,1} where V <: Vec
const ArrayOfVecs{T<:Number} = Array{V,1} where V <: Vec
const MatOrVecs = Union{Mat,ArrayOfVecs}
const anyAttitude = Union{Mat,Vec,DCM,MRP,GRP,quaternion}

"""
    Custom type for quaternions with 2 fields:
    v - the vector part
    s - the scalar part
"""
struct quaternion
    v :: Vec#vector part
    s :: N where {N <: Number} #scalar part
end

struct test
 a :: Nothing
end

function quaternion(a :: Mat,b :: N where {N <: Number})
    if (size(a,1) == 3 & size(a,2) == 1) | (size(a,2) == 3 & size(a,1) == 1)
        return quaternion(a[:],b)
    else
        throw(error("Vector of quaternion part must have length 3"))
    end
end

function quaternion(q :: Vec)
    return quaternion(q[1:3],q[4])
end

"""
    Custom type for generalized Rodrigues parameters with 3 fields:
    p - the 3 element vector specifying the attitude
    a,f - the parameters specifying the exact GRP transformations
"""
struct GRP
    # GRP values
    p :: Vec
    # a=f=1 gives the standard modified rodrigues parameters
    a :: N where {N <: Number}
    f :: N where {N <: Number}
end

"""
    Custom type for modified Rodrigues parameters with one field:
    p - the 3 element vector specifying the attitude
"""
struct MRP
    #Modified Rodrigues Parameters
    # MRP values
    p :: Vec
end

"""
    Custom type for direction cosine matrices with one field:
    A - the DCM represented as a 2D array
"""
struct DCM
    A :: Mat #full attitude matrix
end
