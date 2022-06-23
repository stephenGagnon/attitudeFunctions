const Vec{T<:Number} = AbstractArray{T,1} # need these to handle slices # need T <: Real to handle Forward Diff
const Mat{T<:Number} = AbstractArray{T,2}
const Vecs{T<:Number} = Array{V,1} where V <: Vec
const ArrayOfVecs{T<:Number} = Array{V,1} where V <: Vec
const MatOrVecs = Union{Mat,ArrayOfVecs}

"""
    Custom type for quaternions with 2 fields:
    v - the vector part
    s - the scalar part
"""

struct quaternion{T} #where {T <: Real}
    v :: Vec{T}
    s :: T
end

struct test
 a :: Nothing
end

function quaternion(a :: Mat{T}, b :: T) where {T <: Real}
    if (size(a,1) == 3 & size(a,2) == 1) | (size(a,2) == 3 & size(a,1) == 1)
        return quaternion(a[:],b)
    else
        throw(error("Vector of quaternion part must have length 3"))
    end
end

function quaternion(q :: Vec{T}) where {T <: Real}
    return quaternion(q[1:3],q[4])
end

"""
    Custom type for generalized Rodrigues parameters with 3 fields:
    p - the 3 element vector specifying the attitude
    a,f - the parameters specifying the exact GRP transformations
"""
struct GRP{T} #where {T <: Real}
    # GRP values
    p :: Vec{T}
    # a=f=1 gives the standard modified rodrigues parameters
    a :: Number
    f :: Number
end

"""
    Custom type for modified Rodrigues parameters with one field:
    p - the 3 element vector specifying the attitude
"""
struct MRP{T} #where {T <: Real}
    #Modified Rodrigues Parameters
    # MRP values
    p :: Vec{T}
end

"""
    Custom type for direction cosine matrices with one field:
    A - the DCM represented as a 2D array
"""
struct DCM{T} #where {T <: Real}
    A :: Mat{T} #full attitude matrix
end

"""
    custom type for 2D attitude
    a scalar representing rotation around the out-of-plane axis
    type is mostly for handling the 2D case in generalized functions
"""
struct att2D{T} #where {T <: Real}
    tht :: T
end

"""
    Function to retrieve the data type of attitude value
    Supports all types in the Union 'anyAttitude'
"""
function getDataType(x :: quaternion{T}) where{T}
    return T
end

function getDataType(x :: MRP{T}) where{T}
    return T
end

function getDataType(x :: GRP{T}) where{T}
    return T
end

function getDataType(x :: DCM{T}) where{T}
    return T
end

function getDataType(x :: Vec{T}) where{T}
    return T
end

function getDataType(x :: Mat{T}) where{T}
    return T
end

function getDataType(x :: T) where{T}
    return T
end

const anyAttitude{T} = Union{Mat{T},Vec{T},DCM{T},MRP{T},GRP{T},quaternion{T}, att2D{T}} where {T <: Real}
const arrayofAtts = Array{A,1} where A <: anyAttitude
