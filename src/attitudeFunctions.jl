module attitudeFunctions

using LinearAlgebra
using Random
# using Infiltrator

include("types.jl")
include("attConversions.jl")
include("utilities.jl")
include("attitudeDynamics.jl")
include("quaternionFunctions.jl")

export q2A, p2q, q2p, A2q, p2A, A2p, qprod, qinv, attitudeErrors, randomAtt, quaternion, GRP, MRP, DCM, att2D, anyAttitude, arrayofAtts, any2A, attitude2Array, attitude2Vecs, crossMat, crossMatInv, qdq2w, qPropDisc, qRotate, dDotdp, dAdp, toBodyFrame, _toBodyFrame, dqdp, dAdq, quaternionDistance, sMRP, randomAttPar, dDotdq, _toInertialFrame, attDyn, vecAlignAttGen, randomQGen, attitudeRoughening, randomBoundedAngularVelocity, eye, getDataType, attitudeConverter, getAttParam

end
