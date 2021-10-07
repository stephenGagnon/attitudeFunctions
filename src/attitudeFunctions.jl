module attitudeFunctions

using LinearAlgebra
using Random
# using Infiltrator

include("types.jl")
include("attConversions.jl")
include("utilities.jl")
include("attitudeDynamics.jl")
include("quaternionFunctions.jl")

export q2A, p2q, q2p, A2q, p2A, A2p, qprod, qinv, attitudeErrors, randomAtt,
    quaternion, GRP, MRP, DCM, any2A, attitude2Array, crossMat, crossMatInv,
    qdq2w, qPropDisc, qRotate, dDotdp, dAdp, toBodyFrame, _toBodyFrame, dqdp,
    dAdq, quaternionDistance, sMRP, randomAttPar, dDotdq, _toInertialFrame



# dA[1,1] = [2*q[1]*(-q[1]^2 + q[2]^2 + q[3]^2 + 1 - q[4]^2);
#            2*q[2]*(-q[1]^2 + q[2]^2 + q[3]^2 - (1 + q[4])^2);
#            2*q[3]*(-q[1]^2 + q[2]^2 + q[3]^2 - (1 + q[4])^2)]
#
# dA[1,2] = [2*q[2]*(-2*q[1]^2 + q[4] + 1) + q[1]*q[3]*(2*q[4] + 1);
#            2*q[1]*(-2*q[2]^2 + q[4] + 1) + q[2]*q[3]*(2*q[4] + 1);
#           -2*q[4]*(-2*q[3]^2 + q[4] + 1) + q[3]*(q[3] - 2*q[1]*q[2])]
#
# dA[1,3] = [2*q[3]*(-2*q[1]^2 + q[4] + 1) - q[1]*q[2]*(2*q[4] + 1);
#             2*q[4]*(-2*q[2]^2 + q[4] + 1) - q[2]*(q[2] + 2*q[1]*q[3]);
#             2*q[1]*(-2*q[3]^2 + q[4] + 1) - q[2]*q[3]*(2*q[4] + 1)]
#
# dA[2,1] = [2*q[2]*(-2*q[1]^2 + q[4] + 1) - q[1]*q[3]*(2*q[4] + 1);
#             2*q[1]*(-2*q[2]^2 + q[4] + 1) - q[2]*q[3]*(2*q[4] + 1);
#             2*q[4]*(-2*q[3]^2 + q[4] + 1) - q[3]*(q[3] + 2*q[1]*q[2])]
#
# dA[2,2] = [2*q[1]*(-q[2]^2 + q[1]^2 + q[3]^2 - (1 + q[4])^2);
#             2*q[2]*(-q[2]^2 + q[1]^2 + q[3]^2 + 1 - q[4]^2);
#             2*q[3]*(-q[2]^2 + q[1]^2 + q[3]^2 - (1 + q[4])^2)]
#
# dA[2,3] = [-2*q[4]*(-2*q[1]^2 + q[4] + 1) + q[1]*(q[1] - 2*q[3]*q[2]);
#             2*q[3]*(-2*q[2]^2 + q[4] + 1) + q[1]*q[2]*(2*q[4] + 1);
#             2*q[2]*(-2*q[3]^2 + q[4] + 1) + q[1]*q[3]*(2*q[4] + 1)]
#
# dA[3,1] = [2*q[3]*(-2*q[1]^2 + q[4] + 1) + q[1]*q[2]*(2*q[4] + 1);
#             -2*q[4]*(-2*q[2]^2 + q[4] + 1) + q[2]*(q[2] - 2*q[1]*q[3]);
#             2*q[1]*(-2*q[3]^2 + q[4] + 1) + q[2]*q[3]*(2*q[4] + 1)]
#
# dA[3,2] = [2*q[4]*(-2*q[1]^2 + q[4] + 1) - q[1]*(q[1] + 2*q[3]*q[2]);
#             2*q[3]*(-2*q[2]^2 + q[4] + 1) - q[1]*q[2]*(2*q[4] + 1);
#             2*q[2]*(-2*q[3]^2 + q[4] + 1) - q[1]*q[3]*(2*q[4] + 1)]
#
# dA[3,3] = [2*q[1]*(-q[3]^2 + q[1]^2 + q[2]^2 - (1 + q[4])^2);
#             2*q[2]*(-q[3]^2 + q[1]^2 + q[2]^2 - (1 + q[4])^2);
#             2*q[3]*(-q[3]^2 + q[1]^2 + q[2]^2 + 1 - q[4]^2)]

end
