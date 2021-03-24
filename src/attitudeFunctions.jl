module attitudeFunctions

export q2A, p2q, q2p, A2q, p2A, A2p, qprod, qinv

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

    A = Array{Float64,3}(undef,3,3)

    for i = 1:size(q,2)
        A[:,:,i] = q2A(q[:,i])
    end
    return A
end

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

function q2p(q :: Array{Float64,1}, a, f)
    return f*q[1:3]./(a + q[4])
end

function q2p(q :: Array{Float64,2}, a, f)

    p = zeros(3,size(q,2))
    for i = 1:size(q,2)
        p[:,i] = q2p(q[:,i],a,f)
    end
    return p
end

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

function p2A(p :: Array{Float64,1}, a, f)
    return q2A(p2q(p,a,f))
end

function A2p(A :: Array{Float64,2}, a, f)

    return q2p(A2q(A),a,f)
end

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

function qinv(q :: Array{Float64,1})

    qi = zeros(size(q))
    qi[1:3] = -q[1:3]
    qi[4] = q[4]
    return qi
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

function attitudeErrors(qE,q)

    # angle_errors  Small angle errors between two quaternions
    #       Computes the small angle errors associated with an estimated
    #       quaternion and its true value. The first order approximation is used
    #
    #       [dalpha] = angle_errors(qE, q)
    #       The inputs are:
    #           qE = [4 x m], Estimated Quaterion
    #           q = [4 x m], True Quaternion
    #
    #       The outputs are:
    #           dalpha = [3 x m], [roll pitch yaw]' errors expressed in radians
    #
    # Written by: Chris Nebelecky, converted to Julia by Stephen Gagnon

    if size(q,1) < 4 | size(q,1) > 4
        q = q'
    end
    if size(qE,1) < 4 | size(qE,1) > 4
        qE = qE'
    end
    if size(q,2) == 1
        q = repeat(q,1,size(qE,2))
    end

    dalpha = zeros(4,size(q,2))

    for i = 1:size(q,2)
        m = argmax(q[:,i])
        if !(sign(qE[m,i]) == sign(q[m,i]) == -1) | !(sign(qE[m,i]) == sign(q[m,i]) == 1)
                qE[:,i] = -qE[:,i]
        end
        dalpha[:,i] = 2*qprod(qE[:,i],qinv(q[:,i]))
    end

    #Take only first 3 rows
    dalpha = dalpha[1:3,:]
end



end
