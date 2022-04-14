function qdq2w(q :: Vec, dq :: Vec)

    out = zeros(3,)
    out[1] = 2*( dq[1]*q[4] + dq[2]*q[3] - dq[3]*q[2] - dq[4]*q[1])
    out[2] = 2*(-dq[1]*q[3] + dq[2]*q[4] + dq[3]*q[1] - dq[4]*q[2])
    out[3] = 2*( dq[1]*q[2] - dq[2]*q[1] + dq[3]*q[4] - dq[4]*q[3])

    # E[1] = q[4]
    # E[2] = -q[3]
    # E[3] = q[2]
    #
    # E[4] = q[3]
    # E[5] = q[4]
    # E[6] = -q[1]
    #
    # E[7] = -q[2]
    # E[8] = q[1]
    # E[9] = q[4]
    #
    # E[10] = -q[1]
    # E[11] = -q[2]
    # E[12] = -q[3]

    return out #2 .* E*dq;
end

function qPropDisc(w,q,dt)
    wn = norm(w)
    phi = Array{Float64,1}(undef,3)
    phi[1] =  sin(.5*wn*dt)*w[1]/wn
    phi[2] =  sin(.5*wn*dt)*w[2]/wn
    phi[3] =  sin(.5*wn*dt)*w[3]/wn

    cwn = cos(.5*wn*dt)

    out = Array{Float64,1}(undef,4)
    out[1] =  q[1]*cwn + q[2]*phi[3] - q[3]*phi[2] + q[4]*phi[1]
    out[2] = -q[1]*phi[3] + q[2]*cwn + q[3]*phi[1] + q[4]*phi[2]
    out[3] =  q[1]*phi[2] - q[2]*phi[1] + q[3]*cwn + q[4]*phi[3]
    out[4] = -q[1]*phi[1] - q[2]*phi[2] - q[3]*phi[3] + q[4]*cwn

    # phi = sin(.5*wn)/wn .* w
    # O[1:3, 1:3] = diagm([cwn,cwn,cwn]) - crossMat(phi)
    # O[1:3, 4] = phi
    # O[4, 1:3] = -phi
    # O[4, 4] = cwn
    # return O*q

    return out
end

function crossMat(v :: Vec)
    M = zeros(3,3)
    M[2] = v[3]
    M[3] = -v[2]
    M[4] = -v[3]
    M[6] = v[1]
    M[7] = v[2]
    M[8] = -v[1]
    return M
end

function dAdp(att :: Vec)

    q = p2q(att)

    dA = Array{Array{Float64,1},2}(undef,3,3)

    dqdp_ = dqdp(q)
    dAdq_ = dAdq(q)

    for i = 1:9
        dA[i] = (dAdq_[i]'*dqdp_)[:]
    end
    return dA
end

function dAdq(q)
    dA = Array{Array{Float64,1},2}(undef,3,3)

    dA[1,1] = [2*q[1];-2*q[2];-2*q[3];2*q[4]]
    dA[1,2] = [2*q[2];2*q[1];2*q[4];2*q[3]]
    dA[1,3] = [2*q[3];-2*q[4];2*q[1];-2*q[2]]
    dA[2,1] = [2*q[2];2*q[1];-2*q[4];-2*q[3]]
    dA[2,2] = [-2*q[1];2*q[2];-2*q[3];2*q[4]]
    dA[2,3] = [2*q[4];2*q[3];2*q[2];2*q[1]]
    dA[3,1] = [2*q[3];2*q[4];2*q[1];2*q[2]]
    dA[3,2] = [-2*q[4];2*q[3];2*q[2];-2*q[1]]
    dA[3,3] = [-2*q[1];-2*q[2];2*q[3];2*q[4]]
    return dA
end

function dqdp(q)
    mat = q[1:3]*q[1:3]'-(1+q[4])*[1.0 0 0;0 1 0;0 0 1]
    return -[mat;(1+q[4])*q[1:3]']
end

function dDotdp(v1,v2,p)
    # d = Array{Float64,1}(undef,3)

    dAdp_ = dAdp(p)
    temp = Array{Float64,2}(undef,3,3)
    temp[1,:] = dAdp_[1,1]*v2[1] + dAdp_[1,2]*v2[2] + dAdp_[1,3]*v2[3]
    temp[2,:] = dAdp_[2,1]*v2[1] + dAdp_[2,2]*v2[2] + dAdp_[2,3]*v2[3]
    temp[3,:] = dAdp_[3,1]*v2[1] + dAdp_[3,2]*v2[2] + dAdp_[3,3]*v2[3]
    return (v1'*temp)'
end

function dDotdq(v1,v2,q)
    # d = Array{Float64,1}(undef,4)

    dAdq_ = dAdq(q)
    temp = Array{Float64,2}(undef,3,4)
    temp[1,:] = dAdq_[1,1]*v2[1] + dAdq_[1,2]*v2[2] + dAdq_[1,3]*v2[3]
    temp[2,:] = dAdq_[2,1]*v2[1] + dAdq_[2,2]*v2[2] + dAdq_[2,3]*v2[3]
    temp[3,:] = dAdq_[3,1]*v2[1] + dAdq_[3,2]*v2[2] + dAdq_[3,3]*v2[3]
    return (v1'*temp)'
end

function attDyn(t,x,J,L)

    Xi = Array{typeof(x[1]),2}(undef,4,3)
    Xi[1,1] = x[4]
    Xi[1,2] = -x[3]
    Xi[1,3] = x[2]
    Xi[2,1] = x[3]
    Xi[2,2] = x[4]
    Xi[2,3] = -x[1]
    Xi[3,1] = -x[2]
    Xi[3,2] = x[1]
    Xi[3,3] = x[4]
    Xi[4,:] = -view(x,1:3)


    dq = .5*Xi*view(x,5:7)
    dw = -inv(J)*(crossMat(view(x,5:7))*J*view(x,5:7) + L)

    dx = Array{typeof(x[1]),1}(undef,7)
    dx[1:4] = dq
    dx[5:7] = dw
    return dx
end

function attDyn(t,x,J,Jinv,L)

    Xi = Array{typeof(x[1]),2}(undef,4,3)
    Xi[1,1] = x[4]
    Xi[1,2] = -x[3]
    Xi[1,3] = x[2]
    Xi[2,1] = x[3]
    Xi[2,2] = x[4]
    Xi[2,3] = -x[1]
    Xi[3,1] = -x[2]
    Xi[3,2] = x[1]
    Xi[3,3] = x[4]
    Xi[4,:] = -view(x,1:3)


    dq = .5*Xi*view(x,5:7)
    dw = -Jinv*(crossMat(view(x,5:7))*J*view(x,5:7) + L)

    dx = Array{typeof(x[1]),1}(undef,7)
    dx[1:4] = dq
    dx[5:7] = dw
    return dx
end
