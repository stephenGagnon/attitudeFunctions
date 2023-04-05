function qdq2w(q :: Vec{T}, dq :: Vec) where {T <: Real}
    qout = zeros(T, 3)
    return qdq2w(q :: Vec{T}, dq :: Vec, qout :: Vec)
end

function qdq2w(q :: Vec, dq :: Vec, qout :: Vec)
    qout[1] = 2*( dq[1]*q[4] + dq[2]*q[3] - dq[3]*q[2] - dq[4]*q[1])
    qout[2] = 2*(-dq[1]*q[3] + dq[2]*q[4] + dq[3]*q[1] - dq[4]*q[2])
    qout[3] = 2*( dq[1]*q[2] - dq[2]*q[1] + dq[3]*q[4] - dq[4]*q[3])

    return qout #2 .* E*dq;
end

function qPropDisc(w, q :: Vec{T}, dt) where {T <: Real}
    phi = Array{T,1}(undef,3)
    out = Array{T,1}(undef,4)

    return qPropDisc(w, q , dt, phi, out)
end

function qPropDiscAlt(w, q :: Vec{T}, dt) where {T <: Real}
    out = Array{T,1}(undef,4)
    return qPropDiscAlt(w, q , dt, out)
end

function qPropDisc(w, q , dt, phi, qout)
    wn = norm(w)
    # phi = Array{T,1}(undef,3)
    phi[1] =  sin(.5*wn*dt)*w[1]/wn
    phi[2] =  sin(.5*wn*dt)*w[2]/wn
    phi[3] =  sin(.5*wn*dt)*w[3]/wn

    cwn = cos(.5*wn*dt)

    # out = Array{T,1}(undef,4)
    qout[1] =  q[1]*cwn + q[2]*phi[3] - q[3]*phi[2] + q[4]*phi[1]
    qout[2] = -q[1]*phi[3] + q[2]*cwn + q[3]*phi[1] + q[4]*phi[2]
    qout[3] =  q[1]*phi[2] - q[2]*phi[1] + q[3]*cwn + q[4]*phi[3]
    qout[4] = -q[1]*phi[1] - q[2]*phi[2] - q[3]*phi[3] + q[4]*cwn
    return qout
end

function qPropDiscAlt(w,q,dt,qout)
    qout[1] = q[1] + dt*(q[2]*w[3] - q[3]*w[2] + q[4]*w[1])
    qout[2] = q[2] + dt*(-q[1]*w[3] + q[3]*w[1] + q[4]*w[2])
    qout[3] = q[3] + dt*(q[1]*w[2] - q[2]*w[1] + q[4]*w[3])
    qout[4] = q[4] + dt*(q[1]*w[1] + q[2]*w[2] + q[3]*w[3])
    return qout
end

function qPropDisc!(w, q , dt, phi, qtemp)
    wn = norm(w)
    # phi = Array{T,1}(undef,3)
    phi[1] =  sin(.5*wn*dt)*w[1]/wn
    phi[2] =  sin(.5*wn*dt)*w[2]/wn
    phi[3] =  sin(.5*wn*dt)*w[3]/wn

    cwn = cos(.5*wn*dt)

    # out = Array{T,1}(undef,4)
    qtemp[1] =  q[1]*cwn + q[2]*phi[3] - q[3]*phi[2] + q[4]*phi[1]
    qtemp[2] = -q[1] * phi[3] + q[2] * cwn + q[3] * phi[1] + q[4] * ph
    qtemp[3] = q[1] * phi[2] - q[2] * phi[1] + q[3] * cwn + q[4] * phi
    qtemp[4] = -q[1] * phi[1] - q[2] * phi[2] - q[3] * phi[3] + q[4]

    return q[:] = qtemp
end

function crossMat(v :: Vec{T}) where {T <: Real}
    M = zeros(T, 3, 3)
    M[2] = v[3]
    M[3] = -v[2]
    M[4] = -v[3]
    M[6] = v[1]
    M[7] = v[2]
    M[8] = -v[1]
    return M
end

function dAdtht(tht)
    return [-sin(tht) -cos(tht); cos(tht) -sin(tht)]
end

function dAdp(att :: Vec{T}) where {T <: Real}

    q = p2q(att)

    dA = Array{Array{T,1},2}(undef,3,3)

    dqdp_ = dqdp(q)
    dAdq_ = dAdq(q)

    for i = 1:9
        dA[i] = (dAdq_[i]'*dqdp_)[:]
    end
    return dA
end

function dAdq(q :: Vec{T}) where {T <: Real}
    dA = Array{Array{T,1},2}(undef,3,3)

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

function dDotdp(v1, v2, p :: Vec{T}) where {T <: Real}
    # d = Array{Float64,1}(undef,3)

    dAdp_ = dAdp(p)
    temp = Array{T,2}(undef,3,3)
    temp[1,:] = dAdp_[1,1]*v2[1] + dAdp_[1,2]*v2[2] + dAdp_[1,3]*v2[3]
    temp[2,:] = dAdp_[2,1]*v2[1] + dAdp_[2,2]*v2[2] + dAdp_[2,3]*v2[3]
    temp[3,:] = dAdp_[3,1]*v2[1] + dAdp_[3,2]*v2[2] + dAdp_[3,3]*v2[3]
    return (v1'*temp)'
end

function dDotdq(v1,v2,q :: Vec{T}) where {T <: Real}
    # d = Array{Float64,1}(undef,4)

    dAdq_ = dAdq(q)
    temp = Array{T,2}(undef,3,4)
    temp[1,:] = dAdq_[1,1]*v2[1] + dAdq_[1,2]*v2[2] + dAdq_[1,3]*v2[3]
    temp[2,:] = dAdq_[2,1]*v2[1] + dAdq_[2,2]*v2[2] + dAdq_[2,3]*v2[3]
    temp[3,:] = dAdq_[3,1]*v2[1] + dAdq_[3,2]*v2[2] + dAdq_[3,3]*v2[3]
    return (v1'*temp)'
end

function attDyn(t, x :: Vec{T}, J, L) where {T <: Real}

    Xi = Array{T,2}(undef,4,3)
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

    dx = Array{T,1}(undef,7)
    dx[1:4] = dq
    dx[5:7] = dw
    return dx
end

function attDyn(t, x :: Vec{T}, J, Jinv, L)  where {T <: Real}

    Xi = Array{T,2}(undef,4,3)
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

    dx = Array{T,1}(undef,7)
    dx[1:4] = dq
    dx[5:7] = dw
    return dx
end
