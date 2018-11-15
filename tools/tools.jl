module LinearInterp
"""
linear interpolation without the hot air
I really need this function, and honestly I don't care that much
about its speed or flexibility. I just need sth that works for my
BL formulation
"""
struct Interp1D
    x::Vector{Real}
    y::Vector{Real}
    m::Vector{Real}
    b::Vector{Real}
end

"""
interp1D(x::Array{T, 1}, y::Array{T, 1}) where T<:Real
Linear interpollation function which returns a Interp1D struct
    ```
    x = collect(sort(rand(10)))
    y = sin.(10*x)
    xfine = collect(x[1]:0.01:x[end])
    intp1 = interp1D(x, y)
    yfine = intp1(xfine)
    ```
"""
function interp1D(x::Array{T, 1}, y::Array{T, 1}) where T<:Real
    m = (y[2:end]-y[1:end-1])./(x[2:end]-x[1:end-1])
    b = -m.*x[1:end-1]+y[1:end-1]
    return Interp1D(x, y, m, b)
end

function evaluate(intrp::Interp1D, xval::Real)
    if xval<intrp.x[1]
        return intrp.y[1]
    elseif xval>=intrp.x[end]
        return intrp.y[end]
    end
    for i in 1:length(intrp.x)-1
        if intrp.x[i]<=xval<intrp.x[i+1]
            return intrp.m[i]*xval+intrp.b[i]
        end
    end
end

function evaluate(intrp::Interp1D, xvec::Array{T, 1}) where T<:Real
    [evaluate(intrp, xvec[i]) for i in 1:length(xvec)]
end

(intrp1d::Interp1D)(x::Real) = evaluate(intrp1d, x)
(intrp1d::Interp1D)(x::Array{T, 1}) where T<:Real = evaluate(intrp1d, x)

end # module
