import SparseArrays

struct SparseStorage{Tv,Ti<:Integer}
    I::Vector{Ti}
    J::Vector{Ti}
    V::Vector{Tv}
    klasttouch::Vector{Ti}
    csrrowptr::Vector{Ti}
    csrcolval::Vector{Ti}
    csrnzval::Vector{Tv}
    function SparseStorage{Tv,Ti}(
        m::Integer,
        n::Integer,
        nz::Integer,
    ) where {Tv,Ti<:Integer}
        I = Vector{Ti}(undef, nz)
        J = Vector{Ti}(undef, nz)
        V = Vector{Tv}(undef, nz)
        klasttouch = Vector{Ti}(undef, n)
        csrrowptr = Vector{Ti}(undef, m + 1)
        csrcolval = Vector{Ti}(undef, nz)
        csrnzval = Vector{Tv}(undef, nz)
        new(I, J, V, klasttouch, csrrowptr, csrcolval, csrnzval)
    end
end

SparseStorage(m::Integer, n::Integer, nz::Integer) = SparseStorage{Float64,Int64}(m, n, nz)

function SparseStorage(
    m::Integer,
    n::Integer,
    nz::Integer,
    I::Vector{Ti},
    J::Vector{Ti},
    V::Vector{Tv},
) where {Tv,Ti}
    nz = length(I)
    klasttouch = Vector{Ti}(undef, n)
    csrrowptr = Vector{Ti}(undef, m + 1)
    csrcolval = Vector{Ti}(undef, nz)
    csrnzval = Vector{Tv}(undef, nz)
    SparseStorage(I, J, V, klasstouch, csrrowptr, csrcolval, csrnzval)
end

SparseStorage(
    m::Integer,
    n::Integer,
    nz::Integer,
    I::Vector{Int64},
    J::Vector{Int64},
    V::Vector{Float64},
) = SparseStorage{Float64,Int64}(
    m::Integer,
    n::Integer,
    nz::Integer,
    I::Vector{Int64},
    J::Vector{Int64},
    V::Vector{Float64},
)

sparse!(m::Integer, n::Integer, nz::Integer, SS::SparseStorage) = SparseArrays.sparse!(
    SS.I,
    SS.J,
    SS.V,
    m,
    n,
    (x, y) -> x + y,
    SS.klasttouch,
    SS.csrrowptr,
    SS.csrcolval,
    SS.csrnzval,
)
