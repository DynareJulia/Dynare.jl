function sparse!(I::AbstractVector{Ti}, J::AbstractVector{Ti},
        V::AbstractVector{Tv}, m::Integer, n::Integer, combine, klasttouch::Vector{Tj},
        csrrowptr::Vector{Tj}, csrcolval::Vector{Ti}, csrnzval::Vector{Tv},
        csccolptr::Vector{Ti}, cscrowval::Vector{Ti}, cscnzval::Vector{Tv}) where {Tv,Ti<:Integer,Tj<:Integer}

#    require_one_based_indexing(I, J, V)
#    sparse_check_Ti(m, n, Ti)
#    sparse_check_length("I", I, 0, Tj)
    # Compute the CSR form's row counts and store them shifted forward by one in csrrowptr
    fill!(csrrowptr, Tj(0))
    coolen = length(I)
    min(length(J), length(V)) >= coolen || throw(ArgumentError("J and V need length >= length(I) = $coolen"))
    @inbounds for k in 1:coolen
        Ik = I[k]
        if 1 > Ik || m < Ik
            throw(ArgumentError("row indices I[k] must satisfy 1 <= I[k] <= m"))
        end
        csrrowptr[Ik+1] += Tj(1)
    end
    @show csrrowptr
    # Compute the CSR form's rowptrs and store them shifted forward by one in csrrowptr
    countsum = Tj(1)
    csrrowptr[1] = Tj(1)
    @inbounds for i in 2:(m+1)
        overwritten = csrrowptr[i]
        csrrowptr[i] = countsum
        countsum += overwritten
    end
    @show csrrowptr

    # Counting-sort the column and nonzero values from J and V into csrcolval and csrnzval
    # Tracking write positions in csrrowptr corrects the row pointers
    @inbounds for k in 1:coolen
        Ik, Jk = I[k], J[k]
        if Ti(1) > Jk || Ti(n) < Jk
            throw(ArgumentError("column indices J[k] must satisfy 1 <= J[k] <= n"))
        end
        csrk = csrrowptr[Ik+1]
        @assert csrk >= Tj(1) "index into csrcolval exceeds typemax(Ti)"
        csrrowptr[Ik+1] = csrk + Tj(1)
        csrcolval[csrk] = Jk
        csrnzval[csrk] = V[k]
    end
    @show csrrowptr
    # This completes the unsorted-row, has-repeats CSR form's construction

    # Sweep through the CSR form, simultaneously (1) calculating the CSC form's column
    # counts and storing them shifted forward by one in csccolptr; (2) detecting repeated
    # entries; and (3) repacking the CSR form with the repeated entries combined.
    #
    # Minimizing extraneous communication and nonlocality of reference, primarily by using
    # only a single auxiliary array in this step, is the key to this method's performance.
    fill!(csccolptr, Ti(0))
    fill!(klasttouch, Tj(0))
    writek = Tj(1)
    newcsrrowptri = Ti(1)
    origcsrrowptri = Tj(1)
    origcsrrowptrip1 = csrrowptr[2]
    @inbounds for i in 1:m
        for readk in origcsrrowptri:(origcsrrowptrip1-Tj(1))
            j = csrcolval[readk]
            if klasttouch[j] < newcsrrowptri
                klasttouch[j] = writek
                if writek != readk
                    csrcolval[writek] = j
                    csrnzval[writek] = csrnzval[readk]
                end
                writek += Tj(1)
                csccolptr[j+1] += Ti(1)
            else
                klt = klasttouch[j]
                csrnzval[klt] = combine(csrnzval[klt], csrnzval[readk])
            end
        end
        newcsrrowptri = writek
        origcsrrowptri = origcsrrowptrip1
        origcsrrowptrip1 != writek && (csrrowptr[i+1] = writek)
        i < m && (origcsrrowptrip1 = csrrowptr[i+2])
    end

    # Compute the CSC form's colptrs and store them shifted forward by one in csccolptr
    countsum = Tj(1)
    csccolptr[1] = Ti(1)
    @inbounds for j in 2:(n+1)
        overwritten = csccolptr[j]
        csccolptr[j] = countsum
        countsum += overwritten
        Base.hastypemax(Ti) && (countsum <= typemax(Ti) || throw(ArgumentError("more than typemax(Ti)-1 == $(typemax(Ti)-1) entries")))
    end

    # Now knowing the CSC form's entry count, resize cscrowval and cscnzval if necessary
    cscnnz = countsum - Tj(1)
    resize!(cscrowval, cscnnz)
    resize!(cscnzval, cscnnz)

    # Finally counting-sort the row and nonzero values from the CSR form into cscrowval and
    # cscnzval. Tracking write positions in csccolptr corrects the column pointers.
    @show m
    @show csrrowptr
    @inbounds for i in 1:m
        for csrk in csrrowptr[i]:(csrrowptr[i+1]-Tj(1))
            j = csrcolval[csrk]
            x = csrnzval[csrk]
            csck = csccolptr[j+1]
            #@show j+1, csck + 1
            csccolptr[j+1] = csck + Ti(1)
            cscrowval[csck] = i
            cscnzval[csck] = x
        end
    end
    @show csccolptr
    #(length(colptr) == n + 1 && colptr[end] - 1 == length(rowval) == length(nzval))
    @show n
    @show length(csccolptr)
    @show csccolptr[end] - 1
    @show length(cscrowval)
    @show length(cscnzval)
    SparseMatrixCSC(m, n, csccolptr, cscrowval, cscnzval)
end

I = [1, 3, 4, 1, 2, 4, 2, 4, 3, 5, 1, 3, 2, 4, 6, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 9, 10, 11, 12, 11, 12, 7, 9, 10, 7, 8, 10, 8, 10, 9, 11, 7, 9, 8, 10, 12, 8, 8, 8, 1, 1, 15, 16, 17, 18, 17, 18, 13, 15, 16, 13, 14, 16, 14, 16, 15, 17, 13, 15, 14, 16, 18, 14, 14, 14, 1, 1, 21, 22, 23, 24, 23, 24, 19, 21, 22, 19, 20, 22, 20, 22, 21, 23, 19, 21, 20, 22, 24, 20, 20, 20]
J = [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 8, 12, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 4, 4, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 12, 13, 14, 18, 1, 1, 9, 9, 10, 10, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 18, 19, 20, 24, 1, 1, 15, 15, 16, 16, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 24, 21, 22, 24]
V = [-0.64, 1.0, -1.0, 0.8606811145510832, -13.792569659442691, 1.0, 0.03475000000000006, 1.0, -1.0806825309567203, 1.0, 2.370597639417811, -2.3705976394178108, -11.083604432603583, -0.2770901108150896, 1.0, -0.3564, 13.792569659442691, 10.698449178570607, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.03510101010101016, -0.975, -0.95, -0.025, -0.025, -0.95, -0.64, 1.0, -1.0, 0.8606811145510832, -13.792569659442691, 1.0, 0.03475000000000006, 1.0, -1.0806825309567203, 1.0, 2.370597639417811, -2.3705976394178108, -11.083604432603583, -0.2770901108150896, 1.0, -0.3564, 13.792569659442691, 10.698449178570607, 0.0, 0.0, -0.03510101010101016, -0.975, -0.95, -0.025, -0.025, -0.95, -0.64, 1.0, -1.0, 0.8606811145510832, -13.792569659442691, 1.0, 0.03475000000000006, 1.0, -1.0806825309567203, 1.0, 2.370597639417811, -2.3705976394178108, -11.083604432603583, -0.2770901108150896, 1.0, -0.3564, 13.792569659442691, 10.698449178570607, 0.0, 0.0, -0.03510101010101016, -0.975, -0.95, -0.025, -0.025, -0.95, -0.64, 1.0, -1.0, 0.8606811145510832, -13.792569659442691, 1.0, 0.03475000000000006, 1.0, -1.0806825309567203, 1.0, 2.370597639417811, -2.3705976394178108, -11.083604432603583, -0.2770901108150896, 1.0, 0.5296781221467004, 5.468940420105765, 5.468940420105738]
@assert length(I) == length(J) == length(V)
n = 24
nz = length(V)
klasstouch = Vector{Int64}(undef, nz)
colptr = Vector{Int64}(undef, n + 1)
rowptr = Vector{Int64}(undef, n + 1)
colval = Vector{Int64}(undef, nz)
nzval = Vector{Float64}(undef, nz)

sparse!(
    I,
    J,
    V,
    n,
    n,
    (x, y) -> x + y,
    klasstouch,
    rowptr,
    colval,
    nzval,
    colptr,
    J,
    V,
)
