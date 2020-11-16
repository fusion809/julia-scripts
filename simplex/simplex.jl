include("tableau.jl");

function minNonNeg(arr)
    minVal = Inf;
    index = 0;
    for i in 1:length(arr)
        if 0 <= arr[i] < minVal
            @inbounds minVal = arr[i];
            @inbounds index = i;
        end
    end

    return index;
end

function minRatNegDen(pivRow, zc)
    minVal = Inf;
    index = 0;
    ratio = abs(zc./pivRow);
    for i in 1:length(pivRow)
        if pivRow[i] < 0 && ratio[i] < minVal
            @inbounds minVal = ratio[i];
            @inbounds index = i;
        elseif pivRow[i] >= 0
            @inbounds ratio[i] = Inf;
        end
    end

    return index, ratio;
end
            

function simplex(A, b, cj, x, xB)
    cB, z, zc = calcEntries(A, b, cj, x, xB);
    dims = getDims(A);
    ratio = 0;
    pivColIdx = 0;
    pivRowIdx = 0;
    m = dims.m;
    mn = dims.mn;

    if minimum(b) >= 0 && minimum(zc) < 0
        pivColIdx = argmin(zc)[2];
        pivCol = A[:, pivColIdx];
        leavingVar = x[pivColIdx];
        ratio = b.//pivCol;

        pivRowIdx = minNonNeg(ratio);
        pivRow = A[pivRowIdx, :];
        enteringVar = xB[pivRowIdx];

        pivEl = A[pivRowIdx, pivColIdx];

        A[pivRowIdx, :] /= pivEl;
        b[pivRowIdx] /= pivEl;

        # Update pivot row and column
        pivRow = A[pivRowIdx, :];
        pivCol = A[:, pivColIdx];
        
        for i in 1:m
            if pivRowIdx != i
                @inbounds A[i,:] -= pivCol[i]*pivRow;
                @inbounds b[i] -= pivCol[i]*b[pivRowIdx];
            end
        end

        xB[pivRowIdx] = x[pivColIdx];
    elseif minimum(b) < 0
        pivRowIdx = argmin(b);
        pivRow = A[pivRowIdx, :];
        
        pivColIdx, ratio = minRatNegDen(pivRow, zc);
        pivCol = A[:, pivColIdx];
        
        pivEl = pivCol[pivRowIdx];

        A[pivRowIdx, :] /= pivEl;
        b[pivRowIdx] /= pivEl;

        pivRow = A[pivRowIdx, :];
        pivCol = A[:, pivColIdx];

        for i in 1:m
            if pivRowIdx != i
                @inbounds A[i, :] -= pivCol[i] * pivRow;
                @inbounds b[i] -= pivCol[i] * b[pivRowIdx];
            end
        end

        xB[pivRowIdx] = x[pivColIdx];
    end

    return A, b, cj, x, xB, cB, z, zc, ratio, pivColIdx, pivRowIdx;
end

function simplexIterator(A, b, cj, x, xB)
    cB, z, zc = calcEntries(A, b, cj, x, xB);
    ratio, pivColIdx, pivRowIdx = 0, 0, 0;
    cp("simplex.tex", "tableaux-display.tex", force=true);

    while minimum(zc) < 0 || minimum(b) < 0
        A, b, cj, x, xB, cB, z, zc, ratio, pivColIdx, pivRowIdx = simplex(A, b, cj, x, xB);
        genTableau(A, b, cj, x, xB, cB, z, zc, ratio);
        rowOps(A, pivColIdx, pivRowIdx);
    end

    buildTex("tableaux-display.tex");
    return A, b, cj, x, xB, cB, z, zc, ratio, pivColIdx, pivRowIdx;
end

function basisIndex(x, xB)
    m = length(xB);
    mn = length(x);
    loc = zeros(Int64, m);

    for i in 1:m
        for j in 1:mn
            if xB[i] == x[j]
                @inbounds loc[i] = j;
                break;
            end
        end
    end

    return loc;
end

function getDims(A)
    dims = size(A);
    return (m = dims[1], mn = dims[2]);
end

function calcEntries(A, b, cj, x, xB)
    dims = getDims(A);
    m = dims.m;
    mn = dims.mn;
    zc = zeros(Rational{Int64}, 1, mn);
    loc = basisIndex(x, xB);
    cB = cj[loc];

    if (m != length(xB))
        throw("Length of xB does not match the number of rows in A");
        return;
    end

    zj = cB'*A;
    zc = zj - cj;
    zj = reduce(vcat, zj);
    Base.push!(zj, (cB'b)[1]);
    zj = reduce(hcat, zj);

    return cB, zj, zc;
end