include("insertLine.jl");

function ratToStr(nr)
    num = numerator(nr);
    if denominator(nr) != 1
        str = "\$";
        str *= replace(string(sign(num)), r"\d" => s"");
        str *= "\\dfrac{"
        str *= string(abs(num));
        str *= "}{";
        str *= string(denominator(nr));
        str *= "}\$";
    else
        str = "\$";
        str *= string(num);
        str *= "\$";
    end

    return str;
end

function subscript(str)
    return "\$" * replace(str, r"(\d+)" => s"_\1") * "\$";
end

function genTableau(A, b, cj, x, xB, cB, z, zc, ratio)
    dims = getDims(A);
    m, mn = dims.m, dims.mn;
    isRatCol = (minimum(zc) < 0 && minimum(b) >= 0);
    isRatRow = (minimum(b) < 0);
    file = "tableaux-display.tex";
    lineNr = lineCalc(file);

    # Begin table
    str = "\\begin{table}[ht]\n";
    str *= "\\bgroup\n\\def\\arraystretch{3.0}\n";
    str *= "\\begin{tabular}";

    # First row
    str *= "{";
    str *= "|c"^(mn + 3 + isRatCol);
    str *= "|}\n";
    str *= "\\hline\n";
    str *= "& \$c_j\$ & ";
    for i in cj
        str *= ratToStr(i);
        str *= " & ";
    end
    if isRatCol
        str *= " &  ";
    end
    str *= "\\\\";
    str *= "\\hline\n";

    # Second row
    str *= "\$c_{\\mathbf{B}}\$";
    str *= " & ";
    str *= "\$x_{\\mathbf{B}}\$";
    str *= " & ";
    for i in x
        str *= subscript(i);
        str *= " & ";
    end
    str *= "\$\\mathbf{b}\$";
    if isRatCol
        str *= " & Ratio";
    end
    str *= "\\\\";
    str *= "\\hline\n";

    # Intermediate rows
    for i in 1:m
        str *= ratToStr(cB[i]);
        str *= " & ";
        str *= subscript(xB[i]);
        str *= " & ";
        for j in 1:mn
            str *=  ratToStr(A[i,j]);
            str *= " & ";
        end
        str *= ratToStr(b[i]);
        if isRatCol
            str *= " & ";
            str *= ratToStr(ratio[i]);
        end
        str *= "\\\\";
        str *= "\\hline\n";
    end

    # zj row
    str *= "\\multirow{";
    mergRNo = 2 + isRatRow;
    str *= string(mergRNo);
    str *= "}{*}{}";
    str *= " & ";
    str *= "\$z_j\$";
    str *= " & ";
    for i in 1:mn
        str *= ratToStr(z[i]);
        str *= " & ";
    end
    str *= "\\multirow{2}{*}{";
    str *= ratToStr(z[mn+1]);
    str *= "}";
    str *= " & "^isRatCol;
    str *= "\\\\";

    # zj-cj row
    str *= " & ";
    str *= "\$z_j - c_j\$";
    str *= " & ";
    for i in 1:mn
        str *= ratToStr(zc[i]);
        str *= " & ";
    end
    str *= " & "^isRatCol;
    str *= "\\\\";
    str *= "\\hline";

    # Ratio row
    if isRatRow
        str *= " & ";
        str *= "Ratio";
        str *= " & ";
        for i in 1:mn
            if (ratio[i] != Inf)
                str *= ratToStr(ratio[i]);
                str *= " & ";
            else
                str *= " & ";
            end
        end
        str *= "\\\\";
        str *= "\\hline";
    end
        

    # End table
    str *= "\n\\end{tabular}\n";
    str *= "\\egroup\n";
    str *= "\\end{table}\n";

    # Create file copy and write str to it at lineNr
    insertLine(file, str, lineNr);
    # buildTex(file);
end

function lineCalc(file)
    output = read(`wc -l $file`, String)
    lineNr = parse(Int64, replace(output, r"\s.*" => s"")) - 1;
    return lineNr;
end

function rowOps(A, pivColIdx, pivRowIdx)
    m = getDims(A).m;
    file = "tableaux-display.tex";
    lineNr = lineCalc(file);
    if (pivColIdx != 0)
        pivCol = A[:, pivColIdx];
        pivEl = pivCol[pivRowIdx];
    else
        return;
    end
    
    str = "";

    for i in 1:m
        if pivRowIdx == i
            str *= "\$ \\dfrac{";
            str *= ratToStr(1/pivEl);
            str *= "R_{";
            str *= string(i);
            str *= "} \\rightarrow R_{";
            str *= string(i);
            str *= "}^{'}\$\n";
        elseif pivCol[i] !=0
            str *= "\$ R_{";
            str *= string(i);
            str *= "} - ";
            str *= ratToStr(1/pivCol[i]);
            str *= "R_{";
            str *= string(pivRowIdx);
            str *= "}^{'}";
        end
    end

    insertLine(file, str, lineNr);
end

function buildTex(file)
    pdf = replace(file, ".tex" => ".pdf")
    run(`pdflatex -synctex=1 -interaction=nonstopmode $file`);
    run(`okular $pdf`);
end