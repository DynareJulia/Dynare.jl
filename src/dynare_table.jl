using PrettyTables

function dynare_table_text(data::AbstractVecOrMat{Any}, title::String, column_header::Vector, row_header::Vector, note::String; fmt="%10.4f")
    if length(title) > 0
        pretty_table([title],
                     noheader = true,
                     tf = tf_borderless,
                     highlighters = (hl_row(1, crayon"bold"),
                                     hl_col(1, crayon"bold")),
                     backend = Val(:text))
    end
    formatter = (v, i, j) -> (i > 1 && j > 1) ? round(v, digits=4) : v
    pretty_table(data,
                 noheader = true,
                 formatters = formatter,
                 cell_alignment = Dict( (1, i) => :c
                                        for i = 1 :size(data,1)+1),
                 highlighters = (hl_row(1, crayon"bold"),
                                 hl_col(1, crayon"bold")),
                 hlines = [:begin, 1, :end],
                 body_hlines_format = Tuple('─' for _ = 1:4),
                 vlines = [1],
                 backend = Val(:text))
    if length(note) > 0
        pretty_table([note],
                     noheader = true,
                     tf = tf_borderless,
                     backend = Val(:text))
    end
end

function dynare_table_latex(data::AbstractVecOrMat{Any}, title::String, column_header::Vector, row_header::Vector, note::String; fmt="%10.4f")
    io = IOBuffer()
    write(io, "\\begin{table}[h]\n")
    write(io, "\\centering\n")
    write(io, "\\begin{threeparttable}\n")
    if length(title) > 0
        write(io, "\\caption{$title}\n")
    end
    pretty_table(io,
                 data,
                 noheader = true,
 #                formatters = ft_printf(fmt, 1:size(data,1)+1),
#                 cell_alignment = Dict( (1, i) => :c
#                                        for i = 1 :size(data,1)+1),
#                 tf = latex_simple,
#                 highlighters = (LatexHighlighter((data, i, j) -> (i == 1), ["textbf"]),
#                                 LatexHighlighter((data, i, j) -> (j == 1), ["textbf"])),
#                 vlines = [1],
                 backend = Val(:latex))
    if length(note) > 0
        write(io, "\\begin{tablenotes}\n")
        write(io, "\\item $note\n")
        write(io, "\\end{tablenotes}\n")
    end
    write(io, "\\end{threeparttable}\n")
    write(io, "\\end{table}\n")
    return String(take!(io))
end

function dynare_table(data::AbstractVecOrMat{Any}, title::String, column_header::Vector, row_header::Vector, note::String;
                      fmt::String="%10.4f", backend=Val(:text))
    if backend == Val(:text)
        dynare_table_text(data, title, column_header, row_header, note, fmt=fmt)
    else
        dynare_table_latex(data, title, column_header, row_header, note, fmt=fmt)
    end
end
