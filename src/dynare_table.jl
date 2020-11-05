using PrettyTables

function dynare_table_text(data, title, column_header, row_header, note, fmt)
    if length(title) > 0
        pretty_table([title],
                     noheader = true,
                     tf = borderless,
                     highlighters = (hl_row(1, crayon"bold"),
                                     hl_col(1, crayon"bold")),
                     backend = :text)
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
                 backend = :text)
    if length(note) > 0
        pretty_table([note],
                     noheader = true,
                     tf = borderless,
                     backend = :text)
    end
end

function dynare_table_latex(data, title, column_header, row_header, note, fmt)
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
                 backend = :latex)
    if length(note) > 0
        write(io, "\\begin{tablenotes}\n")
        write(io, "\\item $note\n")
        write(io, "\\end{tablenotes}\n")
    end
    write(io, "\\end{threeparttable}\n")
    write(io, "\\end{table}\n")
    return String(take!(io))
end

function dynare_table(data, title, column_header, row_header, note;
                      fmt="%10.4f", backend=:text)
    if backend == :text
        dynare_table_text(data, title, column_header, row_header, note, fmt)
    else
        dynare_table_latex(data, title, column_header, row_header, note, fmt)
    end
end
