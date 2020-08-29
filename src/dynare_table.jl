using PrettyTables

function dynare_table(data, title, column_header, row_header, note;
                      fmt="%10.4f")
    title = ""
    if length(title) > 0
        pretty_table([title],
                     noheader = true,
                     tf = borderless,
                     highlighters = hl_row(1, crayon"bold"))
    end
    pretty_table(data,
                 noheader = true,
                 formatters = ft_printf(fmt, 1:size(data,1)+1),
                 cell_alignment = Dict( (1, i) => :c
                                        for i = 1 :size(data,1)+1),
                 highlighters = (hl_row(1, crayon"bold"),
                                 hl_col(1, crayon"bold")),
                 hlines = [:begin, 1, :end],
                 body_hlines_format = Tuple('─' for _ = 1:4),
                 vlines = [1])
    if length(note) > 0
        pretty_table([note],
                     noheader = true,
                     tf = borderless)
    end
    println(" ")
end

