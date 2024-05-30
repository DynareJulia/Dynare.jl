import Base.print
using Dates
using Tokenize

struct Graph
    filename::String
end

Base.print(io::IO, s::Graph) = print(
    io,
    "\\begin{centering}\\includegraphics[scale=0.6]{$(s.filename)}\\vspace{30pt}\\end{centering}\n",
)

struct Paragraph
    text::String
end

Base.print(io::IO, s::Paragraph) = print("$(s.text)\n")

struct Table
    string::String
    function Table(data, title, column_header, row_header, note)
        data = vcat(hcat("", column_header), hcat(row_header, data))
        string = dynare_table(data, title, note = note, backend = :latex)
        new(string)
    end
end

Base.print(io::IO, s::Table) = print(io, "$(s.string)\n")


struct Page
    sections::Vector{Any}
end

Page() = Page(Vector{Any}(undef, 0))
Page(s::String) = Page([Paragraph(s)])
Page(t::Table) = Page([t])
Page(g::Graph) = Page([g])

function print(io::IO, p::Page)
    for s in p.sections
        print(io, s)
    end
end

struct Report
    title::String
    subtitle::String
    pages::Vector{Page}
end

"""
  Report(title::String; subtitle::String = "")

initialize empty report

# Keyword arguments
 - `title::String`: Report title [required]
 - `subtitle::String`: Report subtitle
"""
function Report(title::String; subtitle::String = "")
    pages = Vector{Page}(undef, 0)
    Report(title, subtitle, pages)
end

"""
  add_page!(report::Report, page::Page)

adds a page to a report

# Keyword arguments
 - `report::Report`: report
 - `page::Page`: page to be added
"""
function add_page!(report::Report, page::Page)
    push!(report.pages, page)
end

"""
  add_graph!(page::Page, graph::Graph)

adds a graph to a page

# Keyword arguments
 - `page::Page`: page
 - `graph::Graph`: graph to be added
"""
function add_graph!(page::Page, graph::Graph)
    push!(page.sections, graph)
end

"""
  add_model!(page::Page; context::Context = context, lastline::Int = 0, format = 1)
adds the lines of a *.mod file to pages

# Keyword arguments
 - `page::Page`: page [required]
 - `context::Context`: context corresponding to the *.mod file (default: context)
 - `lastline::Int`: last line to be printed
 - `format::Int`: how to display parameter values
                  1: value is written after the parameter name
                  2: value is written below the parameter name
"""
function add_model!(page::Page; context::Context = context, lastline = 0, format = 1)
    model = modelprintout(
        context.modfileinfo.modfilepath,
        context.symboltable,
        context.work.params,
        sqrt.(diag(context.models[1].Sigma_e)),
        lastline,
        format,
    )
    push!(page.sections, model)
end

"""
  add_paragraph!(page::Page, paragraph::String)

adds a graph to a page

# Keyword arguments
 - `page::Page`: page
 - `paragraph::String`: paragraph to be added
"""
function add_paragraph!(page::Page, paragraph::String)
    push!(page.sections, paragraph)
end

"""
  add_table!(page::Page, table::Table)

adds a graph to a page
 
# Keyword arguments
 - `page::Page`: page
 - `table::Table`: table to be added
"""
function add_table!(page::Page, table::Table)
    push!(page.sections, table)
end

"""
  print(report::Report; texfilename::String = "report.tex")

creates a PDF file with the report

Keyword arguments
  - `report::Report`: report name
  - `texfilename::String`: name of the PDF file and of the LaTex file created to produce it (optional, default = `report.tex`)
"""
function print(report::Report; texfilename::String = "report.tex")
    open(texfilename, "w") do io
        print(io, "\\documentclass{report}\n")
        print(io, "\\usepackage{graphicx}\n")
        print(io, "\\usepackage{xcolor}\n")
        print(io, "\\usepackage{stackrel}\n")
        print(io, "\\usepackage{threeparttable}\n")
        print(io, "\\usepackage{listings}\n")
        print(io, "\\lstset{numbers=left}\n")
        print(io, "\\begin{document}\n")
        print(io, "\\vspace*{0.2\\textheight}\n")
        print(io, "\\begin{center}\n")
        print(io, "\\Large\\textbf{$(report.title)}\\\\\n")
        print(io, "\\medskip\n")
        if length(report.subtitle) > 0
            print(io, "\\large $(report.subtitle)\n")
            print(io, "\\medskip\n")
        end
        print(io, "\\end{center}\n")
        print(io, "$(Dates.now())\\\\\n")
        print(io, "\\clearpage\n")
        for (i, page) in enumerate(report.pages)
            print(io, page)
            if i < length(report.pages)
                print(io, "\\newpage\n")
            end
        end
        print(io, "\\end{document}")
    end

    latex = `pdflatex $texfilename`
    run(latex)
    return nothing
end

function modelprintout(
    modname::String,
    symboltable::SymbolTable,
    parameters_value::Vector{Float64},
    sd::Vector{Float64},
    lastline::Integer,
    format::Integer,
)
    out = IOBuffer()
    if lastline > 0
        print(
            out,
            "\\begin{lstlisting}[escapechar = |, breaklines = true, lastline = $(lastline)]\n",
        )
    else
        print(out, "\\begin{lstlisting}[escapechar = |, breaklines = true]\n")
    end
    elements = []
    linenumber = 1
    model_mode = false
    symbols = keys(symboltable)
    for line in readlines(modname * ".mod")
        for token in tokenize(line)
            stringtoken = Tokens.untokenize(token)
            kind = Tokens.kind(token)
            if kind == Tokens.IDENTIFIER
                if stringtoken == "model"
                    model_mode = true
                elseif stringtoken == "end"
                    model_mode = false
                elseif model_mode && stringtoken in symbols
                    k = symboltable[stringtoken].orderintype
                    if is_parameter(stringtoken, symboltable)
                        if format == 1
                            stringtoken = "$(stringtoken)|\$ {\\color{red}\\scriptstyle <$(parameters_value[k])>}\$|"
                        elseif format == 2
                            stringtoken = "|\$ \\stackrel[($(parameters_value[k]))]{}{\\hbox{$(stringtoken)}}\$|"
                        else
                            error("Wrong format value for printing model")
                        end
                    elseif is_exogenous(stringtoken, symboltable)
                        if format == 1
                            stringtoken = "$(stringtoken)|\$ {\\color{red}\\scriptstyle <\\sigma = $(sd[k])>}\$|"
                        end 
                    end
                end
                push!(elements, stringtoken)
            elseif kind == Tokens.WHITESPACE
                subelements = []
                if occursin("\r\n", stringtoken)
                    subelements = split(stringtoken, "\r\n")
                elseif '\n' in stringtoken
                    subelements = split(stringtoken, '\n')
                end
                if length(subelements) > 0
                    for (i, se) in enumerate(subelements)
                        if i == 1
                            continue
                        end
                        push!(elements, se)
                        println(out, join(elements))
                        elements = []
                        linenumber += 1
                    end
                else
                    push!(elements, stringtoken)
                end
            else
                push!(elements, stringtoken)
            end
        end
    end
    print(out, "\\end{lstlisting}\n")
    return String(take!(out))
end
