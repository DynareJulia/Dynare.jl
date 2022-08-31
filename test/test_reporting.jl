include("../src/dynare_table.jl")
include("../src/reporting/report.jl")
using Random
using Test

page1 = Page(
    "This is the first paragraph\\begin{itemize}\\item Item 1 \\item Item 2 \\item Item 3\\end{itemize}",
)
page2 = Page(
    "This is the second paragraph\\begin{itemize}\\item Item 1 \\item Item 2 \\item Item 3\\end{itemize}",
)

Random.seed!(111)
data = randn(4, 3)
title = "Table title"
cheaders = ["col1", "col2", "col3"]
rheaders = ["row1", "row2", "row3", "row4"]
note = "Table note"
table1 = Table(data, title, cheaders, rheaders, note)

target = raw"""
\begin{table}[h]
\centering
\begin{threeparttable}
\caption{Table title}
\begin{tabular}{r|rr}
\hline
\multicolumn{1}{c|}{\textbf{    0.7696}} & \multicolumn{1}{c}{\textbf{   -0.7842}} & \multicolumn{1}{c}{\textbf{    0.1743}} \\
\textbf{   -1.1348} &    -0.2939 &     0.4024 \\
\textbf{    1.3340} &     1.1560 &     0.6967 \\
\textbf{   -0.1824} &    -1.0381 &     1.6320 \\\hline
\end{tabular}
\begin{tablenotes}
\item Table note
\end{tablenotes}
\end{threeparttable}
\end{table}
"""
@test table1.string == target

report = Report("Report title", subtitle = "Subtitle")

add_table!(page2, table1)
add_page!(report, page1)
add_page!(report, page2)

print(report)
