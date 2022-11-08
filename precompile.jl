using PackageCompiler, Dynare
create_sysimage(:Dynare, sysimage_path="sys_dynare.so", precompile_execution_file="precompile_Dynare.jl")