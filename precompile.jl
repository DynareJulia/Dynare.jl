using PackageCompiler, Dynare

sysimage_path = "sys_dynare.so"

if Sys.iswindows()
    sysimage_path = "sys_dynare.dll"
elseif Sys.isapple()
    sysimage_path = "sys_dynare.dylib"
end
    
create_sysimage(:Dynare, sysimage_path=sysimage_path, precompile_execution_file="precompile_Dynare.jl")