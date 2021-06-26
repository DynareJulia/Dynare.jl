using JSON


function t1(mdofilename::String)
    modelstring::String = open(f -> read(f, String), modfilename*"/model/json/modfile.json")
    modeljson = JSON.parse(modelstring)
    return modeljson
end


    
