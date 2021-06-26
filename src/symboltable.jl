for typ in instances(SymbolType)
    for f in fieldnames(DynareSymbol)
        s = Symbol("get_$(lowercase(string(typ)))_$(f)")
        @eval begin
            function $s(symboltable::SymbolTable)
                symbols = collect(values(symboltable))
                subset = filter(s -> s.symboltype::SymbolType == $typ, symbols)
                sorted_index = sortperm(subset, by = v -> v.orderintype)
                names = [s.$f for s in subset[sorted_index]]
                return names
            end
        end
    end
    s = Symbol("is_$(lowercase(string(typ)))") 
    @eval begin
        function $s(name::String, symboltable::SymbolTable)
            return (symboltable[name].symboltype == $typ)
        end
    end
end
