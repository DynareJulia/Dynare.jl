for typ in instances(SymbolType)
    s = Symbol("get_$(lowercase(string(typ)))")
    @eval begin
        function $s(symboltable::SymbolTable)
            subset = [
                (sym[1], sym[2].orderintype) for
                sym in context.symboltable if sym[2].symboltype == $typ
            ]
            return [sym[1] for sym in sort(subset, by = sym -> sym[2])]
        end
    end

    for f in fieldnames(DynareSymbol)
        s = Symbol("get_$(lowercase(string(typ)))_$(f)")
        @eval begin
            function $s(symboltable::SymbolTable)
                symbols = collect(values(symboltable))
                subset = filter(s -> s.symboltype::SymbolType == $typ, symbols)
                sorted_index = sortperm(subset, by = v -> v.orderintype)
                names = [sym.$f for sym in subset[sorted_index]]
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
