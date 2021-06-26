using SnoopCompile
using Dynare
using AbstractTrees
options = Dynare.CommandLineOptions()
tinf = @snoopi_deep Dynare.parser("test/models/example1/example1", options)
itrigs = inference_triggers(tinf)
mtrigs = accumulate_by_source(Method, itrigs)
modtrigs = filtermod(Dynare, mtrigs)
itree = trigger_tree(itrigs)
print_tree(itrigs)
