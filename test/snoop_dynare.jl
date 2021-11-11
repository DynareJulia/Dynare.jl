using Dynare
using SnoopCompile
options = Dynare.CommandLineOptions()
tinf = @snoopi_deep Dynare.parser("test/models/example1/example1", options)
itrigs = inference_triggers(tinf)
mtrigs = accumulate_by_source(Method, itrigs)
modtrigs = filtermod(Dynare, mtrigs)
