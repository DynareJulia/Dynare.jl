DYNARE_ROOT = "/data/projects/dynare/git/preprocessor/src/dynare_m"

dynare_preprocess(modfilename) = run(`$DYNARE_ROOT $modfilename language=julia output=third json=compute`) 

