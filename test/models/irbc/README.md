# Model irbc.mod

Model irbc.mod is a multicountry model with a user defined number of
countries. Each country model has 7 variables. There one common shock
and one shock per country.

Similar equations are repeated for each country using Dynare
macrolanguage (@#for loops)

All countries are identical.

This model can be used to test algorithms for various model sizes.

## Syntax
```
context = @dynare "irbc.mod" "-DN=3";
```
will create and solve a 3 country version.

