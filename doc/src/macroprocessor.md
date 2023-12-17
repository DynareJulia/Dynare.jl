## Macro processing language

It is possible to use "macro" commands in the `.mod` file for performing
tasks such as: including modular source files, replicating blocks of
equations through loops, conditionally executing some code, writing
indexed sums or products inside equations...

The Dynare macro-language provides a new set of *macro-commands* which
can be used in `.mod` files. It features:

-   File inclusion
-   Loops (`for` structure)
-   Conditional inclusion (`if/then/else` structures)
-   Expression substitution

This macro-language is totally independent of the basic Dynare language,
and is processed by a separate component of the Dynare pre-processor.
The macro processor transforms a `.mod` file with macros into a `.mod`
file without macros (doing expansions/inclusions), and then feeds it to
the Dynare parser. The key point to understand is that the macro
processor only does text substitution (like the C preprocessor or the
PHP language). Note that it is possible to see the output of the macro
processor by using the `savemacro` option of the `dynare` command (see
`dyn-invoc`).

The macro processor is invoked by placing *macro directives* in the
`.mod` file. Directives begin with an at-sign followed by a pound sign
(`@#`). They produce no output, but give instructions to the macro
processor. In most cases, directives occupy exactly one line of text. If
needed, two backslashes (`\\`) at the end of the line indicate that the
directive is continued on the next line. Macro directives following `//`
are not interpreted by the macro processor. For historical reasons,
directives in commented blocks, *ie* surrounded by `/*` and `*/`, are
interpreted by the macro processor. The user should not rely on this
behavior. The main directives are:

-   `@#includepath`, paths to search for files that are to be
    included,
-   `@#include`, for file inclusion,
-   `@#define`, for defining a macro processor variable,
-   `@#if, @#ifdef, @#ifndef, @#elseif, @#else, @#endif` for
    conditional statements,
-   `@#for, @#endfor` for constructing loops.

The macro processor maintains its own list of variables (distinct from
model variables and Julia variables). These macro-variables are
assigned using the `@#define` directive and can be of the following
basic types: boolean, real, string, tuple, function, and array (of any
of the previous types).

### Macro expressions

Macro-expressions can be used in two places:

-   Inside macro directives, directly;
-   In the body of the `.mod` file, between an at-sign and curly
    braces: the macro processor will substitute the
    expression with its value

It is possible to construct macro-expressions that can be assigned to
macro-variables or used within a macro-directive. The expressions are
constructed using literals of the basic types (boolean, real, string,
tuple, array), comprehensions, macro-variables, macro-functions, and
standard operators.

!!! note

    Elsewhere in the manual, MACRO_EXPRESSION designates an expression
    constructed as explained in this section.

**Boolean**

The following operators can be used on booleans:

-   Comparison operators: `==, !=`
-   Logical operators: `&&, ||, !`

**Real**

The following operators can be used on reals:

-   Arithmetic operators: `+, -, *, /, ^`
-   Comparison operators: `<, >, <=, >=, ==, !=`
-   Logical operators: `&&, ||, !`
-   Ranges with an increment of `1`: `REAL1:REAL2` (for example, `1:4`
    is equivalent to real array `[1, 2, 3, 4]`).

    4.6 Previously, putting brackets around the arguments to the colon
    operator (e.g. `[1:4]`) had no effect. Now, `[1:4]` will create an
    array containing an array (i.e. `[ [1, 2, 3, 4] ]`).

-   Ranges with user-defined increment: `REAL1:REAL2:REAL3` (for
    example, `6:-2.1:-1` is equivalent to real array
    `[6, 3.9, 1.8, -0.3]`).
-   Functions:
    `max, min, mod, exp, log, log10, sin, cos, tan, asin, acos, atan, sqrt, cbrt, sign, floor, ceil, trunc, erf, erfc, gamma, lgamma, round, normpdf, normcdf`.
    NB `ln` can be used instead of `log`

**String**

String literals have to be enclosed by **double** quotes (like
`"name"`).

The following operators can be used on strings:

-   Comparison operators: `<, >, <=, >=, ==, !=`
-   Concatenation of two strings: `+`
-   Extraction of substrings: if `s` is a string, then `s[3]` is a
    string containing only the third character of `s`, and `s[4:6]`
    contains the characters from 4th to 6th
-   Function: `length`

**Tuple**

Tuples are enclosed by parenthesis and elements separated by commas
(like `(a,b,c)` or `(1,2,3)`).

The following operators can be used on tuples:

-   Comparison operators: `==, !=`
-   Functions: `empty, length`

**Array**

Arrays are enclosed by brackets, and their elements are separated by
commas (like `[1,[2,3],4]` or `["US", "FR"]`).

The following operators can be used on arrays:

-   Comparison operators: `==, !=`
-   Dereferencing: if `v` is an array, then `v[2]` is its 2nd element
-   Concatenation of two arrays: `+`
-   Set union of two arrays: `|`
-   Set intersection of two arrays: `&`
-   Difference `-`: returns the first operand from which the elements
    of the second operand have been removed.
-   Cartesian product of two arrays: `*`
-   Cartesian product of one array N times: `^N`
-   Extraction of sub-arrays: e.g. `v[4:6]`
-   Testing membership of an array: `in` operator (for example: `"b"`
    in `["a", "b", "c"]` returns `1`)
-   Functions: `empty, sum, length`

**Comprehension**

Comprehension syntax is a shorthand way to make arrays from other
arrays. There are three different ways the comprehension syntax can be
employed: [filtering], [mapping], and [filtering
and mapping].

**Filtering**

Filtering allows one to choose those elements from an array for which a
certain condition hold.

*Example*

Create a new array, choosing the even numbers from the array `1:5`:

```
     [ i in 1:5 when mod(i,2) == 0 ]
```

would result in:

```
   [2, 4]
```

**Mapping**

Mapping allows you to apply a transformation to every element of an
array.

*Example*

Create a new array, squaring all elements of the array `1:5`:

```
   [ i^2 for i in 1:5 ]
```

would result in:

```
    [1, 4, 9, 16, 25]
```

**Filtering and Mapping**

Combining the two preceding ideas would allow one to apply a
transformation to every selected element of an array.

*Example*

Create a new array, squaring all even elements of the array `1:5`:

```
   [ i^2 for i in 1:5 when mod(i,2) == 0]
```

would result in:

```
    [4, 16]
```

*Further Examples* :

```
   [ (j, i+1) for (i,j) in (1:2)^2 ]
   [ (j, i+1) for (i,j) in (1:2)*(1:2) when i < j ]
```

would result in:

```
     [(1, 2), (2, 2), (1, 3), (2, 3)]
     [(2, 2)]
```

**Function**

Functions can be defined in the macro processor using the `@#define`
directive (see below). A function is evaluated at the time it is
invoked, not at define time. Functions can be included in expressions
and the operators that can be combined with them depend on their return
type.

**Checking variable type**

Given a variable name or literal, you can check the type it evaluates to
using the following functions: `isboolean`, `isreal`, `isstring`,
`istuple`, and `isarray`.

*Examples*

| **Code**          |     **Output** |
|:----------------  |:-------------- |
| `isboolean(0)`    |     `false`    |
| `isboolean(true)` |     `true`     |
| `isreal("str")`   |     `false`    |


**Casting between types**

Variables and literals of one type can be cast into another type. Some
type changes are straightforward (e.g. changing a [real]{.title-ref} to
a [string]{.title-ref}) whereas others have certain requirements (e.g.
to cast an [array]{.title-ref} to a [real]{.title-ref} it must be a one
element array containing a type that can be cast to [real]{.title-ref}).

*Examples*

|**Code**                      |     **Output** |
|:-----------------------------|:---------------|
|`(bool) -1.1`                 |     `true`	|
|`(bool) 0`                    |     `false`	|
|`(real) "2.2"`                |     `2.2`	|
|`(tuple) [3.3]`               |     `(3.3)`	|
|`(array) 4.4`                 |     `[4.4]`	|
|`(real) [5.5]`                |     `5.5`	|
|`(real) [6.6, 7.7]`           |     `error`	|
|`(real) "8.8 in a string"`    |     `error`	|


Casts can be used in expressions:

*Examples*


|**Code**                 |   **Output** |
|:------------------------|:-------------|
|`(bool) 0 && true`       |   `false`	 |
|`(real) "1" + 2`         |   `3`	 |
|`(string) (3 + 4)`       |   `"7"`	 |
|`(array) 5 + (array) 6`  |   `[5, 6]`	 |


### Macro directives

*Macro Directive*: `@#includepath "PATH"`

*Macro Directive* `@#includepath MACRO_EXPRESSION`

This directive adds the path contained in PATH to the list of those to
search when looking for a `.mod` file specified by `@#include`. If
provided with a MACRO_EXPRESSION argument, the argument must evaluate
to a string. Note that these paths are added *after* any paths passed
using `-I <-I\<\<path\>\>>`{.interpreted-text role="opt"}.

*Example*

```
     @#includepath "/path/to/folder/containing/modfiles"
     @#includepath folders_containing_mod_files
```

*Macro Directive*: `@#include "FILENAME"` 

*Macro Directive*: `@#include MACRO_EXPRESSION`

This directive simply includes the content of another file in its place;
it is exactly equivalent to a copy/paste of the content of the included
file. If provided with a MACRO_EXPRESSION argument, the argument must
evaluate to a string. Note that it is possible to nest includes (i.e. to
include a file from an included file). The file will be searched for in
the current directory. If it is not found, the file will be searched for
in the folders provided by `-I <-I\<\<path\>\>>`{.interpreted-text
role="opt"} and `@#includepath`.

*Example*

```
     @#include "modelcomponent.mod"
     @#include location_of_modfile
```

*Macro Directive*: `@#define MACRO_VARIABLE` 

*Macro Directive*: `@#define MACRO_VARIABLE = MACRO_EXPRESSION`

*Macro Directive*: `@#define MACRO_FUNCTION = MACRO_EXPRESSION`

Defines a macro-variable or macro function.

*Example*

```
     @#define var                      // Equals true
     @#define x = 5                    // Real
     @#define y = "US"                 // String
     @#define v = [ 1, 2, 4 ]          // Real array
     @#define w = [ "US", "EA" ]       // String array
     @#define u = [ 1, ["EA"] ]        // Mixed array
     @#define z = 3 + v[2]             // Equals 5
     @#define t = ("US" in w)          // Equals true
     @#define f(x) = " " + x + y       // Function `f` with argument `x`
                                       // returns the string ' ' + x + 'US'
```

*Example*

```
     @#define x = 1
     @#define y = [ "B", "C" ]
     @#define i = 2
     @#define f(x) = x + " + " + y[i]
     @#define i = 1

     model;
       A = @{y[i] + f("D")};
     end;
```

The latter is strictly equivalent to:

```
     model;
       A = BD + B;
     end;
```

*Macro Directive*: `@#if MACRO_EXPRESSION`

*Macro Directive*: `@#ifdef MACRO_VARIABLE`

*Macro Directive*: `@#ifndef MACRO_VARIABLE` 

*Macro Directive*: `@#elseif MACRO_EXPRESSION` 

*Macro Directive*: `@#else @#endif`

Conditional inclusion of some part of the `.mod` file. The lines between
`@#if`, `@#ifdef`, or `@#ifndef` and the next `@#elseif`, `@#else` or
`@#endif` is executed only if the condition evaluates to `true`.
Following the `@#if` body, you can zero or more `@#elseif` branches. An
`@#elseif` condition is only evaluated if the preceding `@#if` or
`@#elseif` condition evaluated to `false`. The `@#else` branch is
optional and is only evaluated if all `@#if` and `@#elseif` statements
evaluate to false.

Note that when using `@#ifdef`, the condition will evaluate to `true` if
the MACRO_VARIABLE has been previously defined, regardless of its
value. Conversely, `@#ifndef` will evaluate to true if the
MACRO_VARIABLE has not yet been defined.

Note that when using `@#elseif` you can check whether or not a variable
has been defined by using the `defined` operator. Hence, to enter the
body of an `@#elseif` branch if the variable `X` has not been defined,
you would write: `@#elseif !defined(X)`.

Note that if a real appears as the result of the MACRO_EXPRESSION, it
will be interpreted as a boolean; a value of `0` is interpreted as
`false`, otherwise it is interpreted as `true`. Further note that
because of the imprecision of reals, extra care must be taken when
testing them in the MACRO_EXPRESSION. For example, `exp(log(5)) == 5`
will evaluate to `false`. Hence, when comparing real values, you should
generally use a zero tolerance around the value desired, e.g.
`exp(log(5)) > 5-1e-14 && exp(log(5)) < 5+1e-14`

*Example*

Choose between two alternative monetary policy rules using a
macro-variable:

```
     @#define linear_mon_pol = false // 0 would be treated the same
     ...
     model;
     @#if linear_mon_pol
       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
     @#else
       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
     @#endif
     ...
     end;
```

This would result in:

```
     ...
     model;
       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
     ...
     end;
```

*Example*

Choose between two alternative monetary policy rules using a
macro-variable. The only difference between this example and the
previous one is the use of `@#ifdef` instead of `@#if`. Even though
`linear_mon_pol` contains the value `false` because `@#ifdef` only
checks that the variable has been defined, the linear monetary policy
is output:

```
     @#define linear_mon_pol = false // 0 would be treated the same
     ...
     model;
     @#ifdef linear_mon_pol
       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
     @#else
       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
     @#endif
     ...
     end;
```

This would result in:

```
     ...
     model;
       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
     ...
     end;
```

*Macro Directive*: `@#for MACRO_VARIABLE in MACRO_EXPRESSION` 

*Macro Directive*: `@#for MACRO_VARIABLE in MACRO_EXPRESSION when MACRO_EXPRESSION`

*Macro Directive*: `@#for MACRO_TUPLE in MACRO_EXPRESSION` 

*Macro Directive*: `@#for MACRO_TUPLE in MACRO_EXPRESSION when MACRO\_EXPRESSION` 

*Macro Directive*: `@#endfor`

Loop construction for replicating portions of the `.mod` file. Note that
this construct can enclose variable/parameters declaration,
computational tasks, but not a model declaration.

*Example*

```
     model;
     @#for country in [ "home", "foreign" ]
       GDP_@{country} = A * K_@{country}^a * L_@{country}^(1-a);
     @#endfor
     end;
```

The latter is equivalent to:

```
     model;
       GDP_home = A * K_home^a * L_home^(1-a);
       GDP_foreign = A * K_foreign^a * L_foreign^(1-a);
     end;
```

*Example*

```
     model;
     @#for (i, j) in ["GDP"] * ["home", "foreign"]
       @{i}_@{j} = A * K_@{j}^a * L_@{j}^(1-a);
     @#endfor
     end;
```

The latter is equivalent to:

```
     model;
       GDP_home = A * K_home^a * L_home^(1-a);
       GDP_foreign = A * K_foreign^a * L_foreign^(1-a);
     end;
```

*Example*

```
     @#define countries = ["US", "FR", "JA"]
     @#define nth_co = "US"
     model;
     @#for co in countries when co != nth_co
       (1+i_@{co}) = (1+i_@{nth_co}) * E_@{co}(+1) / E_@{co};
     @#endfor
       E_@{nth_co} = 1;
     end;
```

The latter is equivalent to:

```
     model;
       (1+i_FR) = (1+i_US) * E_FR(+1) / E_FR;
       (1+i_JA) = (1+i_US) * E_JA(+1) / E_JA;
       E_US = 1;
     end;
```

*Macro Directive*: `@#echo MACRO_EXPRESSION`

Asks the preprocessor to display some message on standard output. The
argument must evaluate to a string.

*Macro Directive*: `@#error MACRO_EXPRESSION`

Asks the preprocessor to display some error message on standard output
and to abort. The argument must evaluate to a string.

*Macro Directive*: `@#echomacrovars`

*Macro Directive*: `@#echomacrovars MACRO_VARIABLE_LIST`

*Macro Directive*: `@#echomacrovars(save) MACRO_VARIABLE_LIST`

Asks the preprocessor to display the value of all macro variables up
until this point. If the `save` option is passed, then values of the
macro variables are saved to `options_.macrovars_line_<<line_numbers>>`.
If `NAME_LIST` is passed, only display/save variables and functions with
that name.

*Example*

```
     @#define A = 1
     @#define B = 2
     @#define C(x) = x*2
     @#echomacrovars A C D
```

The output of the command above is:

```
     Macro Variables:
       A = 1
     Macro Functions:
       C(x) = (x * 2)
```

### Typical usages

#### Modularization

The `@#include` directive can be used to split `.mod` files into several
modular components.

Example setup:

`modeldesc.mod`

Contains variable declarations, model equations, and shocks
declarations.

`simul.mod`

Includes `modeldesc.mod`, calibrates parameter,s and runs stochastic
simulations.

`estim.mod`

Includes `modeldesc.mod`, declares priors on parameters, and runs
Bayesian estimation.

Dynare can be called on `simul.mod` and `estim.mod` but it makes no
sense to run it on `modeldesc.mod`.

The main advantage is that you don't have to copy/paste the whole model
(at the beginning) or changes to the model (during development).

#### Indexed sums of products

The following example shows how to construct a moving average:

```
    @#define window = 2

    var x MA_x;
    ...
    model;
    ...
    MA_x = @{1/(2*window+1)}*(
    @#for i in -window:window
            +x(@{i})
    @#endfor
           );
    ...
    end;
```

After macro processing, this is equivalent to:

```
    var x MA_x;
    ...
    model;
    ...
    MA_x = 0.2*(
            +x(-2)
            +x(-1)
            +x(0)
            +x(1)
            +x(2)
           );
    ...
    end;
```

#### Multi-country models

Here is a skeleton example for a multi-country model:

```
    @#define countries = [ "US", "EA", "AS", "JP", "RC" ]
    @#define nth_co = "US"

    @#for co in countries
    var Y_@{co} K_@{co} L_@{co} i_@{co} E_@{co} ...;
    parameters a_@{co} ...;
    varexo ...;
    @#endfor

    model;
    @#for co in countries
     Y_@{co} = K_@{co}^a_@{co} * L_@{co}^(1-a_@{co});
    ...
    @#if co != nth_co
     (1+i_@{co}) = (1+i_@{nth_co}) * E_@{co}(+1) / E_@{co}; // UIP relation
    @#else
     E_@{co} = 1;
    @#endif
    @#endfor
    end;
```

#### Endogeneizing parameters

When calibrating the model, it may be useful to consider a parameter as
an endogenous variable (and vice-versa).

For example, suppose production is defined by a CES function:

$$y = \left(\alpha^{1/\xi} \ell^{1-1/\xi}+(1-\alpha)^{1/\xi}k^{1-1/\xi}\right)^{\xi/(\xi-1)}$$

and the labor share in GDP is defined as:

$$\textrm{lab\_rat} = (w \ell)/(p y)$$

In the model, $\alpha$ is a (share) parameter and `lab_rat` is an
endogenous variable.

It is clear that calibrating $\alpha$ is not straightforward; on the
contrary, we have real world data for `lab_rat` and it is clear that
these two variables are economically linked.

The solution is to use a method called *variable flipping*, which
consists in changing the way of computing the steady state. During this
computation, $\alpha$ will be made an endogenous variable and `lab_rat`
will be made a parameter. An economically relevant value will be
calibrated for `lab_rat`, and the solution algorithm will deduce the
implied value for $\alpha$.

An implementation could consist of the following files:

`modeqs.mod`

This file contains variable declarations and model equations. The code
for the declaration of $\alpha$ and `lab_rat` would look like:

```
     @#if steady
       var alpha;
       parameter lab_rat;
     @#else
       parameter alpha;
       var lab_rat;
     @#endif
```

`steady.mod`

This file computes the steady state. It begins with:

```
     @#define steady = 1
     @#include "modeqs.mod"
```

Then it initializes parameters (including `lab_rat`, excluding
$\alpha$), computes the steady state (using guess values for
endogenous, including $\alpha$), then saves values of parameters and
endogenous at steady state in a file, using the
`save_params_and_steady_state` command.

`simul.mod`

This file computes the simulation. It begins with:

```
     @#define steady = 0
     @#include "modeqs.mod"
```

Then it loads values of parameters and endogenous at steady state from
file, using the `load_params_and_steady_state` command, and computes
the simulations.

### Julia loops versus macro processor loops

Suppose you have a model with a parameter $\rho$ and you want to run
simulations for three values: $\rho = 0.8, 0.9,
1$. There are several ways of doing this:

*With a Julia loop*

```
     rhos = [ 0.8, 0.9, 1];
     for i = 1:length(rhos)
       rho = rhos(i);
       stoch_simul(order=1);
     end
```

Here the loop is not unrolled, Julia manages the iterations.
This is interesting when there are a lot of iterations.

*With a macro processor loop (case 1)*

```
     rhos = [ 0.8, 0.9, 1];
     @#for i in 1:3
       rho = rhos(@{i});
       stoch_simul(order=1);
     @#endfor
```

This is very similar to the previous example, except that the loop is
unrolled. The macro processor manages the loop index but not the data
array (`rhos`).

*With a macro processor loop (case 2)*

```
     @#for rho_val in [ 0.8, 0.9, 1]
       rho = @{rho_val};
       stoch_simul(order=1);
     @#endfor
```

The advantage of this method is that it uses a shorter syntax, since
the list of values is directly given in the loop construct. The
inconvenience is that you can not reuse the macro array in
Julia.

## Verbatim inclusion

Pass everything contained within the verbatim block to the
`<mod_file>.m` file.

*Block*: `verbatim ;`

By default, whenever Dynare encounters code that is not understood by
the parser, it is directly passed to the preprocessor output.

In order to force this behavior you can use the `verbatim` block. This
is useful when the code you want passed to the `<mod_file>.m` file
contains tokens recognized by the Dynare preprocessor.

*Example*

```
     verbatim;
     % Anything contained in this block will be passed
     % directly to the <modfile>.m file, including comments
     var = 1;
     end;
```

## Misc commands

*Command*: `save_params_and_steady_state (FILENAME);    

For all parameters, endogenous and exogenous variables, stores their
value in a text file, using a simple name/value associative table.

-   for parameters, the value is taken from the last parameter
    initialization.
-   for exogenous, the value is taken from the last `initval` block.
-   for endogenous, the value is taken from the last steady state
    computation (or, if no steady state has been computed, from the
    last `initval` block).

Note that no variable type is stored in the file, so that the values can
be reloaded with `load_params_and_steady_state` in a setup where the
variable types are different.

The typical usage of this function is to compute the steady-state of a
model by calibrating the steady-state value of some endogenous variables
(which implies that some parameters must be endogeneized during the
steady-state computation).

You would then write a first `.mod` file which computes the steady state
and saves the result of the computation at the end of the file, using
`save_params_and_steady_state`.

In a second file designed to perform the actual simulations, you would
use `load_params_and_steady_state` just after your variable
declarations, in order to load the steady state previously computed
(including the parameters which had been endogeneized during the steady
state computation).

The need for two separate `.mod` files arises from the fact that the
variable declarations differ between the files for steady state
calibration and for simulation (the set of endogenous and parameters
differ between the two); this leads to different `var` and `parameters`
statements.

Also note that you can take advantage of the `@#include` directive to
share the model equations between the two files (see
`macro-proc-lang`).

- `load_params_and_steady_state (FILENAME);`

For all parameters, endogenous and exogenous variables, loads their
value from a file created with `save_params_and_steady_state`.

-   for parameters, their value will be initialized as if they had
    been calibrated in the `.mod` file.
-   for endogenous and exogenous variables, their value will be
    initialized as they would have been from an `initval` block .

This function is used in conjunction with
`save_params_and_steady_state`; see the documentation of that function
for more information.
