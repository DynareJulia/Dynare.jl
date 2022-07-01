::: {.default-domain}
dynare
:::

Time Series
===========

Dynare provides a MATLAB/Octave class for handling time series data,
which is based on a class for handling dates. Dynare also provides a new
type for dates, so that the user does not have to worry about class and
methods for dates. Below, you will first find the class and methods used
for creating and dealing with dates and then the class used for using
time series. Dynare also provides an interface to the X-13 ARIMA-SEATS
seasonal adjustment program produced, distributed, and maintained by the
US Census Bureau (2017).

Dates
-----

::: {.highlight}
matlab
:::

### Dates in a mod file

Dynare understands dates in a mod file. Users can declare annual,
bi-annual, quarterly, or monthly dates using the following syntax:

    1990Y
    1990S2
    1990Q4
    1990M11

Behind the scene, Dynare's preprocessor translates these expressions
into instantiations of the MATLAB/Octave's class `dates` described
below. Basic operations can be performed on dates:

**plus binary operator (+)**

> An integer scalar, interpreted as a number of periods, can be added to
> a date. For instance, if `a = 1950Q1` then `b = 1951Q2` and
> `b = a + 5` are identical.

**plus unary operator (+)**

> Increments a date by one period. `+1950Q1` is identical to `1950Q2`,
> `++++1950Q1` is identical to `1951Q1`.

**minus binary operator (-)**

> Has two functions: difference and subtraction. If the second argument
> is a date, calculates the difference between the first date and the
> secmond date (e.g. `1951Q2-1950Q1` is equal to `5`). If the second
> argument is an integer `X`, subtracts `X` periods from the date (e.g.
> `1951Q2-2` is equal to `1950Q4`).

**minus unary operator (-)**

> Subtracts one period to a date. `-1950Q1` is identical to `1949Q4`.
> The unary minus operator is the reciprocal of the unary plus operator,
> `+-1950Q1` is identical to `1950Q1`.

**colon operator (:)**

> Can be used to create a range of dates. For instance,
> `r = 1950Q1:1951Q1` creates a `dates` object with five elements:
> `1950Q1, 1950Q2, 1950Q3, 1950Q4` and `1951Q1`. By default the
> increment between each element is one period. This default can be
> changed using, for instance, the following instruction:
> `1950Q1:2:1951Q1` which will instantiate a `dates` object with three
> elements: `1950Q1`, `1950Q3` and `1951Q1`.

**horzcat operator (\[,\])**

> Concatenates dates objects without removing repetitions. For instance
> `[1950Q1, 1950Q2]` is a `dates` object with two elements (`1950Q1` and
> `1950Q2`).

**vertcat operator (\[;\])**

> Same as `horzcat` operator.

**eq operator (equal, ==)**

> Tests if two `dates` objects are equal. `+1950Q1==1950Q2` returns
> `true`, `1950Q1==1950Q2` returns `false`. If the compared objects have
> both `n>1` elements, the `eq` operator returns a column vector, `n` by
> `1`, of logicals.

**ne operator (not equal, \~=)**

> Tests if two `dates` objects are not equal. `+1950Q1~=` returns
> `false` while `1950Q1~=1950Q2` returns `true`. If the compared objects
> both have `n>1` elements, the `ne` operator returns an `n` by `1`
> column vector of logicals.

**lt operator (less than, \<)**

> Tests if a `dates` object preceeds another `dates` object. For
> instance, `1950Q1<1950Q3` returns `true`. If the compared objects have
> both `n>1` elements, the `lt` operator returns a column vector, `n` by
> `1`, of logicals.

**gt operator (greater than, \>)**

> Tests if a `dates` object follows another `dates` object. For
> instance, `1950Q1>1950Q3` returns `false`. If the compared objects
> have both `n>1` elements, the `gt` operator returns a column vector,
> `n` by `1`, of logicals.

**le operator (less or equal, \<=)**

> Tests if a `dates` object preceeds another `dates` object or is equal
> to this object. For instance, `1950Q1<=1950Q3` returns `true`. If the
> compared objects have both `n>1` elements, the `le` operator returns a
> column vector, `n` by `1`, of logicals.

**ge operator (greater or equal, \>=)**

> Tests if a `dates` object follows another `dates` object or is equal
> to this object. For instance, `1950Q1>=1950Q3` returns `false`. If the
> compared objects have both `n>1` elements, the `ge` operator returns a
> column vector, `n` by `1`, of logicals.

One can select an element, or some elements, in a `dates` object as he
would extract some elements from a vector in MATLAB/Octave. Let
`a = 1950Q1:1951Q1` be a `dates` object, then `a(1)==1950Q1` returns
`true`, `a(end)==1951Q1` returns `true` and `a(end-1:end)` selects the
two last elements of `a` (by instantiating the `dates` object
`[1950Q4, 1951Q1]`).

Remark: Dynare substitutes any occurrence of dates in the `.mod` file
into an instantiation of the `dates` class regardless of the context.
For instance, `d = 1950Q1` will be translated as `d = dates('1950Q1');`.
This automatic substitution can lead to a crash if a date is defined in
a string. Typically, if the user wants to display a date:

    disp('Initial period is 1950Q1');

Dynare will translate this as:

    disp('Initial period is dates('1950Q1')');

which will lead to a crash because this expression is illegal in MATLAB.
For this situation, Dynare provides the `$` escape parameter. The
following expression:

    disp('Initial period is $1950Q1');

will be translated as:

    disp('Initial period is 1950Q1');

in the generated MATLAB script.

### The dates class {#dates-members}

::: {.dates}

arg freq

:   equal to 1, 2, 4, 12 or 365 (resp. for annual, bi-annual, quarterly,
    monthly, or daily dates).

arg time

:   a `n*1` array of integers, the number of periods since year 0 ().

Each member is private, one can display the content of a member but
cannot change its value directly. Note also that it is not possible to
mix frequencies in a `dates` object: all the elements must have common
frequency.

The `dates` class has the following constructors:

::: {.construct}
dates() dates(FREQ)

Returns an empty `dates` object with a given frequency (if the
constructor is called with one input argument). `FREQ` is a character
equal to 'Y' or 'A' for annual dates, 'S' or 'H' for bi-annual dates,
'Q' for quarterly dates, 'M' for monthly dates, or 'D' for daily dates.
Note that `FREQ` is not case sensitive, so that, for instance, 'q' is
also allowed for quarterly dates. The frequency can also be set with an
integer scalar equal to 1 (annual), 2 (bi-annual), 4 (quarterly), 12
(monthly), or 365 (daily). The instantiation of empty objects can be
used to rename the `dates` class. For instance, if one only works with
quarterly dates, object `qq` can be created as:

    qq = dates('Q')

and a `dates` object holding the date `2009Q2`:

    d0 = qq(2009,2);

which is much simpler if `dates` objects have to be defined
programmatically. For daily dates, we would instantiate an empty daily
dates object as:

    dd = dates('D')

and a `dates` object holding the date `2020-12-31`:

    d1 = dd(2020,12,31);
:::

::: {.construct}
dates(STRING) dates(STRING, STRING, \...)

Returns a `dates` object that represents a date as given by the string
`STRING`. This string has to be interpretable as a date (only strings of
the following forms are admitted: `'1990Y'`, `'1990A'`, `1990S1`,
`1990H1`, `'1990Q1'`, `'1990M2'`, or `'2020-12-31'`), the routine
`isdate` can be used to test if a string is interpretable as a date. If
more than one argument is provided, they should all be dates represented
as strings, the resulting `dates` object contains as many elements as
arguments to the constructor. For the daily dates, the string must be of
the form yyyy-mm-dd with two digits for the months (mm) and days (dd),
even if the number of days or months is smaller than ten (in this case a
leading 0 is required).
:::

::: {.construct}
dates(DATES) dates(DATES, DATES, \...)

Returns a copy of the `dates` object `DATES` passed as input arguments.
If more than one argument is provided, they should all be `dates`
objects. The number of elements in the instantiated `dates` object is
equal to the sum of the elements in the `dates` passed as arguments to
the constructor.
:::

::: {.construct}
dates (FREQ, YEAR, SUBPERIOD\[, S\])

where `FREQ` is a single character ('Y', 'A', 'S', 'H', 'Q', 'M', 'D')
or integer (1, 2, 4, 12, or 365) specifying the frequency, `YEAR` and
`SUBPERIOD` and `S` are `n*1` vectors of integers. Returns a `dates`
object with `n` elements. The last argument, `S`, is only to be used for
daily frequency. If `FREQ` is equal to `'Y'`, `'A'` or `1`, the third
argument is not needed (because `SUBPERIOD` is necessarily a vector of
ones in this case).
:::

*Example*

>     do1 = dates('1950Q1');
>     do2 = dates('1950Q2','1950Q3');
>     do3 = dates(do1,do2);
>     do4 = dates('Q',1950, 1);
>     do5 = dates('D',1973, 1, 25);

A list of the available methods, by alphabetical order, is given below.
Note that by default the methods do not allow in place modifications:
when a method is applied to an object a new object is instantiated. For
instance, to apply the method `multiplybytwo` to an object `X` we write:

    >> X = 2;
    >> Y = X.multiplybytwo();
    >> X

    2

    >> Y

    4

or equivalently:

    >> Y = multiplybytwo(X);

the object `X` is left unchanged, and the object `Y` is a modified copy
of `X` (multiplied by two). This behaviour is altered if the name of the
method is postfixed with an underscore. In this case the creation of a
copy is avoided. For instance, following the previous example, we would
have:

    >> X = 2;
    >> X.multiplybytwo_();
    >> X

    4

Modifying the objects in place, with underscore methods, is particularly
useful if the methods are called in loops, since this saves the object
instantiation overhead.

::: {.datesmethod}
C = append (A, B) C = [append]() (A, B)

Appends `dates` object `B`, or a string that can be interpreted as a
date, to the `dates` object `A`. If `B` is a `dates` object it is
assumed that it has no more than one element.

*Example*

>     >> D = dates('1950Q1','1950Q2');
>     >> d = dates('1950Q3');
>     >> E = D.append(d);
>     >> F = D.append('1950Q3');
>     >> isequal(E,F)
>
>     ans =
>
>          1
>     >> F
>
>     F = <dates: 1950Q1, 1950Q2, 1950Q3>
>
>     >> D
>
>     D = <dates: 1950Q1, 1950Q2>
>
>     >> D.append_('1950Q3')
>
>     ans = <dates: 1950Q1, 1950Q2, 1950Q3>
:::

::: {.datesmethod}
B = char (A)

Overloads the MATLAB/Octave `char` function. Converts a `dates` object
into a character array.

*Example*

>     >> A = dates('1950Q1');
>     > A.char()
>
>     ans =
>
>     '1950Q1'
:::

::: {.datesmethod}
C = colon (A, B) C = colon (A, i, B)

Overloads the MATLAB/Octave colon (`:`) operator. A and B are `dates`
objects. The optional increment `i` is a scalar integer (default value
is `i=1`). This method returns a `dates` object and can be used to
create ranges of dates.

*Example*

>     >> A = dates('1950Q1');
>     >> B = dates('1951Q2');
>     >> C = A:B
>
>     C = <dates: 1950Q1, 1950Q2, 1950Q3, 1950Q4, 1951Q1>
>
>     >> D = A:2:B
>
>     D = <dates: 1950Q1, 1950Q3, 1951Q1>
:::

::: {.datesmethod}
B = copy (A)

Returns a copy of a `dates` object.
:::

::: {.datesmethod}
disp (A)

Overloads the MATLAB/Octave disp function for `dates` object.
:::

::: {.datesmethod}
display (A)

Overloads the MATLAB/Octave display function for `dates` object.

*Example*

>     >> disp(B)
>
>     B = <dates: 1950Q1, 1950Q2, 1950Q3, 1950Q4, 1951Q1, 1951Q2, 1951Q3, 1951Q4, 1952Q1, 1952Q2, 1952Q3>
>
>
>     >> display(B)
>
>     B = <dates: 1950Q1, 1950Q2, ..., 1952Q2, 1952Q3>
:::

::: {.datesmethod}
B = double (A)

Overloads the MATLAB/Octave `double` function. `A` is a `dates` object.
The method returns a floating point representation of a `dates` object,
the integer and fractional parts respectively corresponding to the year
and the subperiod. The fractional part is the subperiod number minus one
divided by the frequency (`1`, `4`, or `12`).

*Example*:

>     >> a = dates('1950Q1'):dates('1950Q4');
>     >> a.double()
>
>     ans =
>
>          1950.00
>          1950.25
>          1950.50
>          1950.75
:::

::: {.datesmethod}
C = eq (A, B)

Overloads the MATLAB/Octave `eq` (equal, `==`) operator. `dates` objects
`A` and `B` must have the same number of elements (say, `n`). The
returned argument is a `n` by `1` vector of logicals. The i-th element
of `C` is equal to `true` if and only if the dates `A(i)` and `B(i)` are
the same.

*Example*

>     >> A = dates('1950Q1','1951Q2');
>     >> B = dates('1950Q1','1950Q2');
>     >> A==B
>
>     ans =
>
>       2x1 logical array
>
>        1
>        0
:::

::: {.datesmethod}
C = ge (A, B)

Overloads the MATLAB/Octave `ge` (greater or equal, `>=`) operator.
`dates` objects `A` and `B` must have the same number of elements (say,
`n`). The returned argument is a `n` by `1` vector of logicals. The i-th
element of `C` is equal to `true` if and only if the date `A(i)` is
posterior or equal to the date `B(i)`.

*Example*

>     >> A = dates('1950Q1','1951Q2');
>     >> B = dates('1950Q1','1950Q2');
>     >> A>=B
>
>     ans =
>
>       2x1 logical array
>
>        1
>        1
:::

::: {.datesmethod}
C = gt (A, B)

Overloads the MATLAB/Octave `gt` (greater than, `>`) operator. `dates`
objects `A` and `B` must have the same number of elements (say, `n`).
The returned argument is a `n` by `1` vector of logicals. The i-th
element of `C` is equal to `1` if and only if the date `A(i)` is
posterior to the date `B(i)`.

*Example*

>     >> A = dates('1950Q1','1951Q2');
>     >> B = dates('1950Q1','1950Q2');
>     >> A>B
>
>     ans =
>
>       2x1 logical array
>
>        0
>        1
:::

::: {.datesmethod}
D = horzcat (A, B, C, \...)

Overloads the MATLAB/Octave `horzcat` operator. All the input arguments
must be `dates` objects. The returned argument is a `dates` object
gathering all the dates given in the input arguments (repetitions are
not removed).

*Example*

>     >> A = dates('1950Q1');
>     >> B = dates('1950Q2');
>     >> C = [A, B];
>     >> C
>
>     C = <dates: 1950Q1, 1950Q2>
:::

::: {.datesmethod}
C = intersect (A, B)

Overloads the MATLAB/Octave `intersect` function. All the input
arguments must be `dates` objects. The returned argument is a `dates`
object gathering all the common dates given in the input arguments. If
`A` and `B` are disjoint `dates` objects, the function returns an empty
`dates` object. Returned dates in `dates` object `C` are sorted by
increasing order.

*Example*

>     >> A = dates('1950Q1'):dates('1951Q4');
>     >> B = dates('1951Q1'):dates('1951Q4');
>     >> C = intersect(A, B);
>     >> C
>
>     C = <dates: 1951Q1, 1951Q2, 1951Q3, 1951Q4>
:::

::: {.datesmethod}
B = isempty (A)

Overloads the MATLAB/Octave `isempty` function.

*Example*

>     >> A = dates('1950Q1');
>     >> A.isempty()
>
>     ans =
>
>       logical
>
>       0
>
>     >> B = dates();
>     >> B.isempty()
>
>     ans =
>
>       logical
>
>       1
:::

::: {.datesmethod}
C = isequal (A, B)

Overloads the MATLAB/Octave `isequal` function.

*Example*

>     >> A = dates('1950Q1');
>     >> B = dates('1950Q2');
>     >> isequal(A, B)
>
>     ans =
>
>       logical
>
>       0
:::

::: {.datesmethod}
C = le (A, B)

Overloads the MATLAB/Octave `le` (less or equal, `<=`) operator. `dates`
objects `A` and `B` must have the same number of elements (say, `n`).
The returned argument is a `n` by `1` vector of logicals. The i-th
element of `C` is equal to `true` if and only if the date `A(i)` is
anterior or equal to the date `B(i)`.

*Example*

>     >> A = dates('1950Q1','1951Q2');
>     >> B = dates('1950Q1','1950Q2');
>     >> A<=B
>
>     ans =
>
>       2x1 logical array
>
>        1
>        0
:::

::: {.datesmethod}
B = length (A)

Overloads the MATLAB/Octave `length` function. Returns the number of
elements in a `dates` object.

*Example*

>     >> A = dates('1950Q1'):dates(2000Q3);
>     >> A.length()
>
>     ans =
>
>        203
:::

::: {.datesmethod}
C = lt (A, B)

Overloads the MATLAB/Octave `lt` (less than, `<`) operator. `dates`
objects `A` and `B` must have the same number of elements (say, `n`).
The returned argument is a `n` by `1` vector of logicals. The i-th
element of `C` is equal to `true` if and only if the date `A(i)` is
anterior or equal to the date `B(i)`.

*Example*

>     >> A = dates('1950Q1','1951Q2');
>     >> B = dates('1950Q1','1950Q2');
>     >> A<B
>
>     ans =
>
>       2x1 logical array
>
>        0
>        0
:::

::: {.datesmethod}
D = max (A, B, C, \...)

Overloads the MATLAB/Octave `max` function. All input arguments must be
`dates` objects. The function returns a single element `dates` object
containing the greatest date.

*Example*

>     >> A = {dates('1950Q2'), dates('1953Q4','1876Q2'), dates('1794Q3')};
>     >> max(A{:})
>
>     ans = <dates: 1953Q4>
:::

::: {.datesmethod}
D = min (A, B, C, \...)

Overloads the MATLAB/Octave `min` function. All input arguments must be
`dates` objects. The function returns a single element `dates` object
containing the smallest date.

*Example*

>     >> A = {dates('1950Q2'), dates('1953Q4','1876Q2'), dates('1794Q3')};
>     >> min(A{:})
>
>     ans = <dates: 1794Q3>
:::

::: {.datesmethod}
C = minus (A, B)

Overloads the MATLAB/Octave `minus` operator (`-`). If both input
arguments are `dates` objects, then number of periods between `A` and
`B` is returned (so that `A+C=B`). If `B` is a vector of integers, the
minus operator shifts the `dates` object by `B` periods backward.

*Example*

>     >> d1 = dates('1950Q1','1950Q2','1960Q1');
>     >> d2 = dates('1950Q3','1950Q4','1960Q1');
>     >> ee = d2-d1
>
>     ee =
>
>          2
>          2
>          0
>
>     >> d1-(-ee)
>
>     ans = <dates: 1950Q3, 1950Q4, 1960Q1>
:::

::: {.datesmethod}
C = mtimes (A, B)

Overloads the MATLAB/Octave `mtimes` operator (`*`). `A` and `B` are
respectively expected to be a `dseries` object and a scalar integer.
Returns `dates` object `A` replicated `B` times.

*Example*

>     >> d = dates('1950Q1');
>     >> d*2
>
>     ans = <dates: 1950Q1, 1950Q1>
:::

::: {.datesmethod}
C = ne (A, B)

Overloads the MATLAB/Octave `ne` (not equal, `~=`) operator. `dates`
objects `A` and `B` must have the same number of elements (say, `n`) or
one of the inputs must be a single element `dates` object. The returned
argument is a `n` by `1` vector of logicals. The i-th element of `C` is
equal to `true` if and only if the dates `A(i)` and `B(i)` are
different.

*Example*

>     >> A = dates('1950Q1','1951Q2');
>     >> B = dates('1950Q1','1950Q2');
>     >> A~=B
>
>     ans =
>
>       2x1 logical array
>
>        0
>        1
:::

::: {.datesmethod}
C = plus (A, B)

Overloads the MATLAB/Octave `plus` operator (`+`). If both input
arguments are `dates` objects, then the method combines `A` and `B`
without removing repetitions. If `B` is a vector of integers, the `plus`
operator shifts the `dates` object by `B` periods forward.

*Example*

>     >> d1 = dates('1950Q1','1950Q2')+dates('1960Q1');
>     >> d2 = (dates('1950Q1','1950Q2')+2)+dates('1960Q1');
>     >> ee = d2-d1;
>
>     ee =
>
>          2
>          2
>          0
>
>     >> d1+ee
>     ans = <dates: 1950Q3, 1950Q4, 1960Q1>
:::

::: {.datesmethod}
C = pop (A) C = pop (A, B) C = [pop]() (A) C = [pop]() (A, B)

Pop method for `dates` class. If only one input is provided, the method
removes the last element of a `dates` object. If a second input argument
is provided, a scalar integer between `1` and `A.length()`, the method
removes element number `B` from `dates` object `A`.

*Example*

>     >> d = dates('1950Q1','1950Q2');
>     >> d.pop()
>
>     ans = <dates: 1950Q1>
>
>     >> d.pop_(1)
>
>     ans = <dates: 1950Q2>
:::

::: {.datesmethod}
C = remove (A, B) C = [remove]() (A, B)

Remove method for `dates` class. Both inputs have to be `dates` objects,
removes dates in `B` from `A`.

*Example*

>     >> d = dates('1950Q1','1950Q2');
>     >> d.remove(dates('1950Q2'))
>
>     ans = <dates: 1950Q1>
:::

::: {.datesmethod}
C = setdiff (A, B)

Overloads the MATLAB/Octave `setdiff` function. All the input arguments
must be `dates` objects. The returned argument is a `dates` object all
dates present in `A` but not in `B`. If `A` and `B` are disjoint `dates`
objects, the function returns `A`. Returned dates in `dates` object `C`
are sorted by increasing order.

*Example*

>     >> A = dates('1950Q1'):dates('1969Q4');
>     >> B = dates('1960Q1'):dates('1969Q4');
>     >> C = dates('1970Q1'):dates('1979Q4');
>     >> setdiff(A, B)
>
>     ans = <dates: 1950Q1, 1950Q2,  ..., 1959Q3, 1959Q4>
>
>     >> setdiff(A, C)
>
>     ans = <dates: 1950Q1, 1950Q2,  ..., 1969Q3, 1969Q4>
:::

::: {.datesmethod}
B = sort (A) B = [sort]() (A)

Sort method for `dates` objects. Returns a `dates` object with elements
sorted by increasing order.

*Example*

>     >> dd = dates('1945Q3','1938Q4','1789Q3');
>     >> dd.sort()
>
>     ans = <dates: 1789Q3, 1938Q4, 1945Q3>
:::

::: {.datesmethod}
B = strings (A)

Converts a `dates` object into a cell of char arrays.

*Example*

>     >> A = dates('1950Q1');
>     >> A = A:A+1;
>     >> A.strings()
>
>       ans =
>
>         1x2 cell array
>
>         {'1950Q1'}    {'1950Q2'}
:::

::: {.datesmethod}
B = subperiod (A)

Returns the subperiod of a date (an integer scalar between 1 and
`A.freq`). This method is not implemented for daily dates.

*Example*

>     >> A = dates('1950Q2');
>     >> A.subperiod()
>
>     ans =
>
>          2
:::

::: {.datesmethod}
B = uminus (A)

Overloads the MATLAB/Octave unary minus operator. Returns a `dates`
object with elements shifted one period backward.

*Example*

    >> dd = dates('1945Q3','1938Q4','1973Q1');
    >> -dd

    ans = <dates: 1945Q2, 1938Q3, 1972Q4>
:::

::: {.datesmethod}
D = union (A, B, C, \...)

Overloads the MATLAB/Octave `union` function. Returns a `dates` object
with elements sorted by increasing order (repetitions are removed, to
keep the repetitions use the `horzcat` or `plus` operators).

*Example*

>     >> d1 = dates('1945Q3','1973Q1','1938Q4');
>     >> d2 = dates('1973Q1','1976Q1');
>     >> union(d1,d2)
>
>     ans = <dates: 1938Q4, 1945Q3, 1973Q1, 1976Q1>
:::

::: {.datesmethod}
B = unique (A) B = [unique]() (A)

Overloads the MATLAB/Octave `unique` function. Returns a `dates` object
with repetitions removed (only the last occurence of a date is kept).

*Example*

>     >> d1 = dates('1945Q3','1973Q1','1945Q3');
>     >> d1.unique()
>
>     ans = <dates: 1973Q1, 1945Q3>
:::

::: {.datesmethod}
B = uplus (A)

Overloads the MATLAB/Octave unary plus operator. Returns a `dates`
object with elements shifted one period ahead.

*Example*

>     >> dd = dates('1945Q3','1938Q4','1973Q1');
>     >> +dd
>
>     ans = <dates: 1945Q4, 1939Q1, 1973Q2>
:::

::: {.datesmethod}
D = vertcat (A, B, C, \...)

Overloads the MATLAB/Octave `horzcat` operator. All the input arguments
must be `dates` objects. The returned argument is a `dates` object
gathering all the dates given in the input arguments (repetitions are
not removed).
:::

::: {.datesmethod}
B = year (A)

Returns the year of a date (an integer scalar between 1 and `A.freq`).

*Example*

>     >> A = dates('1950Q2');
>     >> A.subperiod()
>
>     ans =
>
>            1950
:::
:::

The dseries class {#dseries-members}
-----------------

::: {.dseries}
The MATLAB/Octave `dseries` class handles time series data. As any
MATLAB/Octave statements, this class can be used in a Dynare's mod file.
A `dseries` object has six members:

arg name

:   A `vobs*1` cell of strings or a `vobs*p` character array, the names
    of the variables.

arg tex

:   A `vobs*1` cell of strings or a `vobs*p` character array, the tex
    names of the variables.

arg dates dates

:   An object with `nobs` elements, the dates of the sample.

arg double data

:   A `nobs` by `vobs` array, the data.

arg ops

:   The history of operations on the variables.

arg tags

:   The user-defined tags on the variables.

`data`, `name`, `tex`, and `ops` are private members. The following
constructors are available:

::: {.construct}
dseries () dseries (INITIAL\_DATE)

Instantiates an empty `dseries` object with, if defined, an initial date
given by the single element `dates` object *INITIAL\_DATE.*
:::

::: {.construct}
dseries (FILENAME\[, INITIAL\_DATE\])

Instantiates and populates a `dseries` object with a data file specified
by *FILENAME*, a string passed as input. Valid file types are `.m`,
`.mat`, `.csv` and `.xls/.xlsx` (Octave only supports `.xlsx` files and
the [io](https://octave.sourceforge.io/io/) package from Octave-Forge
must be installed). The extension of the file should be explicitly
provided.

A typical `.m` file will have the following form:

    FREQ__ = 4;
    INIT__ = '1994Q3';
    NAMES__ = {'azert';'yuiop'};
    TEX__ = {'azert';'yuiop'};

    azert = randn(100,1);
    yuiop = randn(100,1);

If a `.mat` file is used instead, it should provide the same
informations, except that the data should not be given as a set of
vectors, but as a single matrix of doubles named `DATA__`. This array
should have as many columns as elements in `NAMES__` (the number of
variables). Note that the `INIT__` variable can be either a `dates`
object or a string which could be used to instantiate the same `dates`
object. If `INIT__` is not provided in the `.mat` or `.m` file, the
initial is by default set equal to `dates('1Y')`. If a second input
argument is passed to the constructor, `dates` object *INITIAL\_DATE*,
the initial date defined in *FILENAME* is reset to *INITIAL\_DATE*. This
is typically usefull if `INIT__` is not provided in the data file.

If an `.xlsx` file is used, the first row should be a header containing
the variable names. The first column may contain date information that
must correspond to a valid date format recognized by Dynare. If such
date information is specified in the first column, its header name must
be left empty.
:::

::: {.construct}
dseries
(DATA\_MATRIX\[,INITIAL\_DATE\[,LIST\_OF\_NAMES\[,TEX\_NAMES\]\]\])
dseries
(DATA\_MATRIX\[,RANGE\_OF\_DATES\[,LIST\_OF\_NAMES\[,TEX\_NAMES\]\]\])

If the data is not read from a file, it can be provided via a
$T \times N$ matrix as the first argument to `dseries` ' constructor,
with $T$ representing the number of observations on $N$ variables. The
optional second argument, *INITIAL\_DATE*, can be either a `dates`
object representing the period of the first observation or a string
which would be used to instantiate a `dates` object. Its default value
is `dates('1Y')`. The optional third argument, *LIST\_OF\_NAMES*, is a
$N \times 1$ cell of strings with one entry for each variable name. The
default name associated with column `i` of *DATA\_MATRIX* is
`Variable_i`. The final argument, *TEX\_NAMES*, is a $N
\times 1$ cell of strings composed of the LaTeX names associated with
the variables. The default LaTeX name associated with column `i` of
*DATA\_MATRIX* is `Variable\_i`. If the optional second input argument
is a range of dates, `dates` object *RANGE\_OF\_DATES*, the number of
rows in the first argument must match the number of elements
*RANGE\_OF\_DATES* or be equal to one (in which case the single
observation is replicated).
:::

::: {.construct}
dseries (TABLE)

Creates a `dseries` object given the MATLAB Table provided as the sole
argument. It is assumed that the first column of the table contains the
dates of the `dseries` and the first row contains the names. This
feature is not available under Octave or MATLAB R2013a or earlier.

*Example*

Various ways to create a `dseries` object:

    do1 = dseries(1999Q3);
    do2 = dseries('filename.csv');
    do3 = dseries([1; 2; 3], 1999Q3, {'var123'}, {'var_{123}'});

    >> do1 = dseries(dates('1999Q3'));
    >> do2 = dseries('filename.csv');
    >> do3 = dseries([1; 2; 3], dates('1999Q3'), {'var123'}, {'var_{123}'});
:::

One can easily create subsamples from a `dseries` object using the
overloaded parenthesis operator. If `ds` is a `dseries` object with $T$
observations and `d` is a `dates` object with $S<T$ elements, such that
$\min(d)$ is not smaller than the date associated to the first
observation in `ds` and $\max(d)$ is not greater than the date
associated to the last observation, then `ds(d)` instantiates a new
`dseries` object containing the subsample defined by `d`.

A list of the available methods, by alphabetical order, is given below.
As in the previous section the in place modifications versions of the
methods are postfixed with an underscore.

::: {.dseriesmethod}
A = abs (B) [abs]() (B)

Overloads the `abs()` function for `dseries` objects. Returns the
absolute value of the variables in dseries `object` `B`.

*Example*

>     >> ts0 = dseries(randn(3,2),'1973Q1',{'A1'; 'A2'},{'A_1'; 'A_2'});
>     >> ts1 = ts0.abs();
>     >> ts0
>
>     ts0 is a dseries object:
>
>            | A1       | A2
>     1973Q1 | -0.67284 | 1.4367
>     1973Q2 | -0.51222 | -0.4948
>     1973Q3 | 0.99791  | 0.22677
>
>     >> ts1
>
>     ts1 is a dseries object:
>
>            | abs(A1) | abs(A2)
>     1973Q1 | 0.67284 | 1.4367
>     1973Q2 | 0.51222 | 0.4948
>     1973Q3 | 0.99791 | 0.22677
:::

::: {.dseriesmethod}
\[A, B\] = align (A, B) [align]() (A, B)

If `dseries` objects `A` and `B` are defined on different time ranges,
this function extends `A` and/or `B` with NaNs so that they are defined
on the same time range. Note that both `dseries` objects must have the
same frequency.

*Example*

>     >> ts0 = dseries(rand(5,1),dates('2000Q1')); % 2000Q1 -> 2001Q1
>     >> ts1 = dseries(rand(3,1),dates('2000Q4')); % 2000Q4 -> 2001Q2
>     >> [ts0, ts1] = align(ts0, ts1);             % 2000Q1 -> 2001Q2
>     >> ts0
>
>     ts0 is a dseries object:
>
>            | Variable_1
>     2000Q1 | 0.81472
>     2000Q2 | 0.90579
>     2000Q3 | 0.12699
>     2000Q4 | 0.91338
>     2001Q1 | 0.63236
>     2001Q2 | NaN
>
>     >> ts1
>
>     ts1 is a dseries object:
>
>            | Variable_1
>     2000Q1 | NaN
>     2000Q2 | NaN
>     2000Q3 | NaN
>     2000Q4 | 0.66653
>     2001Q1 | 0.17813
>     2001Q2 | 0.12801
>
>     >> ts0 = dseries(rand(5,1),dates('2000Q1')); % 2000Q1 -> 2001Q1
>     >> ts1 = dseries(rand(3,1),dates('2000Q4')); % 2000Q4 -> 2001Q2
>     >> align_(ts0, ts1);                         % 2000Q1 -> 2001Q2
>     >> ts1
>
>     ts1 is a dseries object:
>
>            | Variable_1
>     2000Q1 | NaN
>     2000Q2 | NaN
>     2000Q3 | NaN
>     2000Q4 | 0.66653
>     2001Q1 | 0.17813
>     2001Q2 | 0.12801
:::

::: {.dseriesmethod}
C = backcast (A, B\[, diff\]) [backcast]() (A, B\[, diff\])

Backcasts `dseries` object `A` with `dseries` object B\'s growth rates
(except if the last optional argument, `diff`, is true in which case
first differences are used). Both `dseries` objects must have the same
frequency.
:::

::: {.dseriesmethod}
B = baxter\_king\_filter (A, hf, lf, K) [baxter\_king\_filter]() (A, hf,
lf, K)

Implementation of the *Baxter and King* (1999) band pass filter for
`dseries` objects. This filter isolates business cycle fluctuations with
a period of length ranging between `hf` (high frequency) to `lf` (low
frequency) using a symmetric moving average smoother with $2K+1$ points,
so that $K$ observations at the beginning and at the end of the sample
are lost in the computation of the filter. The default value for `hf` is
`6`, for `lf` is `32`, and for `K` is `12`.

*Example*

>     % Simulate a component model (stochastic trend, deterministic
>     % trend, and a stationary autoregressive process).
>     e = 0.2*randn(200,1);
>     u = randn(200,1);
>     stochastic_trend = cumsum(e);
>     deterministic_trend = .1*transpose(1:200);
>     x = zeros(200,1);
>     for i=2:200
>         x(i) = .75*x(i-1) + u(i);
>     end
>     y = x + stochastic_trend + deterministic_trend;
>
>     % Instantiates time series objects.
>     ts0 = dseries(y,'1950Q1');
>     ts1 = dseries(x,'1950Q1'); % stationary component.
>
>     % Apply the Baxter-King filter.
>     ts2 = ts0.baxter_king_filter();
>
>     % Plot the filtered time series.
>     plot(ts1(ts2.dates).data,'-k'); % Plot of the stationary component.
>     hold on
>     plot(ts2.data,'--r');           % Plot of the filtered y.
>     hold off
>     axis tight
>     id = get(gca,'XTick');
>     set(gca,'XTickLabel',strings(ts1.dates(id)));
:::

::: {.dseriesmethod}
B = center (A\[, geometric\]) [center]() (A\[, geometric\])

Centers variables in `dseries` object `A` around their arithmetic means,
except if the optional argument `geometric` is set equal to `true` in
which case all the variables are divided by their geometric means.
:::

::: {.dseriesmethod}
C = chain (A, B) [chain]() (A, B)

Merge two `dseries` objects along the time dimension. The two objects
must have the same number of observed variables, and the initial date in
`B` must not be posterior to the last date in `A`. The returned
`dseries` object, `C`, is built by extending `A` with the cumulated
growth factors of `B`.

*Example*

>     >> ts = dseries([1; 2; 3; 4],dates(`1950Q1'))
>
>     ts is a dseries object:
>
>            | Variable_1
>     1950Q1 | 1
>     1950Q2 | 2
>     1950Q3 | 3
>     1950Q4 | 4
>
>     >> us = dseries([3; 4; 5; 6],dates(`1950Q3'))
>
>     us is a dseries object:
>
>            | Variable_1
>     1950Q3 | 3
>     1950Q4 | 4
>     1951Q1 | 5
>     1951Q2 | 6
>
>     >> chain(ts, us)
>
>     ans is a dseries object:
>
>            | Variable_1
>     1950Q1 | 1
>     1950Q2 | 2
>     1950Q3 | 3
>     1950Q4 | 4
>     1951Q1 | 5
>     1951Q2 | 6
:::

::: {.dseriesmethod}
\[error\_flag, message \] = check (A)

Sanity check of `dseries` object `A`. Returns `1` if there is an error,
`0` otherwise. The second output argument is a string giving brief
informations about the error.
:::

::: {.dseriesmethod}
B = copy (A)

Returns a copy of `A`. If an inplace modification method is applied to
`A`, object `B` will not be affected. Note that if `A` is assigned to
`C`, `C = A`, then any in place modification method applied to `A` will
change `C`.

*Example*

>     >> a = dseries(randn(5,1))
>
>     a is a dseries object:
>
>        | Variable_1
>     1Y | -0.16936
>     2Y | -1.1451
>     3Y | -0.034331
>     4Y | -0.089042
>     5Y | -0.66997
>
>     >> b = copy(a);
>     >> c = a;
>     >> a.abs();
>     >> a.abs_();
>     >> a
>
>     a is a dseries object:
>
>        | Variable_1
>     1Y | 0.16936
>     2Y | 1.1451
>     3Y | 0.034331
>     4Y | 0.089042
>     5Y | 0.66997
>
>     >> b
>
>     b is a dseries object:
>
>        | Variable_1
>     1Y | -0.16936
>     2Y | -1.1451
>     3Y | -0.034331
>     4Y | -0.089042
>     5Y | -0.66997
>
>     >> c
>
>     c is a dseries object:
>
>        | Variable_1
>     1Y | 0.16936
>     2Y | 1.1451
>     3Y | 0.034331
>     4Y | 0.089042
>     5Y | 0.66997
:::

::: {.dseriesmethod}
B = cumprod (A\[, d\[, v\]\]) [cumprod]() (A\[, d\[, v\]\])

Overloads the MATLAB/Octave `cumprod` function for `dseries` objects.
The cumulated product cannot be computed if the variables in `dseries`
object `A` have NaNs. If a `dates` object `d` is provided as a second
argument, then the method computes the cumulated product with the
additional constraint that the variables in the `dseries` object `B` are
equal to one in period `d`. If a single-observation `dseries` object `v`
is provided as a third argument, the cumulated product in `B` is
normalized such that `B(d)` matches `v` (`dseries` objects `A` and `v`
must have the same number of variables).

*Example*

>     >> ts1 = dseries(2*ones(7,1));
>     >> ts2 = ts1.cumprod();
>     >> ts2
>
>     ts2 is a dseries object:
>
>        | cumprod(Variable_1)
>     1Y | 2
>     2Y | 4
>     3Y | 8
>     4Y | 16
>     5Y | 32
>     6Y | 64
>     7Y | 128
>
>     >> ts3 = ts1.cumprod(dates('3Y'));
>     >> ts3
>
>     ts3 is a dseries object:
>
>        | cumprod(Variable_1)
>     1Y | 0.25
>     2Y | 0.5
>     3Y | 1
>     4Y | 2
>     5Y | 4
>     6Y | 8
>     7Y | 16
>
>     >> ts4 = ts1.cumprod(dates('3Y'),dseries(pi));
>     >> ts4
>
>     ts4 is a dseries object:
>
>        | cumprod(Variable_1)
>     1Y | 0.7854
>     2Y | 1.5708
>     3Y | 3.1416
>     4Y | 6.2832
>     5Y | 12.5664
>     6Y | 25.1327
>     7Y | 50.2655
:::

::: {.dseriesmethod}
B = cumsum (A\[, d\[, v\]\]) cumsum (A\[, d\[, v\]\])

Overloads the MATLAB/Octave `cumsum` function for `dseries` objects. The
cumulated sum cannot be computed if the variables in `dseries` object
`A` have NaNs. If a `dates` object `d` is provided as a second argument,
then the method computes the cumulated sum with the additional
constraint that the variables in the `dseries` object `B` are zero in
period `d`. If a single observation `dseries` object `v` is provided as
a third argument, the cumulated sum in `B` is such that `B(d)` matches
`v` (`dseries` objects `A` and `v` must have the same number of
variables).

*Example*

>     >> ts1 = dseries(ones(10,1));
>     >> ts2 = ts1.cumsum();
>     >> ts2
>
>     ts2 is a dseries object:
>
>         | cumsum(Variable_1)
>     1Y  | 1
>     2Y  | 2
>     3Y  | 3
>     4Y  | 4
>     5Y  | 5
>     6Y  | 6
>     7Y  | 7
>     8Y  | 8
>     9Y  | 9
>     10Y | 10
>
>     >> ts3 = ts1.cumsum(dates('3Y'));
>     >> ts3
>
>     ts3 is a dseries object:
>
>         | cumsum(Variable_1)
>     1Y  | -2
>     2Y  | -1
>     3Y  | 0
>     4Y  | 1
>     5Y  | 2
>     6Y  | 3
>     7Y  | 4
>     8Y  | 5
>     9Y  | 6
>     10Y | 7
>
>     >> ts4 = ts1.cumsum(dates('3Y'),dseries(pi));
>     >> ts4
>
>     ts4 is a dseries object:
>
>         | cumsum(Variable_1)
>     1Y  | 1.1416
>     2Y  | 2.1416
>     3Y  | 3.1416
>     4Y  | 4.1416
>     5Y  | 5.1416
>     6Y  | 6.1416
>     7Y  | 7.1416
>     8Y  | 8.1416
>     9Y  | 9.1416
>     10Y | 10.1416
:::

::: {.dseriesmethod}
B = detrend (A, m) [dentrend]() (A, m)

Detrends `dseries` object `A` with a fitted polynomial of order `m`.
Note that each variable is detrended with a different polynomial.
:::

::: {.dseriesmethod}
B = dgrowth (A) [dgrowth]() (A)

Computes daily growth rates.
:::

::: {.dseriesmethod}
B = diff (A) [diff]() (A)

Returns the first difference of `dseries` object `A`.
:::

::: {.datesmethod}
disp (A)

Overloads the MATLAB/Octave disp function for `dseries` object.
:::

::: {.datesmethod}
display (A)

Overloads the MATLAB/Octave display function for `dseries` object.
`display` is the function called by MATLAB to print the content of an
object if a semicolon is missing at the end of a MATLAB statement. If
the `dseries` object is defined over a too large time span, only the
first and last periods will be printed. If the `dseries` object contains
too many variables, only the first and last variables will be printed.
If all the periods and variables are required, the `disp` method should
be used instead.
:::

::: {.dseriesmethod}
C = eq (A, B)

Overloads the MATLAB/Octave `eq` (equal, `==`) operator. `dseries`
objects `A` and `B` must have the same number of observations (say, $T$)
and variables ($N$). The returned argument is a $T \times N$ matrix of
logicals. Element $(i,j)$ of `C` is equal to `true` if and only if
observation $i$ for variable $j$ in `A` and `B` are the same.

*Example*

>     >> ts0 = dseries(2*ones(3,1));
>     >> ts1 = dseries([2; 0; 2]);
>     >> ts0==ts1
>
>     ans =
>
>        3x1 logical array
>
>         1
>         0
>         1
:::

::: {.dseriesmethod}
l = exist (A, varname)

Tests if variable `varname` exists in `dseries` object `A`. Returns
`true` iff variable exists in `A`.

*Example*

>     >> ts = dseries(randn(100,1));
>     >> ts.exist('Variable_1')
>
>     ans =
>
>        logical
>
>         1
>
>     >> ts.exist('Variable_2')
>
>     ans =
>
>        logical
>
>         0
:::

::: {.dseriesmethod}
B = exp (A) [exp]() (A)

Overloads the MATLAB/Octave `exp` function for `dseries` objects.

*Example*

>     >> ts0 = dseries(rand(10,1));
>     >> ts1 = ts0.exp();
:::

::: {.dseriesmethod}
C = extract (A, B\[, \...\])

Extracts some variables from a `dseries` object `A` and returns a
`dseries` object `C`. The input arguments following `A` are strings
representing the variables to be selected in the new `dseries` object
`C`. To simplify the creation of sub-objects, the `dseries` class
overloads the curly braces (`D = extract (A, B, C)` is equivalent to
`D = A{B,C}`) and allows implicit loops (defined between a pair of `@`
symbol, see examples below) or MATLAB/Octave's regular expressions
(introduced by square brackets).

*Example*

> The following selections are equivalent:
>
>     >> ts0 = dseries(ones(100,10));
>     >> ts1 = ts0{'Variable_1','Variable_2','Variable_3'};
>     >> ts2 = ts0{'Variable_@1,2,3@'};
>     >> ts3 = ts0{'Variable_[1-3]$'};
>     >> isequal(ts1,ts2) && isequal(ts1,ts3)
>
>     ans =
>
>        logical
>
>         1
>
> It is possible to use up to two implicit loops to select variables:
>
>     names = {'GDP_1';'GDP_2';'GDP_3'; 'GDP_4'; 'GDP_5'; 'GDP_6'; 'GDP_7'; 'GDP_8'; ...
>         'GDP_9'; 'GDP_10'; 'GDP_11'; 'GDP_12'; ...
>         'HICP_1';'HICP_2';'HICP_3'; 'HICP_4'; 'HICP_5'; 'HICP_6'; 'HICP_7'; 'HICP_8'; ...
>         'HICP_9'; 'HICP_10'; 'HICP_11'; 'HICP_12'};
>
>     ts0 = dseries(randn(4,24),dates('1973Q1'),names);
>     ts0{'@GDP,HICP@_@1,3,5@'}
>
>     ans is a dseries object:
>
>            | GDP_1    | GDP_3     | GDP_5     | HICP_1   | HICP_3   | HICP_5
>     1973Q1 | 1.7906   | -1.6606   | -0.57716  | 0.60963  | -0.52335 | 0.26172
>     1973Q2 | 2.1624   | 3.0125    | 0.52563   | 0.70912  | -1.7158  | 1.7792
>     1973Q3 | -0.81928 | 1.5008    | 1.152     | 0.2798   | 0.88568  | 1.8927
>     1973Q4 | -0.03705 | -0.35899  | 0.85838   | -1.4675  | -2.1666  | -0.62032
:::

::: {.dseriesmethod}
f = firstdate (A)

Returns the first period in `dseries` object `A`.
:::

::: {.dseriesmethod}
f = firstobservedperiod (A)

Returns the first period where all the variables in `dseries` object `A`
are observed (non NaN).
:::

::: {.dseriesmethod}
B = flip (A) [flip]() (A)

Flips the rows in the data member (without changing the periods order).
:::

::: {.dseriesmethod}
f = frequency (B)

Returns the frequency of the variables in `dseries` object `B`.

*Example*

>     >> ts = dseries(randn(3,2),'1973Q1');
>     >> ts.frequency
>
>     ans =
>
>          4
:::

::: {.dseriesmethod}
D = horzcat (A, B\[, \...\])

Overloads the `horzcat` MATLAB/Octave's method for `dseries` objects.
Returns a `dseries` object `D` containing the variables in `dseries`
objects passed as inputs: `A, B, ...` If the inputs are not defined on
the same time ranges, the method adds NaNs to the variables so that the
variables are redefined on the smallest common time range. Note that the
names in the `dseries` objects passed as inputs must be different and
these objects must have common frequency.

*Example*

>     >> ts0 = dseries(rand(5,2),'1950Q1',{'nifnif';'noufnouf'});
>     >> ts1 = dseries(rand(7,1),'1950Q3',{'nafnaf'});
>     >> ts2 = [ts0, ts1];
>     >> ts2
>
>     ts2 is a dseries object:
>
>            | nifnif  | noufnouf | nafnaf
>     1950Q1 | 0.17404 | 0.71431  | NaN
>     1950Q2 | 0.62741 | 0.90704  | NaN
>     1950Q3 | 0.84189 | 0.21854  | 0.83666
>     1950Q4 | 0.51008 | 0.87096  | 0.8593
>     1951Q1 | 0.16576 | 0.21184  | 0.52338
>     1951Q2 | NaN     | NaN      | 0.47736
>     1951Q3 | NaN     | NaN      | 0.88988
>     1951Q4 | NaN     | NaN      | 0.065076
>     1952Q1 | NaN     | NaN      | 0.50946
:::

::: {.dseriesmethod}
B = hpcycle (A\[, lambda\]) [hpcycle]() (A\[, lambda\])

Extracts the cycle component from a `dseries` `A` object using the
*Hodrick and Prescott (1997)* filter and returns a `dseries` object,
`B`. The default value for `lambda`, the smoothing parameter, is `1600`.

*Example*

>     % Simulate a component model (stochastic trend, deterministic
>     % trend, and a stationary autoregressive process).
>     e = 0.2*randn(200,1);
>     u = randn(200,1);
>     stochastic_trend = cumsum(e);
>     deterministic_trend = .1*transpose(1:200);
>     x = zeros(200,1);
>     for i=2:200
>         x(i) = .75*x(i-1) + u(i);
>     end
>     y = x + stochastic_trend + deterministic_trend;
>
>     % Instantiates time series objects.
>     ts0 = dseries(y,'1950Q1');
>     ts1 = dseries(x,'1950Q1'); % stationary component.
>
>     % Apply the HP filter.
>     ts2 = ts0.hpcycle();
>
>     % Plot the filtered time series.
>     plot(ts1(ts2.dates).data,'-k'); % Plot of the stationary component.
>     hold on
>     plot(ts2.data,'--r');           % Plot of the filtered y.
>     hold off
>     axis tight
>     id = get(gca,'XTick');
>     set(gca,'XTickLabel',strings(ts.dates(id)));
:::

::: {.dseriesmethod}
B = hptrend (A\[, lambda\]) [hptrend]() (A\[, lambda\])

Extracts the trend component from a `dseries` A object using the
*Hodrick and Prescott (1997)* filter and returns a `dseries` object,
`B`. Default value for `lambda`, the smoothing parameter, is `1600`.

*Example*

>     % Using the same generating data process
>     % as in the previous example:
>
>     ts1 = dseries(stochastic_trend + deterministic_trend,'1950Q1');
>     % Apply the HP filter.
>     ts2 = ts0.hptrend();
>
>     % Plot the filtered time series.
>     plot(ts1.data,'-k'); % Plot of the nonstationary components.
>     hold on
>     plot(ts2.data,'--r');  % Plot of the estimated trend.
>     hold off
>     axis tight
>     id = get(gca,'XTick');
>     set(gca,'XTickLabel',strings(ts0.dates(id)));
:::

::: {.dseriesmethod}
C = insert (A, B, I)

Inserts variables contained in `dseries` object `B` in `dseries` object
`A` at positions specified by integer scalars in vector `I`, returns
augmented `dseries` object `C`. The integer scalars in `I` must take
values between
`and`A.length()+1`and refers to`A`’s column numbers. The`dseries`objects`A`and`B\`\`
need not be defined over the same time ranges, but it is assumed that
they have common frequency.

*Example*

>     >> ts0 = dseries(ones(2,4),'1950Q1',{'Sly'; 'Gobbo'; 'Sneaky'; 'Stealthy'});
>     >> ts1 = dseries(pi*ones(2,1),'1950Q1',{'Noddy'});
>     >> ts2 = ts0.insert(ts1,3)
>
>     ts2 is a dseries object:
>
>            | Sly | Gobbo | Noddy  | Sneaky | Stealthy
>     1950Q1 | 1   | 1     | 3.1416 | 1      | 1
>     1950Q2 | 1   | 1     | 3.1416 | 1      | 1
>
>     >> ts3 = dseries([pi*ones(2,1) sqrt(pi)*ones(2,1)],'1950Q1',{'Noddy';'Tessie Bear'});
>     >> ts4 = ts0.insert(ts1,[3, 4])
>
>     ts4 is a dseries object:
>
>            | Sly | Gobbo | Noddy  | Sneaky | Tessie Bear | Stealthy
>     1950Q1 | 1   | 1     | 3.1416 | 1      | 1.7725      | 1
>     1950Q2 | 1   | 1     | 3.1416 | 1      | 1.7725      | 1
:::

::: {.dseriesmethod}
B = isempty (A)

Overloads the MATLAB/octave's `isempty` function. Returns `true` if
`dseries` object `A` is empty.
:::

::: {.dseriesmethod}
C = isequal (A, B)

Overloads the MATLAB/octave's `isequal` function. Returns `true` if
`dseries` objects `A` and `B` are identical.
:::

::: {.dseriesmethod}
C = isinf (A)

Overloads the MATLAB/octave's `isinf` function. Returns a logical array,
with element `(i,j)` equal to `true` if and only if variable `j` is
finite in period `A.dates(i)`.
:::

::: {.dseriesmethod}
C = isnan (A)

Overloads the MATLAB/octave's `isnan` function. Returns a logical array,
with element `(i,j)` equal to `true` if and only if variable `j` isn\'t
NaN in period `A.dates(i)`.
:::

::: {.dseriesmethod}
C = isreal (A)

Overloads the MATLAB/octave's `isreal` function. Returns a logical
array, with element `(i,j)` equal to `true` if and only if variable `j`
is real in period `A.dates(i)`.
:::

::: {.dseriesmethod}
B = lag (A\[, p\]) [lag]() (A\[, p\])

Returns lagged time series. Default value of integer scalar `p`, the
number of lags, is `1`.

*Example*

>     >> ts0 = dseries(transpose(1:4), '1950Q1')
>
>     ts0 is a dseries object:
>
>            | Variable_1
>     1950Q1 | 1
>     1950Q2 | 2
>     1950Q3 | 3
>     1950Q4 | 4
>
>     >> ts1 = ts0.lag()
>
>     ts1 is a dseries object:
>
>                | Variable_1
>         1950Q1 | NaN
>         1950Q2 | 1
>         1950Q3 | 2
>         1950Q4 | 3
>
>     >> ts2 = ts0.lag(2)
>
>     ts2 is a dseries object:
>
>            | Variable_1
>     1950Q1 | NaN
>     1950Q2 | NaN
>     1950Q3 | 1
>     1950Q4 | 2
>
>     % dseries class overloads the parenthesis
>     % so that ts.lag(p) can be written more
>     % compactly as ts(-p). For instance:
>
>     >> ts0.lag(1)
>
>     ans is a dseries object:
>
>            | Variable_1
>     1950Q1 | NaN
>     1950Q2 | 1
>     1950Q3 | 2
>     1950Q4 | 3
>
> or alternatively:
>
>     >> ts0(-1)
>
>     ans is a dseries object:
>
>            | Variable_1
>     1950Q1 | NaN
>     1950Q2 | 1
>     1950Q3 | 2
>     1950Q4 | 3
:::

::: {.dseriesmethod}
l = lastdate (B)

Returns the last period in `dseries` object `B`.

*Example*

>     >> ts = dseries(randn(3,2),'1973Q1');
>     >> ts.lastdate()
>
>     ans = <dates: 1973Q3>
:::

::: {.dseriesmethod}
f = lastobservedperiod (A)

Returns the last period where all the variables in `dseries` object `A`
are observed (non NaN).
:::

::: {.dseriesmethod}
B = lead (A\[, p\]) [lead]() (A\[, p\])

Returns lead time series. Default value of integer scalar `p`, the
number of leads, is `1`. As in the `lag` method, the `dseries` class
overloads the parenthesis so that `ts.lead(p)` is equivalent to `ts(p)`.

*Example*

>     >> ts0 = dseries(transpose(1:4),'1950Q1');
>     >> ts1 = ts0.lead()
>
>     ts1 is a dseries object:
>
>            | Variable_1
>     1950Q1 | 2
>     1950Q2 | 3
>     1950Q3 | 4
>     1950Q4 | NaN
>
>     >> ts2 = ts0(2)
>
>     ts2 is a dseries object:
>
>            | Variable_1
>     1950Q1 | 3
>     1950Q2 | 4
>     1950Q3 | NaN
>     1950Q4 | NaN

*Remark*

The overloading of the parenthesis for `dseries` objects, allows to
easily create new `dseries` objects by copying/pasting equations
declared in the `model` block. For instance, if an Euler equation is
defined in the `model` block:

    model;
    ...
    1/C - beta/C(1)*(exp(A(1))*K^(alpha-1)+1-delta) ;
    ...
    end;

and if variables
`,`A`and`K`are defined as`dseries`objects, then by writing::      Residuals = 1/C - beta/C(1)*(exp(A(1))*K^(alpha-1)+1-delta) ;  outside of the`model`block, we create a new`dseries`object, called`Residuals`, for the residuals of the Euler equation (the conditional expectation of the equation defined in the`model\`\`
block is zero, but the residuals are non zero).
:::

::: {.dseriesmethod}
B = lineartrend (A)

Returns a linear trend centered on 0, the length of the trend is given
by the size of `dseries` object `A` (the number of periods).

*Example*

>     >> ts = dseries(ones(3,1));
>     >> ts.lineartrend()
>
>     ans =
>
>          -1
>           0
>           1
:::

::: {.dseriesmethod}
B = log (A) [log]() (A)

Overloads the MATLAB/Octave `log` function for `dseries` objects.

*Example*

>     >> ts0 = dseries(rand(10,1));
>     >> ts1 = ts0.log();
:::

::: {.dseriesmethod}
B = mdiff (A) [mdiff]() (A) B = mgrowth (A) [mgrowth]() (A)

Computes monthly differences or growth rates of variables in `dseries`
object `A`.
:::

::: {.dseriesmethod}
B = mean (A\[, geometric\])

Overloads the MATLAB/Octave `mean` function for `dseries` objects.
Returns the mean of each variable in `dseries` object `A`. If the second
argument is `true` the geometric mean is computed, otherwise (default)
the arithmetic mean is reported.
:::

::: {.dseriesmethod}
C = merge (A, B\[, legacy\])

Merges two `dseries` objects `A` and `B` in `dseries` object `C`.
Objects `A` and `B` need to have common frequency but can be defined on
different time ranges. If a variable, say `x`, is defined both in
`dseries` objects `A` and `B`, then the `merge` will select the variable
`x` as defined in the second input argument, `B`, except for the NaN
elements in `B` if corresponding elements in `A` (ie same periods) are
well defined numbers. This behaviour can be changed by setting the
optional argument `legacy` equal to true, in which case the second
variable overwrites the first one even if the second variable has NaNs.

*Example*

>     >> ts0 = dseries(rand(3,2),'1950Q1',{'A1';'A2'})
>
>     ts0 is a dseries object:
>
>            | A1      | A2
>     1950Q1 | 0.96284 | 0.5363
>     1950Q2 | 0.25145 | 0.31866
>     1950Q3 | 0.34447 | 0.4355
>
>     >> ts1 = dseries(rand(3,1),'1950Q2',{'A1'})
>
>     ts1 is a dseries object:
>
>            | A1
>     1950Q2 | 0.40161
>     1950Q3 | 0.81763
>     1950Q4 | 0.97769
>
>     >> merge(ts0,ts1)
>
>     ans is a dseries object:
>
>            | A1      | A2
>     1950Q1 | 0.96284 | 0.5363
>     1950Q2 | 0.40161 | 0.31866
>     1950Q3 | 0.81763 | 0.4355
>     1950Q4 | 0.97769 | NaN
>
>      >> merge(ts1,ts0)
>
>      ans is a dseries object:
>
>            | A1      | A2
>     1950Q1 | 0.96284 | 0.5363
>     1950Q2 | 0.25145 | 0.31866
>     1950Q3 | 0.34447 | 0.4355
>     1950Q4 | 0.97769 | NaN
:::

::: {.dseriesmethod}
C = minus (A, B)

Overloads the MATLAB/Octave `minus` (`-`) operator for `dseries`
objects, element by element subtraction. If both `A` and `B` are
`dseries` objects, they do not need to be defined over the same time
ranges. If `A` and `B` are `dseries` objects with $T_A$ and $T_B$
observations and $N_A$ and $N_B$ variables, then $N_A$ must be equal to
$N_B$ or $1$ and $N_B$ must be equal to $N_A$ or $1$. If $T_A=T_B$,
`isequal(A.init,B.init)` returns `1` and $N_A=N_B$, then the `minus`
operator will compute for each couple $(t,n)$, with $1\le t\le T_A$ and
$1\le n\le N_A$, `C.data(t,n)=A.data(t,n)-B.data(t,n)`. If $N_B$ is
equal to $1$ and $N_A>1$, the smaller `dseries` object (`B`) is
"broadcast" across the larger `dseries` (`A`) so that they have
compatible shapes, the `minus` operator will subtract the variable
defined in `B` from each variable in `A`. If `B` is a double scalar,
then the method `minus` will subtract `B` from all the
observations/variables in `A`. If `B` is a row vector of length $N_A$,
then the `minus` method will subtract `B(i)` from all the observations
of variable `i`, for $i=1,...,N_A$. If `B` is a column vector of length
$T_A$, then the `minus` method will subtract `B` from all the variables.

*Example*

>     >> ts0 = dseries(rand(3,2));
>     >> ts1 = ts0{'Variable_2'};
>     >> ts0-ts1
>
>     ans is a dseries object:
>
>        | Variable_1 | Variable_2
>     1Y | -0.48853   | 0
>     2Y | -0.50535   | 0
>     3Y | -0.32063   | 0
>
>     >> ts1
>
>     ts1 is a dseries object:
>
>        | Variable_2
>     1Y | 0.703
>     2Y | 0.75415
>     3Y | 0.54729
>
>     >> ts1-ts1.data(1)
>
>     ans is a dseries object:
>
>        | Variable_2
>     1Y | 0
>     2Y | 0.051148
>     3Y | -0.15572
>
>     >> ts1.data(1)-ts1
>
>     ans is a dseries object:
>
>        | Variable_2
>     1Y | 0
>     2Y | -0.051148
>     3Y | 0.15572
:::

::: {.dseriesmethod}
C = mpower (A, B)

Overloads the MATLAB/Octave `mpower` (`^`) operator for `dseries`
objects and computes element-by-element power. `A` is a `dseries` object
with `N` variables and `T` observations. If `B` is a real scalar, then
`mpower(A,B)` returns a `dseries` object `C` with
`C.data(t,n)=A.data(t,n)^C`. If `B` is a `dseries` object with `N`
variables and `T` observations then `mpower(A,B)` returns a `dseries`
object `C` with `C.data(t,n)=A.data(t,n)^C.data(t,n)`.

*Example*

>     >> ts0 = dseries(transpose(1:3));
>     >> ts1 = ts0^2
>
>     ts1 is a dseries object:
>
>        | Variable_1
>     1Y | 1
>     2Y | 4
>     3Y | 9
>
>     >> ts2 = ts0^ts0
>
>     ts2 is a dseries object:
>
>        | Variable_1
>     1Y | 1
>     2Y | 4
>     3Y | 27
:::

::: {.dseriesmethod}
C = mrdivide (A, B)

Overloads the MATLAB/Octave `mrdivide` (`/`) operator for `dseries`
objects, element by element division (like the `./` MATLAB/Octave
operator). If both `A` and `B` are `dseries` objects, they do not need
to be defined over the same time ranges. If `A` and `B` are `dseries`
objects with $T_A$ and $T_B$ observations and $N_A$ and $N_B$ variables,
then $N_A$ must be equal to $N_B$ or $1$ and $N_B$ must be equal to
$N_A$ or $1$. If $T_A=T_B$, `isequal(A.init,B.init)` returns `1` and
$N_A=N_B$, then the `mrdivide` operator will compute for each couple
$(t,n)$, with $1\le t\le T_A$ and $1\le n\le
N_A$, `C.data(t,n)=A.data(t,n)/B.data(t,n)`. If $N_B$ is equal to $1$
and $N_A>1$, the smaller `dseries` object (`B`) is "broadcast" across
the larger `dseries` (`A`) so that they have compatible shapes. In this
case the `mrdivide` operator will divide each variable defined in A by
the variable in B, observation per observation. If B is a double scalar,
then `mrdivide` will divide all the observations/variables in `A` by
`B`. If `B` is a row vector of length $N_A$, then `mrdivide` will divide
all the observations of variable `i` by `B(i)`, for $i=1,...,N_A$. If
`B` is a column vector of length $T_A$, then `mrdivide` will perform a
division of all the variables by `B`, element by element.

*Example*

>     >> ts0 = dseries(rand(3,2))
>
>     ts0 is a dseries object:
>
>        | Variable_1 | Variable_2
>     1Y | 0.72918    | 0.90307
>     2Y | 0.93756    | 0.21819
>     3Y | 0.51725    | 0.87322
>
>     >> ts1 = ts0{'Variable_2'};
>     >> ts0/ts1
>
>     ans is a dseries object:
>
>        | Variable_1 | Variable_2
>     1Y | 0.80745    | 1
>     2Y | 4.2969     | 1
>     3Y | 0.59235    | 1
:::

::: {.dseriesmethod}
C = mtimes (A, B)

Overloads the MATLAB/Octave `mtimes` (`*`) operator for `dseries`
objects and the Hadammard product (the .\* MATLAB/Octave operator). If
both `A` and `B` are `dseries` objects, they do not need to be defined
over the same time ranges. If `A` and `B` are `dseries` objects with
$T_A$ and $_B$ observations and $N_A$ and $N_B$ variables, then $N_A$
must be equal to $N_B$ or $1$ and $N_B$ must be equal to $N_A$ or $1$.
If $T_A=T_B$, `isequal(A.init,B.init)` returns `1` and $N_A=N_B$, then
the `mtimes` operator will compute for each couple $(t,n)$, with
$1\le t\le T_A$ and $1\le n\le
N_A$, `C.data(t,n)=A.data(t,n)*B.data(t,n)`. If $N_B$ is equal to $1$
and $N_A>1$, the smaller `dseries` object (`B`) is "broadcast" across
the larger `dseries` (`A`) so that they have compatible shapes, `mtimes`
operator will multiply each variable defined in `A` by the variable in
`B`, observation per observation. If `B` is a double scalar, then the
method `mtimes` will multiply all the observations/variables in `A` by
`B`. If `B` is a row vector of length $N_A$, then the `mtimes` method
will multiply all the observations of variable `i` by `B(i)`, for
$i=1,...,N_A$. If `B` is a column vector of length $T_A$, then the
`mtimes` method will perform a multiplication of all the variables by
`B`, element by element.
:::

::: {.dseriesmethod}
B = nanmean (A\[, geometric\])

Overloads the MATLAB/Octave `nanmean` function for `dseries` objects.
Returns the mean of each variable in `dseries` object `A` ignoring the
NaN values. If the second argument is `true` the geometric mean is
computed, otherwise (default) the arithmetic mean is reported.
:::

::: {.dseriesmethod}
B = nanstd (A\[, geometric\])

Overloads the MATLAB/Octave `nanstd` function for `dseries` objects.
Returns the standard deviation of each variable in `dseries` object `A`
ignoring the NaN values. If the second argument is `true` the geometric
std is computed, default value of the second argument is `false`.
:::

::: {.dseriesmethod}
C = ne (A, B)

Overloads the MATLAB/Octave `ne` (not equal, `~=`) operator. `dseries`
objects `A` and `B` must have the same number of observations (say, $T$)
and variables ($N$). The returned argument is a $T$ by $N$ matrix of
zeros and ones. Element $(i,j)$ of `C` is equal to `1` if and only if
observation $i$ for variable $j$ in `A` and `B` are not equal.

*Example*

>     >> ts0 = dseries(2*ones(3,1));
>     >> ts1 = dseries([2; 0; 2]);
>     >> ts0~=ts1
>
>     ans =
>
>       3x1 logical array
>
>        0
>        1
>        0
:::

::: {.dseriesmethod}
B = nobs (A)

Returns the number of observations in `dseries` object `A`.

*Example*

>     >> ts0 = dseries(randn(10));
>     >> ts0.nobs
>
>     ans =
>
>         10
:::

::: {.dseriesmethod}
B = onesidedhpcycle (A\[, lambda\[, init\]\]) [onesidedhpcycle]() (A\[,
lambda\[, init\]\])

Extracts the cycle component from a `dseries` `A` object using a one
sided HP filter (with a Kalman filter) and returns a `dseries` object,
`B`. The default value for `lambda`, the smoothing parameter, is `1600`.
By default, if `ìnit` is not provided, the initial value is based on the
first two observations.
:::

::: {.dseriesmethod}
B = onesidedhptrend (A\[, lambda\[, init\]\]) [onesidedhptrend]() (A\[,
lambda\[, init\]\])

Extracts the trend component from a `dseries` `A` object using a one
sided HP filter (with a Kalman filter) and returns a `dseries` object,
`B`. The default value for `lambda`, the smoothing parameter, is `1600`.
By default, if `ìnit` is not provided, the initial value is based on the
first two observations.
:::

::: {.dseriesmethod}
h = plot (A) h = plot (A, B) h = plot (A\[, \...\]) h = plot (A, B\[,
\...\])

Overloads MATLAB/Octave's `plot` function for `dseries` objects. Returns
a MATLAB/Octave plot handle, that can be used to modify the properties
of the plotted time series. If only one `dseries` object, `A`, is passed
as argument, then the plot function will put the associated dates on the
x-abscissa. If this `dseries` object contains only one variable,
additional arguments can be passed to modify the properties of the plot
(as one would do with the MATLAB/Octave's version of the plot function).
If `dseries` object `A` contains more than one variable, it is not
possible to pass these additional arguments and the properties of the
plotted time series must be modified using the returned plot handle and
the MATLAB/Octave `set` function (see example below). If two `dseries`
objects, `A` and `B`, are passed as input arguments, the plot function
will plot the variables in `A` against the variables in `B` (the number
of variables in each object must be the same otherwise an error is
issued). Again, if each object contains only one variable, additional
arguments can be passed to modify the properties of the plotted time
series, otherwise the MATLAB/Octave `set` command has to be used.

*Example*

> Define a `dseries` object with two variables (named by default
> `Variable_1` and `Variable_2`):
>
>     >> ts = dseries(randn(100,2),'1950Q1');
>
> The following command will plot the first variable in `ts`:
>
>     >> plot(ts{'Variable_1'},'-k','linewidth',2);
>
> The next command will draw all the variables in `ts` on the same
> figure:
>
>     >> h = plot(ts);
>
> If one wants to modify the properties of the plotted time series (line
> style, colours, \...), the set function can be used (see MATLAB's
> documentation):
>
>     >> set(h(1),'-k','linewidth',2);
>     >> set(h(2),'--r');
>
> The following command will plot `Variable_1` against
> `exp(Variable_1)`:
>
>     >> plot(ts{'Variable_1'},ts{'Variable_1'}.exp(),'ok');
>
> Again, the properties can also be modified using the returned plot
> handle and the `set` function:
>
>     >> h = plot(ts, ts.exp());
>     >> set(h(1),'ok');
>     >> set(h(2),'+r');
:::

::: {.dseriesmethod}
C = plus (A, B)

Overloads the MATLAB/Octave `plus` (`+`) operator for `dseries` objects,
element by element addition. If both `A` and `B` are `dseries` objects,
they do not need to be defined over the same time ranges. If `A` and `B`
are `dseries` objects with $T_A$ and $T_B$ observations and $N_A$ and
$N_B$ variables, then $N_A$ must be equal to $N_B$ or $1$ and $N_B$ must
be equal to $N_A$ or $1$. If $T_A=T_B$, `isequal(A.init,B.init)` returns
`1` and $N_A=N_B$, then the `plus` operator will compute for each couple
$(t,n)$, with $1\le t\le T_A$ and $1\le n\le N_A$,
`C.data(t,n)=A.data(t,n)+B.data(t,n)`. If $N_B$ is equal to $1$ and
$N_A>1$, the smaller `dseries` object (`B`) is "broadcast" across the
larger `dseries` (`A`) so that they have compatible shapes, the plus
operator will add the variable defined in `B` to each variable in `A`.
If `B` is a double scalar, then the method `plus` will add `B` to all
the observations/variables in `A`. If `B` is a row vector of length
$N_A$, then the `plus` method will add `B(i)` to all the observations of
variable `i`, for $i=1,...,N_A$. If `B` is a column vector of length
$T_A$, then the `plus` method will add `B` to all the variables.
:::

::: {.dseriesmethod}
C = pop (A\[, B\]) [pop]() (A\[, B\])

Removes variable `B` from `dseries` object `A`. By default, if the
second argument is not provided, the last variable is removed.

*Example*

>     >> ts0 = dseries(ones(3,3));
>     >> ts1 = ts0.pop('Variable_2');
>
>     ts1 is a dseries object:
>
>        | Variable_1 | Variable_3
>     1Y | 1          | 1
>     2Y | 1          | 1
>     3Y | 1          | 1
:::

::: {.dseriesmethod}
A = projection (A, info, periods)

Projects variables in dseries object `A`. `info` is is a $n \times 3$
cell array. Each row provides informations necessary to project a
variable. The first column contains the name of variable (row char
array). the second column contains the name of the method used to
project the associated variable (row char array), possible values are
`'Trend'`, `'Constant'`, and `'AR'`. Last column provides quantitative
information about the projection. If the second column value is
`'Trend'`, the third column value is the growth factor of the
(exponential) trend. If the second column value is `'Constant'`, the
third column value is the level of the variable. If the second column
value is `'AR'`, the third column value is the autoregressive parameter.
The variables can be projected with an AR(p) model, if the third column
contains a 1×p vector of doubles. The stationarity of the AR(p) model is
not tested. The case of the constant projection, using the last value of
the variable, is covered with \'Trend\' and a growth factor equal to 1,
or \'AR\' with an autoregressive parameter equal to one (random walk).
This projection routine only deals with exponential trends.

*Example*

>     >> data = ones(10,4);
>     >> ts = dseries(data, '1990Q1', {'A1', 'A2', 'A3', 'A4'});
>     >> info = {'A1', 'Trend', 1.2; 'A2', 'Constant', 0.0; 'A3', 'AR', .5; 'A4', 'AR', [.4, -.2]};
>     >> ts.projection(info, 10);
:::

::: {.dseriesmethod}
B = qdiff (A) B = qgrowth (A) [qdiff]() (A) [qgrowth]() (A)

Computes quarterly differences or growth rates.

*Example*

>     >> ts0 = dseries(transpose(1:4),'1950Q1');
>     >> ts1 = ts0.qdiff()
>
>     ts1 is a dseries object:
>
>            | Variable_1
>     1950Q1 | NaN
>     1950Q2 | 1
>     1950Q3 | 1
>     1950Q4 | 1
>
>     >> ts0 = dseries(transpose(1:6),'1950M1');
>     >> ts1 = ts0.qdiff()
>
>     ts1 is a dseries object:
>
>             | Variable_1
>     1950M1  | NaN
>     1950M2  | NaN
>     1950M3  | NaN
>     1950M4  | 3
>     1950M5  | 3
>     1950M6  | 3
:::

::: {.dseriesmethod}
C = remove (A, B) [remove]() (A, B)

Alias for the `pop` method with two arguments. Removes variable `B` from
`dseries` object `A`.

*Example*

>     >> ts0 = dseries(ones(3,3));
>     >> ts1 = ts0.remove('Variable_2');
>
>     ts1 is a dseries object:
>
>        | Variable_1 | Variable_3
>     1Y | 1          | 1
>     2Y | 1          | 1
>     3Y | 1          | 1
>
> A shorter syntax is available: `remove(ts,'Variable_2')` is equivalent
> to `ts{'Variable_2'} = []` (`[]` can be replaced by any empty object).
> This alternative syntax is useful if more than one variable has to be
> removed. For instance:
>
>     ts{'Variable_@2,3,4@'} = [];
>
> will remove `Variable_2`, `Variable_3` and `Variable_4` from `dseries`
> object `ts` (if these variables exist). Regular expressions cannot be
> used but implicit loops can.
:::

::: {.dseriesmethod}
B = rename (A, oldname, newname) [rename]() (A, oldname, newname)

Rename variable `oldname` to `newname` in `dseries` object `A`. Returns
a `dseries` object. If more than one variable needs to be renamed, it is
possible to pass cells of char arrays as second and third arguments.

*Example*

>     >> ts0 = dseries(ones(2,2));
>     >> ts1 = ts0.rename('Variable_1','Stinkly')
>
>     ts1 is a dseries object:
>
>        | Stinkly | Variable_2
>     1Y | 1       | 1
>     2Y | 1       | 1
:::

::: {.dseriesmethod}
C = rename (A, newname) [rename]() (A, newname)

Replace the names in `A` with those passed in the cell string array
`newname`. `newname` must have the same number of elements as `dseries`
object `A` has variables. Returns a `dseries` object.

*Example*

>     >> ts0 = dseries(ones(2,3));
>     >> ts1 = ts0.rename({'TinkyWinky','Dipsy','LaaLaa'})
>
>     ts1 is a dseries object:
>
>        | TinkyWinky | Dipsy | LaaLaa
>     1Y | 1          | 1     | 1
>     2Y | 1          | 1     | 1
:::

::: {.dseriesmethod}
A = resetops (A, ops)

Redefine `ops` member.
:::

::: {.dseriesmethod}
A = resetags (A, ops)

Redefine `tags` member.
:::

::: {.dseriesmethod}
B = round (A\[, n\]) [round]() (A\[, n\])

Rounds to the nearest decimal or integer. `n` is the precision parameter
(number of decimals), default value is 0 meaning that that by default
the method rounds to the nearest integer.

*Example*

>     >> ts = dseries(pi)
>
>     ts is a dseries object:
>
>        | Variable_1
>     1Y | 3.1416
>
>     >> ts.round_();
>     >> ts
>
>     ts is a dseries object:
>
>        | Variable_1
>     1Y | 3
:::

::: {.dseriesmethod}
save (A, basename\[, format\])

Overloads the MATLAB/Octave `save` function and saves `dseries` object
`A` to disk. Possible formats are `mat` (this is the default), `m`
(MATLAB/Octave script), and `csv` (MATLAB binary data file). The name of
the file without extension is specified by `basename`.

*Example*

>     >> ts0 = dseries(ones(2,2));
>     >> ts0.save('ts0', 'csv');
>
> The last command will create a file ts0.csv with the following
> content:
>
>     ,Variable_1,Variable_2
>     1Y,               1,               1
>     2Y,               1,               1
>
> To create a MATLAB/Octave script, the following command:
>
>     >> ts0.save('ts0','m');
>
> will produce a file ts0.m with the following content:
>
>     % File created on 14-Nov-2013 12:08:52.
>
>     FREQ__ = 1;
>     INIT__ = ' 1Y';
>
>     NAMES__ = {'Variable_1'; 'Variable_2'};
>     TEX__ = {'Variable_{1}'; 'Variable_{2}'};
>     OPS__ = {};
>     TAGS__ = struct();
>
>     Variable_1 = [
>                   1
>                   1];
>
>     Variable_2 = [
>                   1
>                   1];
>
> The generated (`csv`, `m`, or `mat`) files can be loaded when
> instantiating a `dseries` object as explained above.
:::

::: {.dseriesmethod}
B = set\_names(A, s1, s2, \...)

Renames variables in `dseries` object `A` and returns a `dseries` object
`B` with new names `s1`, `s2`, \... The number of input arguments after
the first one (`dseries` object `A`) must be equal to `A.vobs` (the
number of variables in `A`). `s1` will be the name of the first variable
in `B`, `s2` the name of the second variable in `B`, and so on.

*Example*

>     >> ts0 = dseries(ones(1,3));
>     >> ts1 = ts0.set_names('Barbibul',[],'Barbouille')
>
>     ts1 is a dseries object:
>
>        | Barbibul | Variable_2 | Barbouille
>     1Y | 1        | 1          | 1
:::

::: {.dseriesmethod}
\[T, N \] = size(A\[, dim\])

Overloads the MATLAB/Octave's `size` function. Returns the number of
observations in `dseries` object `A` (i.e. `A.nobs`) and the number of
variables (i.e. `A.vobs`). If a second input argument is passed, the
`size` function returns the number of observations if `dim=1` or the
number of variables if `dim=2` (for all other values of `dim` an error
is issued).

*Example*

>     >> ts0 = dseries(ones(1,3));
>     >> ts0.size()
>
>     ans =
>
>          1     3
:::

::: {.dseriesmethod}
B = std (A\[, geometric\])

Overloads the MATLAB/Octave `std` function for `dseries` objects.
Returns the standard deviation of each variable in `dseries` object `A`.
If the second argument is `true` the geometric standard deviation is
computed (default value of the second argument is `false`).
:::

::: {.dseriesmethod}
B = subsample (A, d1, d2)

Returns a subsample, for periods between `dates` `d1` and `d2`. The same
can be achieved by indexing a `dseries` object with a `dates` object,
but the `subsample` method is easier to use programmatically.

*Example*

>     >> o = dseries(transpose(1:5));
>     >> o.subsample(dates('2y'),dates('4y'))
>
>     ans is a dseries object:
>
>        | Variable_1
>     2Y | 2
>     3Y | 3
>     4Y | 4
:::

::: {.dseriesmethod}
A = tag (A, a\[, b, c\])

Add a tag to a variable in `dseries` object `A`.

*Example*

>     >> ts = dseries(randn(10, 3));
>     >> tag(ts, 'type');             % Define a tag name.
>     >> tag(ts, 'type', 'Variable_1', 'Stock');
>     >> tag(ts, 'type', 'Variable_2', 'Flow');
>     >> tag(ts, 'type', 'Variable_3', 'Stock');
:::

::: {.dseriesmethod}
B = tex\_rename (A, name, newtexname) B = tex\_rename (A, newtexname)
[tex\_rename]() (A, name, newtexname) [tex\_rename]() (A, newtexname)

Redefines the tex name of variable `name` to `newtexname` in `dseries`
object `A`. Returns a `dseries` object.

With only two arguments `A` and `newtexname`, it redefines the tex names
of the `A` to those contained in `newtexname`. Here, `newtexname` is a
cell string array with the same number of entries as variables in `A`.
:::

::: {.dseriesmethod}
B = uminus(A)

Overloads `uminus` (`-`, unary minus) for `dseries` object.

*Example*

>     >> ts0 = dseries(1)
>
>     ts0 is a dseries object:
>
>        | Variable_1
>     1Y | 1
>
>     >> ts1 = -ts0
>
>     ts1 is a dseries object:
>
>        | Variable_1
>     1Y | -1
:::

::: {.dseriesmethod}
D = vertcat (A, B\[, \...\])

Overloads the `vertcat` MATLAB/Octave method for `dseries` objects. This
method is used to append more observations to a `dseries` object.
Returns a `dseries` object `D` containing the variables in `dseries`
objects passed as inputs. All the input arguments must be `dseries`
objects with the same variables defined on different time ranges.

*Example*

>     >> ts0 = dseries(rand(2,2),'1950Q1',{'nifnif';'noufnouf'});
>     >> ts1 = dseries(rand(2,2),'1950Q3',{'nifnif';'noufnouf'});
>     >> ts2 = [ts0; ts1]
>
>     ts2 is a dseries object:
>
>            | nifnif   | noufnouf
>     1950Q1 | 0.82558  | 0.31852
>     1950Q2 | 0.78996  | 0.53406
>     1950Q3 | 0.089951 | 0.13629
>     1950Q4 | 0.11171  | 0.67865
:::

::: {.dseriesmethod}
B = vobs (A)

Returns the number of variables in `dseries` object `A`.

*Example*

>     >> ts0 = dseries(randn(10,2));
>     >> ts0.vobs
>
>     ans =
>
>         2
:::

::: {.dseriesmethod}
B = ydiff (A) B = ygrowth (A) [ydiff]() (A) [ygrowth]() (A)

Computes yearly differences or growth rates.
:::
:::

X-13 ARIMA-SEATS interface {#x13-members}
--------------------------

::: {.x13}
The x13 class provides a method for each X-13 command as documented in
the X-13 ARIMA-SEATS reference manual ([x11]{.title-ref},
[automdl]{.title-ref}, [estimate]{.title-ref}, \...), options can then
be passed by key/value pairs. The `x13` class has 22 members:

arg y

:   `dseries` object with a single variable.

arg x

:   `dseries` object with an arbitrary number of variables (to be used
    in the REGRESSION block).

arg arima

:   structure containing the options of the ARIMA model command.

arg automdl

:   structure containing the options of the ARIMA model selection
    command.

arg regression

:   structure containing the options of the Regression command.

arg estimate

:   structure containing the options of the estimation command.

arg transform

:   structure containing the options of the transform command.

arg outlier

:   structure containing the options of the outlier command.

arg forecast

:   structure containing the options of the forecast command.

arg check

:   structure containing the options of the check command.

arg x11

:   structure containing the options of the X11 command.

arg force

:   structure containing the options of the force command.

arg history

:   structure containing the options of the history command.

arg metadata

:   structure containing the options of the metadata command.

arg identify

:   structure containing the options of the identify command.

arg pickmdl

:   structure containing the options of the pickmdl command.

arg seats

:   structure containing the options of the seats command.

arg slidingspans

:   structure containing the options of the slidingspans command.

arg spectrum

:   structure containing the options of the spectrum command.

arg x11regression

:   structure containing the options of the x11Regression command.

arg results

:   structure containing the results returned by x13.

arg commands

:   cell array containing the list of commands.

All these members are private. The following constructors are available:

::: {.construct}
x13 (y)

Instantiates an `x13` object with [dseries]{.title-ref} object `y`. The
`dseries` object passed as an argument must contain only one variable,
the one we need to pass to X-13.
:::

::: {.construct}
x13 (y, x)

Instantiates an `x13` object with [dseries]{.title-ref} objects `y` and
`x`. The first `dseries` object passed as an argument must contain only
one variable, the second `dseries` object contains the exogenous
variables used by some of the X-13 commands. Both objects must be
defined on the same time span.
:::

The Following methods allow to set sequence of X-13 commands, write an
[.spc]{.title-ref} file and run the X-13 binary:

*Example*

>     >> ts = dseries(rand(100,1),'1999M1');
>     >> o = x13(ts);
>
>     >> o.x11('save','(d11)');
>     >> o.automdl('savelog','amd','mixed','no');
>     >> o.outlier('types','all','save','(fts)');
>     >> o.check('maxlag',24,'save','(acf pcf)');
>     >> o.estimate('save','(mdl est)');
>     >> o.forecast('maxlead',18,'probability',0.95,'save','(fct fvr)');
>
>     >> o.run();
:::

Miscellaneous
-------------

### Time aggregation

> A set of functions allows to convert time series to lower frequencies:
>
> > -   `dseries2M` converts daily time series object to monthly time
> >     series object.
> > -   `dseries2Q` converts daily or monthly time series object to
> >     quarterly time series object.
> > -   `dseries2S` converts daily, monthly, or quarterly time series
> >     object to bi-annual time series object.
> > -   `dseries2Y` converts daily, monthly, quarterly, or bi-annual
> >     time series object to annual time series object.
>
> All these routines have two mandatory input arguments: the first one
> is a `dseries` object, the second one the name (row char array) of the
> aggregation method. Possible values for the second argument are:
>
> > -   `arithmetic-average` (for growth rates),
> > -   `geometric-average` (for growth factors),
> > -   `sum` (for flow variables), and
> > -   `end-of-period` (for stock variables).
>
> *Example*
>
> >     >> ts = dseries(rand(12,1),'2000M1')
> >
> >     ts is a dseries object:
> >
> >             | Variable_1
> >     2000M1  | 0.55293
> >     2000M2  | 0.14228
> >     2000M3  | 0.38036
> >     2000M4  | 0.39657
> >     2000M5  | 0.57674
> >     2000M6  | 0.019402
> >     2000M7  | 0.57758
> >     2000M8  | 0.9322
> >     2000M9  | 0.10687
> >     2000M10 | 0.73215
> >     2000M11 | 0.97052
> >     2000M12 | 0.60889
> >
> >     >> ds = dseries2Y(ts, 'end-of-period')
> >
> >     ds is a dseries object:
> >
> >           | Variable_1
> >     2000Y | 0.60889

### Create time series with a univariate model

> It is possible to expand a `dseries` object recursively with the
> `from` command. For instance to create a `dseries` object containing
> the simulation of an ARMA(1,1) model:
>
> >     >> e = dseries(randn(100, 1), '2000Q1', 'e', '\varepsilon');
> >     >> y = dseries(zeros(100, 1), '2000Q1', 'y');
> >     >> from 2000Q2 to 2024Q4 do y(t)=.9*y(t-1)+e(t)-.4*e(t-1);
> >     >> y
> >
> >     y is a dseries object:
> >
> >            | y
> >     2000Q1 | 0
> >     2000Q2 | -0.95221
> >     2000Q3 | -0.6294
> >     2000Q4 | -1.8935
> >     2001Q1 | -1.1536
> >     2001Q2 | -1.5905
> >     2001Q3 | 0.97056
> >     2001Q4 | 1.1409
> >     2002Q1 | -1.9255
> >     2002Q2 | -0.29287
> >            |
> >     2022Q2 | -1.4683
> >     2022Q3 | -1.3758
> >     2022Q4 | -1.2218
> >     2023Q1 | -0.98145
> >     2023Q2 | -0.96542
> >     2023Q3 | -0.23203
> >     2023Q4 | -0.34404
> >     2024Q1 | 1.4606
> >     2024Q2 | 0.901
> >     2024Q3 | 2.4906
> >     2024Q4 | 0.79661
>
> The expression following the `do` keyword can be any univariate
> equation, the only constraint is that the model cannot have leads. It
> can be a static equation, or a very nonlinear backward equation with
> an arbitrary number of lags. The `from` command must be followed by a
> range, which is separated from the (recursive) expression to be
> evaluated by the `do` command.
