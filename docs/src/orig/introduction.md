::: {.default-domain}
dynare
:::

Introduction
============

What is Dynare?
---------------

Dynare is a software platform for handling a wide class of economic
models, in particular dynamic stochastic general equilibrium (DSGE) and
overlapping generations (OLG) models. The models solved by Dynare
include those relying on the *rational expectations* hypothesis, wherein
agents form their expectations about the future in a way consistent with
the model. But Dynare is also able to handle models where expectations
are formed differently: on one extreme, models where agents perfectly
anticipate the future; on the other extreme, models where agents have
limited rationality or imperfect knowledge of the state of the economy
and, hence, form their expectations through a learning process. In terms
of types of agents, models solved by Dynare can incorporate consumers,
productive firms, governments, monetary authorities, investors and
financial intermediaries. Some degree of heterogeneity can be achieved
by including several distinct classes of agents in each of the
aforementioned agent categories.

Dynare offers a user-friendly and intuitive way of describing these
models. It is able to perform simulations of the model given a
calibration of the model parameters and is also able to estimate these
parameters given a dataset. In practice, the user will write a text file
containing the list of model variables, the dynamic equations linking
these variables together, the computing tasks to be performed and the
desired graphical or numerical outputs.

A large panel of applied mathematics and computer science techniques are
internally employed by Dynare: multivariate nonlinear solving and
optimization, matrix factorizations, local functional approximation,
Kalman filters and smoothers, MCMC techniques for Bayesian estimation,
graph algorithms, optimal control, ...

Various public bodies (central banks, ministries of economy and finance,
international organisations) and some private financial institutions use
Dynare for performing policy analysis exercises and as a support tool
for forecasting exercises. In the academic world, Dynare is used for
research and teaching purposes in postgraduate macroeconomics courses.

Dynare is a free software, which means that it can be downloaded free of
charge, that its source code is freely available, and that it can be
used for both non-profit and for-profit purposes. Most of the source
files are covered by the GNU General Public Licence (GPL) version 3 or
later (there are some exceptions to this, see the file license.txt in
Dynare distribution). It is available for the Windows, macOS, and Linux
platforms and is fully documented through a reference manual. Part of
Dynare is programmed in C++, while the rest is written using the
[MATLAB](https://www.mathworks.com/products/matlab/) programming
language. The latter implies that commercially-available MATLAB software
is required in order to run Dynare. However, as an alternative to
MATLAB, Dynare is also able to run on top of [GNU
Octave](https://www.octave.org/) (basically a free clone of MATLAB):
this possibility is particularly interesting for students or
institutions who cannot afford, or do not want to pay for, MATLAB and
are willing to bear the concomitant performance loss.

The development of Dynare is mainly done at
[CEPREMAP](https://www.cepremap.fr/) by a core team of researchers who
devote part of their time to software development. Increasingly, the
developer base is expanding, as tools developed by researchers outside
of CEPREMAP are integrated into Dynare. Financial support is provided by
CEPREMAP, Banque de France and DSGE-net (an international research
network for DSGE modeling).

Interaction between developers and users of Dynare is central to the
project. A [web forum](https://forum.dynare.org/) is available for users
who have questions about the usage of Dynare or who want to report bugs.
Current known and fixed bugs are listed on the [Dynare
wiki](https://git.dynare.org/Dynare/dynare/wikis). Issues or whishes can
be reported on our [Git
repository](https://git.dynare.org/Dynare/dynare). Training sessions are
given through the Dynare Summer School, which is organized every year
and is attended by about 40 people. Finally, priorities in terms of
future developments and features to be added are decided in cooperation
with the institutions providing financial support.

Documentation sources
---------------------

The present document is the reference manual for Dynare. It documents
all commands and features in a systematic fashion.

Other useful sources of information include the [Dynare
wiki](https://git.dynare.org/Dynare/dynare/wikis) and the [Dynare
forums](https://forum.dynare.org/).

Citing Dynare in your research
------------------------------

You should cite Dynare if you use it in your research. The recommended
way todo this is to cite the present manual, as:

> Stéphane Adjemian, Houtan Bastani, Michel Juillard, Frédéric Karamé,
> Ferhat Mihoubi, Willi Mutschler, Johannes Pfeifer, Marco Ratto,
> Normann Rion and Sébastien Villemot (2022), "Dynare: Reference Manual,
> Version 5," *Dynare Working Papers*, 72, CEPREMAP

For convenience, you can copy and paste the following into your BibTeX
file:

> ``` {.sourceCode .bibtex}
> @TechReport{Adjemianetal2022,
>   author      = {Adjemian, St\'ephane and Bastani, Houtan and
>                  Juillard, Michel and Karam\'e, Fr\'ederic and
>                  Mihoubi, Ferhat and Mutschler, Willi
>                  and Pfeifer, Johannes and Ratto, Marco and
>                  Rion, Normann and Villemot, S\'ebastien},
>   title       = {Dynare: Reference Manual Version 5},
>   year        = {2022},
>   institution = {CEPREMAP},
>   type        = {Dynare Working Papers},
>   number      = {72},
> }
> ```

If you want to give a URL, use the address of the Dynare website:
<https://www.dynare.org>.
