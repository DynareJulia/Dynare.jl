# The Dynare Julia Reference Manual
# Introduction 

DynareJulia is a rewriting of Dynare (https://www.dynare.org) that was initially written in
Gauss in 1994 and rewritten in Matlab around 2000.

Dynare provides several algorithms to work with Dynamic Stochastic
General Equilibrium (DSGE) models often used in macroeconomics. Among
other features, it helps
 - solving such models,
 - simulating them,
 - estimating the parameters,
 - making forecasts.
 
The user of the package writes a text file, usually with an `.mod`
extension, that contains the equations of the model and the
computation tasks. Then, DynareJulia compiles the model and runs the computations.

DynareJulia honors a subset of commands valid in DynareMatlab. Tell us if one of your favorite command or option is missing.

For many computing tasks, DynareJulia provides also Julia functions that can be used in the `*.mod` file or issued interactively after having run the `*.mod` file. These Julia functions use keyword arguments for the options and you need only to enter them if you want to change the default value. The keyword arguments without a default value are required arguments. In the sections of this documentation, the Dynare Commands are presented first, then the Julia functions.

Dynare has benefited from many contributions over the years. Here is a list of the contributors:

## The Dynare Team

Currently the development team of Dynare is composed of:

-   Stéphane Adjemian (Le Mans Université, Gains)
-   Michel Juillard (Banque de France)
-   Sumudu Kankanamge (Toulouse School of Economics and CEPREMAP)
-   Frédéric Karamé (Le Mans Université, Gains and CEPREMAP)
-   Junior Maih (Norges Bank)
-   Willi Mutschler (University of Tübingen)
-   Johannes Pfeifer (Universität der Bundeswehr München)
-   Marco Ratto (European Commission, Joint Research Centre - JRC)
-   Normann Rion (CY Cergy Paris Université and CEPREMAP)
-   Sébastien Villemot (CEPREMAP)

The following people contribute or have contributed to  DynareJulia
-   Satyanarayana Bade
-   Petre Caraiani
-   Lilith Hafner
-   Michel Juillard
-   Félix Ordoñez
-   Louis Ponet
-   Rohit Singh Rathaur
-   Dawie van Lill

The following people used to be members of the Dynare team:

-   Houtan Bastani
-   Abdeljabar Benzougar
-   Alejandro Buesa
-   Fabrice Collard
-   Assia Ezzeroug
-   Dóra Kocsis
-   Stéphane Lhuissier
-   Ferhat Mihoubi
-   George Perendia

Copyright © 1996-2023, Dynare Team.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3 or
any later version published by the Free Software Foundation; with no
Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

A copy of the license can be found at
<https://www.gnu.org/licenses/fdl.txt>.



