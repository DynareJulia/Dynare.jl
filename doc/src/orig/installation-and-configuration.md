::: {.default-domain}
dynare
:::

Installation and configuration
==============================

Software requirements
---------------------

Packaged versions of Dynare are available for Windows (8.1, 10 and 11),
several GNU/Linux distributions (Debian, Ubuntu, Linux Mint, Arch
Linux), macOS (12 "Monterey"), and FreeBSD. Dynare should work on other
systems, but some compilation steps are necessary in that case.

In order to run Dynare, you need one of the following:

-   MATLAB, any version ranging from 8.3 (R2014a) to 9.12 (R2022a);
-   GNU Octave, any version ranging from 6.2.0 to 7.1.0, with the
    statistics package from
    [Octave-Forge](https://octave.sourceforge.io/). Note however that
    the Dynare installer for Windows requires a more specific version of
    Octave, as indicated on the download page.

The following optional extensions are also useful to benefit from extra
features, but are in no way required:

-   If under MATLAB: the Optimization Toolbox, the Statistics Toolbox,
    the Control System Toolbox;
-   If under Octave, the following
    [Octave-Forge](https://octave.sourceforge.io/) packages:
    `optim, io, control`.

Installation of Dynare
----------------------

After installation, Dynare can be used in any directory on your
computer. It is best practice to keep your model files in directories
different from the one containing the Dynare toolbox. That way you can
upgrade Dynare and discard the previous version without having to worry
about your own files.

### On Windows

Execute the automated installer called `dynare-x.y-win.exe` (where `x.y`
is the version number), and follow the instructions. The default
installation directory is `c:\dynare\x.y`.

After installation, this directory will contain several sub-directories,
among which are `matlab`, `mex` and `doc`.

The installer will also add an entry in your Start Menu with a shortcut
to the documentation files and uninstaller.

Note that you can have several versions of Dynare coexisting (for
example in `c:\dynare`), as long as you correctly adjust your path
settings (see see `words-warning`{.interpreted-text role="ref"}).

Also note that it is possible to do a silent installation, by passing
the `/S` flag to the installer on the command line. This can be useful
when doing an unattended installation of Dynare on a computer pool.

### On GNU/Linux

On Debian, Ubuntu and Linux Mint, the Dynare package can be installed
with: `apt install dynare`. This will give a fully-functional Dynare
installation usable with Octave. If you have MATLAB installed, you
should also do: `apt install dynare-matlab` (under Debian, this package
is in the `contrib` section). Documentation can be installed with
`apt install dynare-doc`. The status of those packages can be checked at
those pages:

-   [Package status in Debian](https://packages.debian.org/sid/dynare)
-   [Package status in
    Ubuntu](https://launchpad.net/ubuntu/+source/dynare)
-   [Package status in Linux
    Mint](https://community.linuxmint.com/software/view/dynare)

On Arch Linux, the Dynare package is not in the official repositories,
but is available in the [Arch User
Repository](https://wiki.archlinux.org/index.php/Arch_User_Repository).
The needed sources can be downloaded from the [package status in Arch
Linux](https://aur.archlinux.org/packages/dynare/).

Dynare will be installed under `/usr/lib/dynare`. Documentation will be
under `/usr/share/doc/dynare` (only on Debian, Ubuntu and Linux Mint).

### On macOS

#### With MATLAB

To install Dynare for use with MATLAB, execute the automated installer
called `dynare-x.y.pkg` (where *x.y* is the version number), and follow
the instructions. The default installation directory is
`/Applications/Dynare/x.y`. After installation, this directory will
contain several sub-directories, among which are `matlab`, `mex`, and
`doc`.

Note that several versions of Dynare can coexist (by default in
`/Applications/Dynare`), as long as you correctly adjust your path
settings (see `words-warning`{.interpreted-text role="ref"}).

By default, the installer installs a version of GCC (for use with
`use_dll`{.interpreted-text role="opt"}) in the installation directory,
under the `.brew` folder. To do so, it also installs a version of
[Homebrew](https://brew.sh) in the same folder and Xcode Command Line
Tools (this is an Apple product) in a system folder.

All of this requires a bit of time and hard disk space. The amount of
time it takes will depend on your computing power and internet
connection. To reduce the time the Dynare installer takes, you can
install Xcode Command Line Tools yourself (see
`prerequisites-macos`{.interpreted-text role="ref"}). Dynare, Homebrew,
and GCC use about 600 MB of disk space while the Xcode Command Line
Tools require about 400 MB.

If you do not use the `use_dll`{.interpreted-text role="opt"} option,
you have the choice to forgo the installation of GCC and hence Dynare
will only take about 50 MB of disk space.

#### With Octave

We don't provide Dynare packages for macOS with Octave support, but
there is a Dynare package with Octave support in
[Homebrew](https://brew.sh).

Once [Homebrew](https://brew.sh) is installed, run a terminal and
install Dynare (and Octave) by typing the following:

    brew install dynare

Then open Octave by running the following in the same terminal:

    octave --gui

Finally, at the Octave prompt, install some add-ons (you only have to do
it once):

    octave:1> pkg install -forge io statistics control struct optim

### On FreeBSD

A [FreeBSD port for Dynare](https://www.freshports.org/science/dynare/)
is available. It can be installed with:

    pkg install dynare

### For other systems

You need to download Dynare source code from the [Dynare
website](https://www.dynare.org/) and unpack it somewhere.

Then you will need to recompile the pre-processor and the dynamic
loadable libraries. Please refer to
[README.md](https://git.dynare.org/Dynare/dynare/blob/master/README.md).

Compiler installation {#compil-install}
---------------------

### Prerequisites on Windows

There are no prerequisites on Windows. Dynare now ships a compilation
environment that can be used with the `use_dll`{.interpreted-text
role="opt"} option.

### Prerequisites on GNU/Linux

Users of MATLAB under GNU/Linux need a working compilation environment
installed. Under Debian, Ubuntu or Linux Mint, it can be installed via
`apt install build-essential`.

Users of Octave under GNU/Linux should install the package for MEX file
compilation (under Debian, Ubuntu or Linux Mint, it can be done via
`apt install liboctave-dev`).

### Prerequisites on macOS {#prerequisites-macos}

#### With MATLAB

Dynare now ships a compilation environment that can be used with the
`use_dll`{.interpreted-text role="opt"} option. To install this
environment correctly, the Dynare installer ensures that the Xcode
Command Line Tools (an Apple product) have been installed on a system
folder. To install the Xcode Command Line Tools yourself, simply type
`xcode-select --install` into the Terminal
(`/Applications/Utilities/Terminal.app`) prompt.

#### With Octave

The compiler can be installed via [Homebrew](https://brew.sh). In a
terminal, run:

    brew install gcc-11

Configuration
-------------

### For MATLAB

::: {.highlight}
matlab
:::

You need to add the `matlab` subdirectory of your Dynare installation to
MATLAB path. You have two options for doing that:

-   Using the `addpath` command in the MATLAB command window:

    Under Windows, assuming that you have installed Dynare in the
    standard location, and replacing `x.y` with the correct version
    number, type:

        >> addpath c:/dynare/x.y/matlab

    Under GNU/Linux, type:

        >> addpath /usr/lib/dynare/matlab

    Under macOS, assuming that you have installed Dynare in the standard
    location, and replacing `x.y` with the correct version number, type:

        >> addpath /Applications/Dynare/x.y/matlab

    MATLAB will not remember this setting next time you run it, and you
    will have to do it again.

-   Via the menu entries:

    Select the "Set Path" entry in the "File" menu, then click on "Add
    Folder...", and select the `matlab` subdirectory of 'your Dynare
    installation. Note that you *should not* use "Add with
    Subfolders...". Apply the settings by clicking on "Save". Note that
    MATLAB will remember this setting next time you run it.

### For Octave

You need to add the `matlab` subdirectory of your Dynare installation to
Octave path, using the `addpath` at the Octave command prompt.

Under Windows, assuming that you have installed Dynare in the standard
location, and replacing "*x.y*" with the correct version number, type:

    octave:1> addpath c:/dynare/x.y/matlab

Under Debian, Ubuntu or Linux Mint, there is no need to use the
`addpath` command; the packaging does it for you. Under Arch Linux, you
need to do:

    octave:1> addpath /usr/lib/dynare/matlab

Under macOS, assuming you have installed Dynare via
[Homebrew](https://brew.sh):

    octave:1> addpath /usr/local/lib/dynare/matlab

If you don't want to type this command every time you run Octave, you
can put it in a file called `.octaverc` in your home directory (under
Windows this will generally be `c:\Users\USERNAME` while under macOS it
is `/Users/USERNAME/`). This file is run by Octave at every startup.

### Some words of warning {#words-warning}

You should be very careful about the content of your MATLAB or Octave
path. You can display its content by simply typing `path` in the command
window.

The path should normally contain system directories of MATLAB or Octave,
and some subdirectories of your Dynare installation. You have to
manually add the `matlab` subdirectory, and Dynare will automatically
add a few other subdirectories at runtime (depending on your
configuration). You must verify that there is no directory coming from
another version of Dynare than the one you are planning to use.

You have to be aware that adding other directories (on top of the dynare
folders) to your MATLAB or Octave path can potentially create problems
if any of your M-files have the same name as a Dynare file. Your routine
would then override the Dynare routine, making Dynare unusable.

::: {.warning}
::: {.admonition-title}
Warning
:::

Never add all the subdirectories of the `matlab` folder to the MATLAB or
Octave path. You must let Dynare decide which subdirectories have to be
added to the MATLAB or Octave path. Otherwise, you may end up with a non
optimal or un-usable installation of Dynare.
:::
