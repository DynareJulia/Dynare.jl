::: {.default-domain}
dynare
:::

The configuration file {#conf-file}
======================

The configuration file is used to provide Dynare with information not
related to the model (and hence not placed in the model file). At the
moment, it is only used when using Dynare to run parallel computations.

On Linux and macOS, the default location of the configuration file is
`$HOME/.dynare`, while on Windows it is `%APPDATA%\dynare.ini`
(typically `c:\Users\USERNAME\AppData\dynare.ini`). You can specify a
non standard location using the `conffile` option of the `dynare`
command (see `dyn-invoc`{.interpreted-text role="ref"}).

The parsing of the configuration file is case-sensitive and it should
take the following form, with each option/choice pair placed on a
newline:

    [command0]
    option0 = choice0
    option1 = choice1

    [command1]
    option0 = choice0
    option1 = choice1

The configuration file follows a few conventions (self-explanatory
conventions such as `USER_NAME` have been excluded for concision):

`COMPUTER_NAME`

> Indicates the valid name of a server (e.g. `localhost`,
> `server.cepremap.org`) or an IP address.

`DRIVE_NAME`

> Indicates a valid drive name in Windows, without the trailing colon
> (e.g. `C`).

`PATH`

> Indicates a valid path in the underlying operating system (e.g.
> `/home/user/dynare/matlab/`).

`PATH_AND_FILE`

> Indicates a valid path to a file in the underlying operating system
> (e.g. `/usr/local/MATLAB/R2010b/bin/matlab`).

`BOOLEAN`

> Is `true` or `false`.

Dynare Configuration
--------------------

This section explains how to configure Dynare for general processing.
Currently, there is only one option available.

::: {.confblock}
\[hooks\]

This block can be used to specify configuration options that will be
used when running Dynare.

*Options*

::: {.option}
GlobalInitFile = PATH\_AND\_FILE

The location of the global initialization file to be run at the end of
`global_initialization.m`.
:::

*Example*

>     [hooks]
>     GlobalInitFile = /home/usern/dynare/myInitFile.m
:::

::: {.confblock}
\[paths\]

This block can be used to specify paths that will be used when running
dynare.

*Options*

::: {.option}
Include = PATH

A colon-separated path to use when searching for files to include via
`@#include`. Paths specified via `-I
<-I\<\<path\>\>>`{.interpreted-text role="opt"} take priority over paths
specified here, while these paths take priority over those specified by
`@#includepath`.
:::

*Example*

>     [paths]
>     Include = /path/to/folder/containing/modfiles:/path/to/another/folder
:::

Parallel Configuration {#paral-conf}
----------------------

This section explains how to configure Dynare for parallelizing some
tasks which require very little inter-process communication.

The parallelization is done by running several MATLAB or Octave
processes, either on local or on remote machines. Communication between
leader and follower processes are done through SMB on Windows and SSH on
UNIX. Input and output data, and also some short status messages, are
exchanged through network filesystems. Currently the system works only
with homogenous grids: only Windows or only Unix machines.

The following routines are currently parallelized:

> -   the posterior sampling algorithms when using multiple chains;
> -   the Metropolis-Hastings diagnostics;
> -   the posterior IRFs;
> -   the prior and posterior statistics;
> -   some plotting routines.

Note that creating the configuration file is not enough in order to
trigger parallelization of the computations: you also need to specify
the `parallel` option to the `dynare` command. For more details, and for
other options related to the parallelization engine, see
`dyn-invoc`{.interpreted-text role="ref"}.

You also need to verify that the following requirements are met by your
cluster (which is composed of a leader and of one or more followers):

For a Windows grid:

> -   a standard Windows network (SMB) must be in place;
> -   the
>     [PsTools](https://technet.microsoft.com/sysinternals/pstools.aspx)
>     suite must be installed in the path of the leader Windows machine;
> -   the Windows user on the leader machine has to be user of any other
>     follower machine in the cluster, and that user will be used for
>     the remote computations.
> -   detailed step-by-step setup instructions can be found in
>     `win-ssg`{.interpreted-text role="ref"}.

For a UNIX grid:

> -   SSH must be installed on the leader and on the follower machines;
> -   SSH keys must be installed so that the SSH connection from the
>     leader to the follower can be done without passwords, or using an
>     SSH agent.

::: {.warning}
::: {.admonition-title}
Warning
:::

Compatibility considerations between leader and follower

It is highly recommended to use the same version of Dynare on both the
leader and all followers. Different versions regularly cause problems
like zero acceptance rates during estimation. When upgrading to a newer
Dynare version do not forget to adjust the `DynarePath`.
:::

We now turn to the description of the configuration directives. Note
that comments in the configuration file can be provided by separate
lines starting with a hashtag (\#).

::: {.confblock}
\[cluster\]

When working in parallel, `[cluster]` is required to specify the group
of computers that will be used. It is required even if you are only
invoking multiple processes on one computer.

*Options*

::: {.option}
Name = CLUSTER\_NAME

The reference name of this cluster.
:::

::: {.option}
Members = NODE\_NAME\[(WEIGHT)\] NODE\_NAME\[(WEIGHT)\] \...

A list of nodes that comprise the cluster with an optional computing
weight specified for that node. The computing weight indicates how much
more powerful one node is with respect to the others (e.g.
`n1(2) n2(1) n3(3)` means that `n1` is two times more powerful than `n2`
whereas `n3` is three times more powerful than `n2`). Each node is
separated by at least one space and the weights are in parenthesis with
no spaces separating them from their node.
:::

*Example*

>     [cluster]
>     Name = c1
>     Members = n1 n2 n3
>
>     [cluster]
>     Name = c2
>     Members = n1(4) n2 n3
:::

::: {.confblock}
\[node\]

When working in parallel, `[node]` is required for every computer that
will be used. The options that are required differ, depending on the
underlying operating system and whether you are working locally or
remotely.

*Options*

::: {.option}
Name = NODE\_NAME

The reference name of this node.
:::

::: {.option}
CPUnbr = INTEGER \| \[INTEGER:INTEGER\]

If just one integer is passed, the number of processors to use. If a
range of integers is passed, the specific processors to use (processor
counting is defined to begin at one as opposed to zero). Note that using
specific processors is only possible under Windows; under Linux and
macOS, if a range is passed the same number of processors will be used
but the range will be adjusted to begin at one.
:::

::: {.option}
ComputerName = COMPUTER\_NAME

The name or IP address of the node. If you want to run locally, use
`localhost` (case-sensitive).
:::

::: {.option}
Port = INTEGER

The port number to connect to on the node. The default is empty, meaning
that the connection will be made to the default SSH port (22).
:::

::: {.option}
UserName = USER\_NAME

The username used to log into a remote system. Required for remote runs
on all platforms.
:::

::: {.option}
Password = PASSWORD

The password used to log into the remote system. Required for remote
runs originating from Windows.
:::

::: {.option}
RemoteDrive = DRIVE\_NAME

The drive to be used for remote computation. Required for remote runs
originating from Windows.
:::

::: {.option}
RemoteDirectory = PATH

The directory to be used for remote computation. Required for remote
runs on all platforms.
:::

::: {.option}
DynarePath = PATH

The path to the matlab subdirectory within the Dynare installation
directory. The default is the empty string.
:::

::: {.option}
MatlabOctavePath = PATH\_AND\_FILE

The path to the MATLAB or Octave executable. The default value is
`matlab`.
:::

::: {.option}
NumberOfThreadsPerJob = INTEGER

For Windows nodes, sets the number of threads assigned to each remote
MATLAB/Octave run. The default value is 1.
:::

::: {.option}
SingleCompThread = BOOLEAN

Whether or not to disable MATLAB's native multithreading. The default
value is `false`. Option meaningless under Octave.
:::

::: {.option}
OperatingSystem = OPERATING\_SYSTEM

The operating system associated with a node. Only necessary when
creating a cluster with nodes from different operating systems. Possible
values are `unix` or `windows`. There is no default value.
:::

*Example*

>     [node]
>     Name = n1
>     ComputerName = localhost
>     CPUnbr = 1
>
>     [node]
>     Name = n2
>     ComputerName = dynserv.cepremap.org
>     CPUnbr = 5
>     UserName = usern
>     RemoteDirectory = /home/usern/Remote
>     DynarePath = /home/usern/dynare/matlab
>     MatlabOctavePath = matlab
>
>     [node]
>     Name = n3
>     ComputerName = dynserv.dynare.org
>     Port = 3333
>     CPUnbr = [2:4]
>     UserName = usern
>     RemoteDirectory = /home/usern/Remote
>     DynarePath = /home/usern/dynare/matlab
>     MatlabOctavePath = matlab
:::

Windows Step-by-Step Guide {#win-ssg}
--------------------------

This section outlines the steps necessary on most Windows systems to set
up Dynare for parallel execution.

> 1.  Write a configuration file containing the options you want. A
>     mimimum working example setting up a cluster consisting of two
>     local CPU cores that allows for e.g. running two Monte Carlo
>     Markov Chains in parallel is shown below.
> 2.  Save the configuration file somwhere. The name and file ending do
>     not matter if you are providing it with the `conffile` command
>     line option. The only restrictions are that the path must be a
>     valid filename, not contain non-alpha-numeric characters, and not
>     contain any whitespaces. For the configuration file to be
>     accessible without providing an explicit path at the command line,
>     you must save it under the name `dynare.ini` into your user
>     account's `Application Data` folder.
> 3.  Install
>     [PSTools](https://technet.microsoft.com/sysinternals/pstools.aspx)
>     to your system, e.g. into `C:\PSTools.`
> 4.  Set the Windows System Path to the `PSTools` folder (e.g. using
>     something along the line of pressing Windows Key+Pause to open the
>     System Configuration, then go to Advanced -\> Environment
>     Variables -\> Path).
> 5.  Restart your computer to make the path change effective.
> 6.  Open MATLAB and type into the command window:
>
>         !psexec
>
>     This executes the `psexec.exe` from PSTools on your system and
>     shows whether Dynare will be able to locate it. If MATLAB
>     complains at this stage, you did not correctly set your Windows
>     system path for the `PSTools` folder.
>
> 7.  If `psexec.exe` was located in the previous step, a popup will
>     show up, asking for confirmation of the license agreement. Confirm
>     this copyright notice of `psexec` (this needs to be done only
>     once). After this, Dynare should be ready for parallel execution.
> 8.  Call Dynare on your mod-file invoking the `parallel` option and
>     providing the path to your configuration file with the `conffile`
>     option (if you did not save it as `%APPDATA%\dynare.ini` in step 2
>     where it should be detected automatically):
>
>         dynare ls2003 parallel conffile='C:\Users\Dynare~1\parallel\conf_file.ini'
>
> Please keep in mind that no white spaces or names longer than 8
> characters are allowed in the `conffile` path. The 8-character
> restriction can be circumvented by using the tilde Windows path
> notation as in the above example.

*Example*:

    #cluster needs to always be defined first
    [cluster]
    #Provide a name for the cluster
    Name=Local
    #declare the nodes being member of the cluster
    Members=n1

    #declare nodes (they need not all be part of a cluster)
    [node]
    #name of the node
    Name=n1
    #name of the computer (localhost for the current machine)
    ComputerName=localhost
    #cores to be included from this node
    CPUnbr=[1:2]
    #path to matlab.exe; on Windows, the MATLAB bin folder is in the system path
    #so we only need to provide the name of the exe file
    MatlabOctavePath=matlab
    #Dynare path you are using
    DynarePath=C:/dynare/4.7.0/matlab
