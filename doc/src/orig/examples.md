::: {.default-domain}
dynare
:::

Examples
========

Dynare comes with a database of example `.mod` files, which are designed
to show a broad range of Dynare features, and are taken from academic
papers for most of them. You should have these files in the `examples`
subdirectory of your distribution.

Here is a short list of the examples included. For a more complete
description, please refer to the comments inside the files themselves.

`ramst.mod`

> An elementary real business cycle (RBC) model, simulated in a
> deterministic setup.

`example1.mod` `example2.mod`

> Two examples of a small RBC model in a stochastic setup, presented in
> *Collard (2001)* (see the file `guide.pdf` which comes with Dynare).

`example3.mod`

> A small RBC model in a stochastic setup, presented in *Collard
> (2001)*. The steady state is solved analytically using the
> `steady_state_model` block (see `steady_state_model`{.interpreted-text
> role="bck"}).

`fs2000.mod`

> A cash in advance model, estimated by *Schorfheide (2000)*. The file
> shows how to use Dynare for estimation.

`fs2000_nonstationary.mod`

> The same model than `fs2000.mod`, but written in non-stationary form.
> Detrending of the equations is done by Dynare.

`bkk.mod`

> Multi-country RBC model with time to build, presented in *Backus,
> Kehoe and Kydland (1992)*. The file shows how to use Dynare's macro
> processor.

`agtrend.mod`

> Small open economy RBC model with shocks to the growth trend,
> presented in *Aguiar and Gopinath (2004)*.

`Gali_2015.mod`

> Basic New Keynesian model of *Galí (2015)*, Chapter 3 showing how to
> i) use \"system prior\"-type prior restrictions as in *Andrle and
> Plašil (2018)* and ii) run prior/posterior-functions.

`NK_baseline.mod`

> Baseline New Keynesian Model estimated in *Fernández-Villaverde
> (2010)*. It demonstrates how to use an explicit steady state file to
> update parameters and call a numerical solver.

`Occbin_example.mod`

> RBC model with two occasionally binding constraints. Demonstrates how
> to set up Occbin.

`Ramsey_Example.mod`

> File demonstrating how to conduct optimal policy experiments in a
> simple New Keynesian model either under commitment (Ramsey) or using
> optimal simple rules (OSR)

`Ramsey_steady_file.mod`

> File demonstrating how to conduct optimal policy experiments in a
> simple New Keynesian model under commitment (Ramsey) with a
> user-defined conditional steady state file
