# Results generated by Dynare

Dynare stores results in the `context` structure. The main entry
points are
- `trends`: `context.results.model_results[1].trends`
- `LREresults`: `context.results.model_results[1].linearrationalresults`
| Object             | Storage                        | Generated by     | Accessor | Remark |
|--------------------|--------------------------------|------------------|----------|--------|
| steady state       | trends.endogenous_steady_state | Steady_State.jl  |          |        |
|                    | trends.exogenous_steady_state  |                  |          |        |
| 1st order solution | LREresults.g1                  | perturbations.jl |          |        |
|                    | LREresults.g1_1                |                  |          |        |
|                    | LREresults.g1_2                |                  |          |        |
|                    | LREresults.gns1                |                  |          |        |
|                    | LREresults.gs1                 |                  |          |        |
|                    | LREresults.hns1                |                  |          |        |
|                    | LREresults.hs1                 |                  |          |        |
