# Demonstrated features 


|                         | GAB | KTC | LMA | TPF | MVC | HMI | EDN | CHT |
|-------------------------|-----|-----|-----|-----|-----|-----|-----|-----|
| lowMach                 |     |     |  x  |     |  x  |     |     |     |
| varying p0th            |     |     |     |     |  x  |     |     |     |
| recycling               |     |     |     |  x  |     |     |     |     |
| source term             |  x  |  x  |  x  |  x  |  x  |     |     |     |
| scalar transport        |  x  |  x  |  x  |     |  x  |     |     |  x  |
| variable props          |  x  |  x  |  x  |     |  x  |     |     |     |
| Dirichlet BC            |  x  |     |  x  |     |     |  x  |  x  |  x  |
| Flux (Neumann) BC       |  x  |     |     |     |     |     |     |  x  |
| user BC data (usrwrk)   |  x  |     |     |     |     |     |     |     |
| traction BC             |  x  |     |     |     |     |     |     |     |
| sym BC                  |     |  x  |     |     |  x  |  x  |     |     |
| turbulent outflow       |     |     |     |     |     |  x  |     |     |
| variable dt             |  x  |     |     |  x  |     |     |     |     |
| OIFS                    |  x  |     |     |  x  |     |     |     |     |
| Lagrangian particles    |     |     |     |     |     |  x  |     |     |
| Point interpolation     |     |     |     |     |     |  x  |     |     |
| usrchk postprocessing   |  x  |     |     |     |     |     |     |     |
| runtime averages        |  x  |     |     |     |     |     |     |     |
| planar average          |  x  |     |     |     |     |     |     |     |
| conjugate heat transfer |     |     |     |     |     |     |     |  x  |
| RANS (k-tau)            |     |  x  |     |     |     |     |     |     |
| drag                    |     |  x  |     |     |     |     |     |     |
| mesh manipulation       |     |     |     |     |     |     |     |     |
| moving mesh (ALE)       |     |     |     |     |  x  |  x  |     |     |
| constant flow rate      |     |  x  |     |     |     |     |     |     |
| overset grids (neknek)  |     |     |     |     |     |     |  x  |     |
| par casedata            |  x  |     |     |  x  |  x  |  x  |  x  |     |
| nek data exchange       |  x  |     |     |     |     |     |     |     |
| hpf-RT                  |  x  |     |     |  x  |     |  x  |     |     |
| avm (scalar)            |     |     |     |  x  |     |     |     |     |
| time subcycling         |     |     |     |     |     |     |  x  |     |


### Ledgend
`GAB`: gabls1

`KTC`: ktauChannel

`LMA`: lowMach

`TPF`: turbPipe

`MVC`: mv_cyl

`HMI`: hemi

`EDN`: eddyNekNek                

`CHT`: conj_ht                
