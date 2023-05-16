# Demonstrated features 

|                          | GAB | KTC | LMA | TPF | MVC | HMI | EDN | CHT |
|--------------------------|-----|-----|-----|-----|-----|-----|-----|-----|
| lowMach                  |     |     |  x  |     |  x  |     |     |     |
| varying p0th             |     |     |     |     |  x  |     |     |     |
| inflow recycling         |     |     |     |  x  |     |     |     |     |
| user source term         |  x  |  x  |  x  |  x  |  x  |     |     |     |
| scalar transport         |  x  |  x  |  x  |     |  x  |     |     |  x  |
| variable props           |  x  |  x  |  x  |     |  x  |     |     |     |
| Dirichlet BC             |  x  |     |  x  |     |     |  x  |  x  |  x  |
| flux (Neumann) BC        |  x  |     |     |     |     |     |     |  x  |
| traction BC              |  x  |     |     |     |     |     |     |     |
| sym BC                   |     |  x  |     |     |  x  |  x  |     |     |
| user BC data (usrwrk)    |  x  |     |     |     |     |     |     |     |
| turbulent outflow (Dong) |     |     |     |     |     |  x  |     |     |
| variable dt              |  x  |     |     |  x  |     |     |     |     |
| OIFS/subcycling          |  x  |     |     |  x  |     |     |  x  |     |
| Lagrangian particles     |     |     |     |     |     |  x  |     |     |
| point interpolation      |  x  |  x  |     |     |     |  x  |     |     |
| extract line (hpts)      |     |  x  |     |     |     |     |     |     |
| runtime averages         |  x  |     |     |     |     |     |     |     |
| planar average           |  x  |     |     |     |     |     |     |     |
| conjugate heat transfer  |     |     |     |     |     |     |     |  x  |
| RANS (k-tau)             |     |  x  |     |     |     |     |     |     |
| surfaceIntegral          |     |  x  |     |     |     |     |     |     |
| viscous drag             |     |  x  |     |     |     |     |     |     |
| mesh manipulation        |  x  |     |     |     |     |     |     |     |
| moving mesh (ALE)        |     |     |     |     |  x  |  x  |     |     |
| constant flow rate       |     |  x  |     |     |     |     |     |     |
| overset grids (neknek)   |     |     |     |     |     |     |  x  |     |
| par casedata             |  x  |     |     |  x  |  x  |  x  |  x  |     |
| nek data exchange        |  x  |     |     |     |     |     |     |     |
| hpf-RT                   |  x  |     |     |  x  |     |  x  |     |     |
| avm (scalar)             |     |     |     |  x  |     |     |     |     |
| predictor-corrector iter |     |     |     |     |     |     |  x  |     |
| usrchk postprocessing    |     |     |  x  |     |     |     |     |     |

### Ledgend
`GAB`: gabls1

`KTC`: ktauChannel

`LMA`: lowMach

`TPF`: turbPipe

`MVC`: mv_cyl

`HMI`: hemi

`EDN`: eddyNekNek                

`CHT`: conj_ht                
