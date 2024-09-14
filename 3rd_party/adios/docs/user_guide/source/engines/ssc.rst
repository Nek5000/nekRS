**************************
SSC Strong Staging Coupler
**************************

The SSC engine is designed specifically for strong code coupling. Currently SSC only supports fixed IO pattern, which means once the first step is finished, users are not allowed to write or read a data block with a *start* and *count* that have not been written or read in the first step. SSC uses a combination of one sided MPI and two sided MPI methods. In any cases, all user applications are required to be launched within a single mpirun or mpiexec command, using the MPMD mode.

The SSC engine takes the following parameters:

1. ``OpenTimeoutSecs``: Default **10**. Timeout in seconds for opening a stream. The SSC engine's open function will block until the RendezvousAppCount is reached, or timeout, whichever comes first. If it reaches the timeout, SSC will throw an exception.

2. ``Threading``: Default **False**. SSC will use threads to hide the time cost for metadata manipulation and data transfer when this parameter is set to **true**. SSC will check if MPI is initialized with multi-thread enabled, and if not, then SSC will force this parameter to be **false**. Please do NOT enable threading when multiple I/O streams are opened in an application, as it will cause unpredictable errors. This parameter is only effective when writer definitions and reader selections are NOT locked. For cases definitions and reader selections are locked, SSC has a more optimized way to do data transfers, and thus it will not use this parameter.

=============================== ================== ================================================
 **Key**                         **Value Format**   **Default** and Examples
=============================== ================== ================================================
 OpenTimeoutSecs                        integer            **10**, 2, 20, 200
 Threading                              bool               **false**, true
=============================== ================== ================================================


