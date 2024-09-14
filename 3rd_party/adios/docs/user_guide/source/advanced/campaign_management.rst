####################
Campaign Management
####################

The campaign management in ADIOS2 is for collecting basic information and metadata about a collection of ADIOS2 output files, from a single application run or multiple runs. The campaign archive is a single file (.ACA) that can be transferred to other locations. The campaign file can be opened by ADIOS2 and all the metadata can be processed (including the values of GlobalValue and LocalValue variables, or min/max of each Arrays at each step and decomposition/min/max of each block in an Array at each step). However, Get() operations will only succeed to read actual data of the arrays, if the data belonging to the campaign is either local or some mechanism for remote data access to the location of the data is set up in advance.

.. warning::

    In 2.10, Campaign Management is just a first prototype and is included only for evaluation purposes. It will change substantially in the future and campaign files produced by this version will unlikely to be supported going forward. 

The idea
========
Applications produce one or more output files in a single run. Subsequent analysis and visualization runs produce more output files. Campaign is a data organization concept one step higher than a file. A campaign archive includes information about multiple files, including the scalar variable's values and the min/max of arrays and the location of the data files (host and directory information). A science project can agree on how to organize their campaigns, i.e., how to name them, what files to include in a single campaign archive, how to distribute them, how to name the hosts where the actual data resides. 

Example
-------

The Gray-Scott example, that is included with ADIOS2, in `examples/simulation/gray-scott`, has two programs, Gray-Scott and PDF-Calc. The first one produces the main output `gs.bp` which includes the main 3D variables `U` and `V`, and a checkpoint file `ckpt.bp` with a single step in it. PDF-Calc processes the main output and produces histograms on 2D slices of U and V (`U/bins` and `U/pdf`) in `pdf.bp`. A campaign can include all the three output files as they logically belong together. 

.. code-block:: bash
		
    # run application as usual
    $ mpirun -n 4 adios2_simulations_gray-scott settings-files.json
    $ ls -d *.bp
    ckpt.bp gs.bp

    $ adios2_campaign_manager.py create demoproject/frontier_gray-scott_100
    
    $ mpirun -n 3 adios2_simulations_gray-scott_pdf-calc gs.bp pdf.bp 1000
    $ ls -d *.bp
    ckpt.bp gs.bp pdf.bp

    $ adios2_campaign_manager.py update demoproject/frontier_gray-scott_100

    $ adios2_campaign_manager.py info  demoproject/frontier_gray-scott_100
    info archive
    ADIOS Campaign Archive, version 1.0, created on 2024-04-01 10:44:11.644942
    hostname = OLCF   longhostname = frontier.olcf.ornl.gov
        dir = /lustre/orion/csc143/proj-shared/demo/gray-scott
            dataset = ckpt.bp     created on 2024-04-01 10:38:19
            dataset = gs.bp     created on 2024-04-01 10:38:17
            dataset = pdf.bp     created on 2024-04-01 10:38:08

    # The campaign archive is small compared to the data it points to 
    $ du -sh *bp
    7.9M    ckpt.bp
    40M     gs.bp
    9.9M    pdf.bp

    $ du -sh /lustre/orion/csc143/proj-shared/adios-campaign-store/demoproject/frontier_gray-scott_100.aca
    97K     /lustre/orion/csc143/proj-shared/adios-campaign-store/demoproject/frontier_gray-scott_100.aca

    # ADIOS can list the content of the campaign archive
    $ bpls -l demoproject/frontier_gray-scott_100
        double   ckpt.bp/U      {4, 34, 34, 66} = 0.171103 / 1
        double   ckpt.bp/V      {4, 34, 34, 66} = 1.71085e-19 / 0.438921
        int32_t  ckpt.bp/step   scalar = 700
        double   gs.bp/U        10*{64, 64, 64} = 0.090778 / 1
        double   gs.bp/V        10*{64, 64, 64} = 8.24719e-63 / 0.515145
        int32_t  gs.bp/step     10*scalar = 100 / 1000
        double   pdf.bp/U/bins  10*{1000} = 0.0908158 / 0.999938
        double   pdf.bp/U/pdf   10*{64, 1000} = 0 / 4096
        double   pdf.bp/V/bins  10*{1000} = 8.24719e-63 / 0.514267
        double   pdf.bp/V/pdf   10*{64, 1000} = 0 / 4096
        int32_t  pdf.bp/step    10*scalar = 100 / 1000

    # scalar over steps is available in metadata
    $ bpls -l demoproject/frontier_gray-scott_100 -d pdf.bp/step -n 10
      int32_t  pdf.bp/step    10*scalar = 100 / 1000
        (0)    100 200 300 400 500 600 700 800 900 1000

    # Array decomposition including min/max are available in metadata
    $ bpls -l demoproject/frontier_gray-scott_100 -D gs.bp/V
      double   gs.bp/V        10*{64, 64, 64} = 8.24719e-63 / 0.515145
        step 0:
          block 0: [ 0:63,  0:31,  0:31] = 8.24719e-63 / 0.410653
          block 1: [ 0:63, 32:63,  0:31] = 8.24719e-63 / 0.410652
          block 2: [ 0:63,  0:31, 32:63] = 8.24719e-63 / 0.410653
          block 3: [ 0:63, 32:63, 32:63] = 8.24719e-63 / 0.410653
        ...
        step 9:
          block 0: [ 0:63,  0:31,  0:31] = 3.99908e-09 / 0.441847
          block 1: [ 0:63, 32:63,  0:31] = 3.99931e-09 / 0.44192
          block 2: [ 0:63,  0:31, 32:63] = 3.99928e-09 / 0.441813
          block 3: [ 0:63, 32:63, 32:63] = 3.99899e-09 / 0.441796

    # Array data is only available if data is local
    $ ./bin/bpls -l demoproject/frontier_gray-scott_100 -d pdf.bp/U/bins
      double   pdf.bp/U/bins  10*{1000} = 0.0908158 / 0.999938
        (0,  0)    0.93792 0.937982 0.938044 0.938106 0.938168 0.93823 0.938292 0.938354 0.938416 0.938479
        ...
        (9,990)    0.990306 0.991157 0.992007 0.992858 0.993708 0.994559 0.995409 0.99626 0.99711 0.997961


Setup
=====

There are three paths/names important in the campaign setup. 

- `hostname` is the name detected by the adios2_campaign_manager when creating a campaign archive, however, it is better to define a specific name the project agrees upon (e.g. OLCF, NERSC, ALCF) that identifies the generic location of the data and then use that name later to specify the modes of remote data access (not available in this release).

- `campaignstorepath` is the directory where all the campaign archives are stored. This should be shared between project members in a center, and a private one on every member's laptop. It is up to the project to determine what file sharing / synchronization mechanism to use to sync this directories. `Rclone is a great command-line tool <https://rclone.org>`_ to sync the campaign store with many cloud-based file sharing services and cloud instances.

- `cachepath` is the directory where ADIOS can unpack metadata from the campaign archive so that ADIOS engines can read them as if they were entirely local datasets. The cache only contains the metadata for now but in the future data that have already been retrieved by previous read requests will be stored here as well. 


Use `~/.config/adios2/adios2.yaml` to specify these options. 

.. code-block:: bash
		
    $ cat ~/.config/adios2/adios2.yaml

    Campaign:
      active: true
      hostname: OLCF
      campaignstorepath: /lustre/orion/csc143/proj-shared/adios-campaign-store
      cachepath: /lustre/orion/csc143/proj-shared/campaign-cache
      verbose: 0

    $ ls -R ~/dropbox/adios-campaign-store
    /lustre/orion/csc143/proj-shared/adios-campaign-store/demoproject:
    frontier_gray-scott_100.aca

    $ adios2_campaign_manager.py list
    demoproject/frontier_gray-scott_100.aca


Remote access
=============
For now, we have one way to access data, through SSH port forwarding and running a remote server program to read in data on the remote host and to send back the data to the local ADIOS program. `adios2_remote_server` is included in the adios installation. You need to use the one built on the host.

Launch the server by SSH-ing to the remote machine, and specifying the `26200` port for fowarding. For example:

.. code-block:: bash
		
    $ ssh -L 26200:dtn.olcf.ornl.gov:26200 -l <username> dtn.olcf.ornl.gov "<path_to_adios_install>/bin/adios2_remote_server -v "

Assuming the campaign archive was synced to a local machine's campaign store under `csc143/demoproject`, now we can retrieve data:

.. code-block:: bash

    $ adios2_campaign_manager.py list
    csc143/demoproject/frontier_gray-scott_100.aca

    $ bpls -l csc143/demoproject/frontier_gray-scott_100 
      double   ckpt.bp/U      {4, 34, 34, 66} = 0.171103 / 1
      ...
      double   pdf.bp/U/bins  10*{1000} = 0.0908158 / 0.999938

    # metadata is extracted to the local cachepath
    $ du -sh /tmp/campaign/OLCF/csc143/demoproject/frontier_gray-scott_100.aca/*
    20K     /tmp/campaign/OLCF/csc143/demoproject/frontier_gray-scott_100.aca/ckpt.bp
    40K     /tmp/campaign/OLCF/csc143/demoproject/frontier_gray-scott_100.aca/gs.bp
    32K     /tmp/campaign/OLCF/csc143/demoproject/frontier_gray-scott_100.aca/pdf.bp

    # data is requested from the remote server
    # read 16 values (4x4x4) from U from last step, from offset 30,30,30
    $ bpls -l csc143/demoproject/frontier_gray-scott_100  -d gs.bp/U -s "-1,30,30,30" -c "1,4,4,4" -n 4
    double   gs.bp/U        10*{64, 64, 64}
      slice (9:9, 30:33, 30:33, 30:33)
      (9,30,30,30)    0.89189 0.899854 0.899854 0.891891
      (9,30,31,30)    0.899851 0.908278 0.908278 0.899852
      (9,30,32,30)    0.899849 0.908276 0.908277 0.899851
      (9,30,33,30)    0.891885 0.899848 0.899849 0.891886
      (9,31,30,30)    0.89985 0.908276 0.908276 0.899849
      (9,31,31,30)    0.908274 0.916977 0.916977 0.908274
      (9,31,32,30)    0.908273 0.916976 0.916976 0.908273
      (9,31,33,30)    0.899844 0.908271 0.908271 0.899844
      (9,32,30,30)    0.89985 0.908276 0.908275 0.899848
      (9,32,31,30)    0.908274 0.916976 0.916976 0.908272
      (9,32,32,30)    0.908272 0.916975 0.916974 0.908271
      (9,32,33,30)    0.899844 0.90827 0.90827 0.899842
      (9,33,30,30)    0.89189 0.899851 0.899851 0.891886
      (9,33,31,30)    0.89985 0.908275 0.908275 0.899847
      (9,33,32,30)    0.899848 0.908274 0.908273 0.899845
      (9,33,33,30)    0.891882 0.899845 0.899844 0.89188

Requirements
============
The Campaign Manager uses SQlite3 and ZLIB for its operations, and Python3 3.8 or higher for the `adios2_campaign_manager` tool. Check `bpls -Vv` to see if `CAMPAIGN` is in the list of "Available features".

Limitations
===========

- The Campaign Reader engine only supports ReadRandomAccess mode, not step-by-step reading. Campaign management will need to change in the future to support sorting the steps from different outputs to a coherent order. 
- Updates to moving data for other location is not supported yet
