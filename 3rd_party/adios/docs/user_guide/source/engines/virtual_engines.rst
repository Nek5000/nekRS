##############################
Supported Virtual Engine Names
##############################

This section provides a description of the Virtual Engines that can be used to set up an actual Engine with specific parameters. 
These virtual names are used for beginner users to simplify the selection of an engine and its parameters. 
The following I/O uses cases are supported by virtual engine names:

1. ``File``: File I/O (Default engine).

   This sets up the I/O for files. If the file name passed in Open() ends with ".bp", then the BP5 engine will be used starting in v2.9.0.
   If it ends with ".h5", the HDF5 engine will be used. For old .bp files (BP version 3 format), the BP3 engine 
   will be used for reading (v2.4.0 and below). 

2. ``FileStream``: Online processing via files.

   This allows a Consumer to concurrently read the data while the Producer is writing new output steps into it. The Consumer will
   wait for the appearance of the file itself in Open() (for up to one hour) and wait for the appearance of new steps in the file
   (in BeginStep() up to the specificed timeout in that function). 

3. ``InSituAnalysis``: Streaming data to another application. 

   This sets up ADIOS for transferring data from a Producer to a Consumer application. The Producer and Consumer are synchronized
   at Open(). The Consumer will receive every single output step from the Producer, therefore, the Producer will
   block on output if the Consumer is slow.

4. ``InSituVisualization``:: Streaming data to another application without waiting for consumption.

   This sets up ADIOS for transferring data from a Producer to a Consumer without ever blocking the Producer. The Producer will
   throw away all output steps that are not immediately requested by a Consumer. It will also not wait for a Consumer to connect. 
   This kind of streaming is great for an interactive visualization session where the user wants to see the most current state of the 
   application.

5. ``CodeCoupling``:: Streaming data between two applications for code coupling. 

   Producer and Consumer are waiting for each other in Open() and every step must be consumed. 
   Currently, this is the same as in situ analysis.

These virtual engine names are used to select a specific engine and its parameters. In practice, after selecting the virtual engine name, 
one can modify the settings by adding/overwriting parameters. Eventually, a seasoned user would use the actual Engine names and parameterize 
it for the specific run. 


Virtual Engine Setups
---------------------

These are the actual settings in ADIOS when a virtual engine is selected. The parameters below can be modified before the Open call. 

1. ``File``. Refer to the parameter settings for these engines of
   ``BP5``, ``BP4``, ``BP3`` and ``HDF5`` engines earlier in this section. 

2. ``FileStream``. The engine is ``BP5``. The parameters are set to:

============================== ===================== ===========================================================
 **Key**                       **Value Format**      **Default** and Examples
============================== ===================== ===========================================================
 OpenTimeoutSecs                float                 **3600**  (wait for up to an hour)
 BeginStepPollingFrequencySecs  float                 **1**     (poll the file system with 1 second frequency
============================== ===================== ===========================================================

3. ``InSituAnalysis``. The engine is ``SST``. The parameters are set to:

============================== ===================== ===========================================================
 **Key**                       **Value Format**      **Default** and Examples
============================== ===================== ===========================================================
RendezvousReaderCount          integer               **1**      (Producer waits for the Consumer in Open)
QueueLimit                     integer               **1**      (only buffer one step)
QueueFullPolicy                string                **Block**  (wait for the Consumer to get every step)
FirstTimestepPrecious          bool                  false      (SST default)
AlwaysProvideLatestTimestep    bool                  false      (SST default)
============================== ===================== ===========================================================

4. ``InSituVisualization``. The engine is ``SST``. The parameters are set to:

============================== ===================== ===========================================================
 **Key**                       **Value Format**      **Default** and Examples
============================== ===================== ===========================================================
RendezvousReaderCount          integer               **0**       (Producer does NOT wait for Consumer in Open)
QueueLimit                     integer               **3**       (buffer first step + last two steps)
QueueFullPolicy                string                **Discard** (slow Consumer will miss out on steps)
FirstTimestepPrecious          bool                  **true**    (First step is kept around for late Consumers)
AlwaysProvideLatestTimestep    bool                  false       (SST default)
============================== ===================== ===========================================================


5. ``Code Coupling``. The engine is ``SST``. The parameters are set to:

============================== ===================== ===========================================================
 **Key**                       **Value Format**      **Default** and Examples
============================== ===================== ===========================================================
RendezvousReaderCount          integer               **1**      (Producer waits for the Consumer in Open)
QueueLimit                     integer               **1**      (only buffer one step)
QueueFullPolicy                string                **Block**  (wait for the Consumer to get every step)
FirstTimestepPrecious          bool                  false      (SST default)
AlwaysProvideLatestTimestep    bool                  false      (SST default)
============================== ===================== ===========================================================






