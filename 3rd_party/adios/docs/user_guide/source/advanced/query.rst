#################
ADIOS2 query API
#################

The query API in ADIOS2 allows a client to pass a query in XML or json format,
and get back a list of blocks or sub-blocks that contains hits.
Both BP4 and BP5 engines are supported.  


The interface
=============
User is expected to pass a query file (configFile), and init a read engine (engine)
to construct a query and evaluate using the engine.
(note that the engine and query should be using the same ADIOS IO)

.. code-block:: c++
		
    class QueryWorker
    {
    public:
        // configFile has query, can be either xml or json
        QueryWorker(const std::string &configFile, adios2::Engine &engine);

	     // touched_blocks is a list of regions specified by (start, count),
	     // that contains data that satisfies the query file
        void GetResultCoverage(std::vector<adios2::Box<adios2::Dims>> &touched_blocks);
    ... 
    }

A Sample Compound Query  
-----------------------

This query targets a 1D variable "doubleV", data of interest is (x  > 6.6) or (x < -0.17) or (2.8 < x < 2.9) 
In addition, this query also specied an output region [start=5,count=80]. 


.. code-block:: xml

	<adios-query>
  	  <io name="query">
   	  <var name="doubleV">
      	    <boundingbox  start="5" count="80"/>
            <op value="OR">
              <range  compare="GT" value="6.6"/>
              <range  compare="LT" value="-0.17"/>
              <op value="AND">
                 <range  compare="LT" value="2.9"/>
                 <range  compare="GT" value="2.8"/>
              </op>
            </op>
          </var>
          </io>
        </adios-query>
		

Code EXAMPLES:
==============
C++:
----
.. code-block:: c++
		
    while (reader.BeginStep() == adios2::StepStatus::OK)
    {
        adios2::QueryWorker w = adios2::QueryWorker(queryFile, reader);
        w.GetResultCoverage(touched_blocks);
	
        std::cout << " ... now can read out touched blocks ... size=" << touched_blocks.size()
                  << std::endl;
    }


The Full C++ example is here:
    https://github.com/ornladios/ADIOS2/blob/master/examples/query/test.cpp
    

Python:
-------

.. code-block:: python
	
	while (reader.BeginStep() == adios2.StepStatus.OK):
        # say only rank 0 wants to process result
        var = [queryIO.InquireVariable("T")]

        if (rank == 0):
            touched_blocks = w.GetResult()
            doAnalysis(reader, touched_blocks, var)

Full python example is here:
	https://github.com/ornladios/ADIOS2/blob/master/testing/adios2/bindings/python/TestQuery.py

	This example generates data, the query file (in xml) and runs the query, all in python.

