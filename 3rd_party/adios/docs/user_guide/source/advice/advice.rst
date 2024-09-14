******
Advice
******

This list is similar to the Advice sections for each chapter in `The C++ Programming Language, Fourth Edition by Bjarne Stroustrup <http://www.stroustrup.com/4th.html>`_
The goal is to provide specific advice and good practices about the use of ADIOS2 in other projects. 

1. Use ``MPI_COMM_SELF`` to run MPI compiled versions of ADIOS 2 in "serial" mode

2. Use a runtime configuration file in the ``ADIOS`` constructor or ``adios2_init`` when targeting multiple engines

3. Check object validity when developing (similar to ``fstream``):

   -  C++: `operator bool`
         
         .. code-block:: c++ 
            
            if(var) engine.Put(var, data);
         
   -  C: NULL pointer 
         
         .. code-block:: c 
         
            if(var) adios2_put(engine, var, data, adios2_mode_deferred);
         
   -  Python: v2 ``__nonzero__`` v3 ``__bool__``. Note: do not check for None object
         .. code-block:: python
         
            if(var) engine.Put(var, data);
   
   -  Fortran: ``type%valid``
         
         .. code-block:: fortran
         
            if( adios%valid .eqv. .true. ) then
               adios2_declare_io(adios, io, "foo")
            end if
         
         
4. C++11 and Python: use ``try-catch`` (``try-except`` in Python) blocks to handle exceptions from ADIOS 2

5. C++11: use fixed-width types (``int8_t``, ``uint32_t``) for portability

6. Define your data structure: set of variables and attributes before developing. Data hierarchies/models can be built on top of ADIOS 2.

7. Read the documentation for :ref:`Supported Engines` before targeting development for a particular engine

8. MPI development: treat ``ADIOS`` constructor/destructor (``adios2_init``/``adios2_finalize``) and Engine ``Open`` and ``Close`` always as collective functions. For the most part, ADIOS 2 API functionality is local, but other Engine functions might follow other rules, :ref:`Supported Engines`.  

9. Use `Remove` functions carefully. They create dangling objects/pointers.

10. Thread-safety: treat ADIOS 2 as NOT thread-safe. Either use a mutex or only handle I/O from a master thread. ADIOS 2 is about performance, adding I/O serial algorithm operations into a parallel execution block may reduce parallel portions from Amdahl's Law. 

11. Prefer the high-level Python and C++ APIs for simple tasks that do not require performance. The more familiar Write/Read overloads for File I/O return native data constructs (``std::vector`` and ``numpy arrays``) immediately for a requested selection. ``open`` only explores the metadata index.

12. C++: prefer templates to ``void*`` to increase compile-time safety. Use ``IO::InquireVariableType("variableName")`` and ``adios2::GetType<T>()`` to cast upfront to a ``Variable<T>``. C++17 has ``std::any`` as an alternative to ``void*``. ADIOS 2 follows closely the STL model.

13. Understand ``Put`` and ``Get`` memory contracts from :ref:`Engine`

14. Prefer Put/Get ``Deferred`` mode, treat ``Sync`` as a special mode

15. ``Put Span``: create all spans in a step before populating them. Spans follow the same iterator invalidation rules as ``std::vector``, so use ``span.data()`` to always keep the span pointer up-to-date 

16. Always populate data before calling ``Put`` in deferred mode,
    and do not change it between ``Put`` and ``EndStep``, or ``Close``

17. Never call ``PerformPuts`` right before ``EndStep``.  This was a
    code pattern that had no adverse effects with the BP3/4 file
    engines and is present in some older code, but was never
    beneficial.
    
18. Use ``BeginStep`` and ``EndStep`` to write code that is portable
    across all ADIOS 2 Engine types: file and streaming.

19. Always use ``Close`` for every call to ``Open``.

20. C, Fortran: always call ``adios2_finalize`` for every call to ``adios2_init`` to avoid memory leaks.

21. Reminder: C++, C, Python: Row-Major, while Fortran: Column-Major. ADIOS 2 will handle interoperability between ordering. Remember that :ref:`bpls : Inspecting Data` is always a Row-Major reader so Fortran reader need to swap dimensions seen in bpls.  bpls: (slow, ...., fast) -> Fortran(fast,...,slow).

22. Fortran API: use the type members (``var%valid``, ``var%name``, etc.) to get extra type information.

23. Fortran C interoperability: Fortran bindings support the majority of applications using Fortran 90. We currently don't support the ``ISO_C_BINDING`` interoperability module in Fortran 2003. 

24. Always keep the ``IO`` object self-contained keeping its own set of ``Variables``, ``Attributes`` and ``Engines``. Do not combine Variables with multiple Engines or multiple modes, unless it's 100% guaranteed to be safe in your program avoiding Variable access conflicts.

25. Developers: explore the testing infrastructure ``ADIOS2/testing`` in ADIOS 2 as a starting point for using ADIOS 2 in your own testing environment. 

26. Become a super-user of :ref:`bpls : Inspecting Data` to analyze datasets generated by ADIOS 2.
 
