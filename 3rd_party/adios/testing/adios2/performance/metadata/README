
We'd like to emphasize metadata processing in these tests and hopefully isolate the following costs -

Writer side:
	metadata composition	- anything done to combine metadata rank contributions

Reader side:
	metadata installation	- creation of variables, blocks info 
	metadata traversal	- generating read requests, etc.


plan:

(Likely most effective for SST, but hopefully useful for other engines.)

For measuring writer-side metadata composition.
  - in sst (other staging engines?) measure writer completion time with zero readers (metadata composition happens, but data is never delivered anywhere).
  - for file engines, just make data small to de-emphasize it and measure write time (to /dev/null ?)


for measuring reader-side metadata installation time.
  -  in sst, open both sides of the stream, but delay reader some amount of time so that all metadata arrives, then measure time for BeginStep/EndStep without any Get()s.
     Must turn off preload to avoid data interference, and delay should be sufficiently long (whatever that means)
  -  in file engine.  All parsing/installation *might* happen in Open(), so we should include that?  Obviously without a delay, because the file is already written


for measuring reader-side traversal
  - in sst, same as above, but add Get()s and disable actual reads of data (obviously don't check for data correctness)

