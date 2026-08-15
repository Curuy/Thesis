178,0
S,SleepMS,Sleeps for n milliseconds,0,1,0,0,1,0,0,0,0,-1,,-38,-38,-38,-38,-38,-38
S,StartDistributedWorkers,"Start N DistributedWorker instances via system call, where the worker specifications appear in the file `file`",0,2,0,0,1,0,0,0,0,-1,,0,0,-1,,-38,-38,-38,-38,-38,-38
S,ParallelMap,"Return the Sequence [fun(x) : x in L] using a parallel implementation. Note that the intrinsic `fun` must be available at magma startup, and cannot refer to any variables assigned in the current session. Additionally, ParallelTools must be Attached on startup",0,2,0,0,0,0,0,0,0,82,,0,0,-1,,82,-38,-38,-38,-38,-38
S,testSleep,,0,1,0,0,0,0,0,0,0,148,,148,-38,-38,-38,-38,-38
S,KillProcessOnPort,Kills all processes listening on the given port,0,1,0,0,1,0,0,0,0,-1,,-38,-38,-38,-38,-38,-38
S,IsWorkerProcess,Determines if this file has been launched as a worker,0,0,0,0,0,0,0,36,-38,-38,-38,-38,-38
S,DefaultParallelSetup,Return the default host/port/socket information,0,0,0,0,0,0,0,298,148,-38,-38,-38,-38
S,DoWorkThenDie,Launch a worker instance of the function `fun`. Connects to a currently open manager process. Then terminate magma,0,1,0,0,1,0,0,0,0,-1,,-38,-38,-38,-38,-38,-38
S,ActivateWorkers,"Runs the code in `file` in parallel, with workers receiving elements of `tasks` as input",0,2,0,0,0,0,0,0,0,-1,,0,0,-1,,-1,-38,-38,-38,-38,-38
