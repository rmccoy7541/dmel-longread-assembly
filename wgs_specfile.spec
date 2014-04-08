# assembly options

unitigger=bogart
merSize=31
merThreshold=auto*2
ovlMinLen=800
obtErrorRate=0.03
obtErrorLimit=4.5
ovlErrorRate=0.03
utgErrorRate=0.015
utgGraphErrorRate=0.015
utgGraphErrorLimit=0
utgMergeErrorRate=0.03
utgMergeErrorLimit=0

# memory management options

frgCorrBatchSize = 400000
frgCorrThreads   = 6
ovlCorrConcurrency = 6
ovlCorrBatchSize = 400000

merylMemory      = -segments 6 -threads 6
merylThreads     = 6

ovlThreads	  = 6
ovlHashBits      = 25
ovlHashBlockLength = 180000000

merOverlapperThreads         = 6
merOverlapperSeedBatchSize   = 400000
merOverlapperExtendBatchSize = 400000

cnsConcurrency   = 6

fakeUIDs = 1
