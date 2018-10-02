main :: IO ()
main = do
    putStrLn "test suite under construction."


-- test functions for running in main

-- 180329 parse and trim as PairedAln sets
runPrimerTrimmingTest :: Opts -> IO [AlignedRead]
runPrimerTrimmingTest args = do
    (fmp, rmp) <- createprimerbedmaps args
    trimdalns <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              P..| CA.conduitParserEither parsePairedAlnsOrHdr
              P..| P.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimerPairsE fmp rmp)
              P..| P.mapC flattenPairedAln
              P..| P.concatC
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.sinkList
