module Main where

import Lib
import Control.Monad
import Control.Applicative
import Options.Applicative
import Data.Semigroup ((<>))
import Data.List
import qualified Data.ByteString.Char8 as B
import Data.Maybe
import qualified Data.Map.Strict as M
import qualified Data.IntMap.Strict as I
import Data.Either (isRight, rights)
import qualified Conduit as C
import Data.Conduit.TQueue
import qualified Data.Conduit.Binary as CB
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Data.Conduit.Attoparsec as CA
import Control.Concurrent hiding (yield)
import Control.Concurrent.Async
import Control.Concurrent.STM
import Control.Concurrent.STM.TBMQueue
import Control.Concurrent.STM.TBMChan
import Control.Monad.Trans.Resource (register)


-- main
main :: IO ()
main = do
    let opts = info (helper <*> optargs)
            (fullDesc <> progDesc
                        "Trim PCR primer sequences from aligned reads (multi-threaded test version)"
                      <> header
                        "primerclip -- Swift Biosciences Accel-Amplicon targeted panel primer trimming tool v0.3.5")
    args <- execParser opts
    runPrimerTrimmingPEpar args -- 181205 test parallel execution performance
    {--
    runstats <- case (sereads args) of
                    True  -> runPrimerTrimmingSE args
                    False -> runPrimerTrimmingPE args
    --}
    putStrLn "primer trimming complete."
    -- writeRunStats (outfilename args) runstats -- 180226
-- end main

-- 180329 parse and trim as PairedAln sets
runPrimerTrimmingPE :: Opts -> IO RunStats
runPrimerTrimmingPE args = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- C.runConduitRes
              $ C.sourceFile (insamfile args)
              C..| CA.conduitParserEither parsePairedAlnsOrHdr
              C..| C.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              C..| C.concatC
              C..| C.mapC (trimprimerPairsE fmp rmp)
              C..| C.mapC flattenPairedAln
              C..| C.concatC
              C..| C.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              C..| C.getZipSink
                       (C.ZipSink (printAlnStreamToFile (outfilename args))
                                *> calcRunStats) -- 180226 --}
    return runstats

-- 181125 parse and trim single-end read alignments
runPrimerTrimmingSE :: Opts -> IO RunStats
runPrimerTrimmingSE args = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- C.runConduitRes
              $ C.sourceFile (insamfile args)
              C..| CA.conduitParserEither parseSingleAlnsOrHdr
              C..| C.mapC rightOrDefaultSingle -- convert parse fails to defaultAlignment
              C..| C.concatC
              C..| C.mapC (trimprimersE fmp rmp)
              C..| C.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              C..| C.getZipSink
                       (C.ZipSink (printAlnStreamToFile (outfilename args))
                                *> calcRunStats) -- 180226 --}
    return runstats

-- 181205 test parallel execution using a TBMQueue and async
runPrimerTrimmingPEpar :: Opts -> IO ()
runPrimerTrimmingPEpar args = do
    -- (fmp, rmp) <- createprimerbedmaps args
    q <- atomically $ newTBMQueue 10
    _ <- async $ C.runResourceT $ do
        _ <- register $ atomically $ closeTBMQueue q
        C.runConduitRes $ C.sourceFile (insamfile args) C..| sinkTBMQueue q
    _ <- replicateConcurrently_ 8 (runPrimerTrimPEchan args q) -- TODO: print run stats to log file
    return ()

-- runPrimerTrimPEchan :: Opts -> IO RunStats
runPrimerTrimPEchan args q = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- C.runConduitRes
              $ sourceTBMQueue q
              C..| CA.conduitParserEither parsePairedAlnsOrHdr
              C..| C.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              C..| C.concatC
              C..| C.mapC (trimprimerPairsE fmp rmp)
              C..| C.mapC flattenPairedAln
              C..| C.concatC
              C..| C.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              C..| C.getZipSink
                       (C.ZipSink (printAlnStreamToFile (outfilename args))
                                *> calcRunStats) -- 180226 --}
    return runstats
