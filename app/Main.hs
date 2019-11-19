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
import qualified Conduit as P
import qualified Data.Conduit.Binary as CB
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Data.Conduit.Attoparsec as CA
import System.IO (stdin, stdout, stderr)

-- main
main :: IO ()
main = do
    let opts = info (helper <*> optargs)
            (fullDesc <> progDesc
                        "Trim PCR primer sequences from aligned reads"
                      <> header
                        "primerclip -- Swift Biosciences Accel-Amplicon targeted panel primer trimming tool v0.3.8")
    args <- execParser opts
    -- {--
    runstats <- case (sereads args) of
                    True  -> runPrimerTrimmingSE args
                    False -> runPrimerTrimmingPE args
    --}
    -- runstats <- selectRunmode args
    -- B.hPutStrLn stderr "primer trimming complete."
    writeRunStats (outfilename args) runstats -- 180226
-- end main

{--
-- 191119 allow for optional use of stdin without breaking existing command
-- line conventions for this tool
selectRunmode :: Opts -> IO RunStats
selectRunmode args
    | isSingleEnd && isStdIn     = hrunPrimerTrimmingSE args
    | isSingleEnd && isFilenames = runPrimerTrimmingSE args
    | isPairedEnd && isStdIn     = hrunPrimerTrimmingPE args
    | isPairedEnd && isFilenames = runPrimerTrimmingPE args
    | otherwise                  = runPrimerTrimmingPE args
        where isSingleEnd = sereads args
              isPairedEnd = not isSingleEnd
              isStdIn = case (insamfile args) of
                "-" -> True
                _   -> False
              isFilenames = not isStdIn

-- 191119 stdin version
hrunPrimerTrimmingPE :: Opts -> IO RunStats
hrunPrimerTrimmingPE args = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- P.runConduitRes
              $ P.sourceIOHandle (return stdin)
              P..| CA.conduitParserEither parseSAMtoPairedAlns -- parsePairedAlnsOrHdr
              P..| P.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimerPairsE fmp rmp)
              P..| P.mapC flattenPairedAln
              P..| P.concatC
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.getZipSink
                       (P.ZipSink hprintAlnStreamToFile *> calcRunStats) -- 191119
    return runstats
--}

-- 180329 parse and trim as PairedAln sets
runPrimerTrimmingPE :: Opts -> IO RunStats
runPrimerTrimmingPE args = do
    (fmp, rmp) <- createprimerbedmaps args
    let readFunc = case (insamfile args) of
                        "-"   -> P.sourceIOHandle (return stdin)
                        fname -> P.sourceFile fname
        writeFunc = case (outfilename args) of
                        "-"   -> hprintAlnStreamToFile
                        fname -> printAlnStreamToFile fname
    runstats <- P.runConduitRes
              $ readFunc -- P.sourceFile (insamfile args)
              P..| CA.conduitParserEither parseSAMtoPairedAlns -- parsePairedAlnsOrHdr
              P..| P.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimerPairsE fmp rmp)
              P..| P.mapC flattenPairedAln
              P..| P.concatC
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.getZipSink
                       (P.ZipSink writeFunc *> calcRunStats)
    return runstats

-- 191119 std and stdout
runPrimerTrimmingSE :: Opts -> IO RunStats
runPrimerTrimmingSE args = do
    (fmp, rmp) <- createprimerbedmaps args
    let readFunc = case (insamfile args) of
                        "-"   -> P.sourceIOHandle (return stdin)
                        fname -> P.sourceFile fname
        writeFunc = case (outfilename args) of
                        "-"   -> hprintAlnStreamToFile
                        fname -> printAlnStreamToFile fname
    runstats <- P.runConduitRes
              $ readFunc -- $ P.sourceIOHandle (return stdin)
              P..| CA.conduitParserEither parseSingleAlnsOrHdr
              P..| P.mapC rightOrDefaultSingle -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimersE fmp rmp)
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.getZipSink
                       (P.ZipSink writeFunc *> calcRunStats) -- 191119
    return runstats

{--
-- 181125 parse and trim single-end read alignments
runPrimerTrimmingSE :: Opts -> IO RunStats
runPrimerTrimmingSE args = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              P..| CA.conduitParserEither parseSingleAlnsOrHdr
              P..| P.mapC rightOrDefaultSingle -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimersE fmp rmp)
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.getZipSink
                       (P.ZipSink (printAlnStreamToFile (outfilename args))
                                *> calcRunStats) -- 180226
    return runstats
--}
