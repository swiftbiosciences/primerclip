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

-- main
main :: IO ()
main = do
    let opts = info (helper <*> optargs)
            (fullDesc <> progDesc
                        "Trim PCR primer sequences from aligned reads"
                      <> header
<<<<<<< HEAD
                        "primerclip -- Swift Biosciences Accel-Amplicon targeted panel primer trimming tool v0.3.3")
    args <- execParser opts
    runstats <- runPrimerTrimming args
=======
                        "primerclip -- Swift Biosciences Accel-Amplicon targeted panel primer trimming tool for single-end reads v0.3.5")
    args <- execParser opts
    runstats <- case (sereads args) of
                    True  -> runPrimerTrimmingSE args
                    False -> runPrimerTrimmingPE args
>>>>>>> 4ade68d... add single-end trim feature and command option flag to corrected CIGAR mods
    putStrLn "primer trimming complete."
    writeRunStats (outfilename args) runstats -- 180226
-- end main

-- 180329 parse and trim as PairedAln sets
runPrimerTrimmingPE :: Opts -> IO RunStats
runPrimerTrimmingPE args = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              P..| CA.conduitParserEither parsePairedAlnsOrHdr
              P..| P.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimerPairsE fmp rmp)
              P..| P.mapC flattenPairedAln
              P..| P.concatC
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.getZipSink
                       (P.ZipSink (printAlnStreamToFile (outfilename args))
                                *> calcRunStats) -- 180226 --}
    return runstats

<<<<<<< HEAD
=======
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
                                *> calcRunStats) -- 180226 --}
    return runstats
>>>>>>> 4ade68d... add single-end trim feature and command option flag to corrected CIGAR mods

