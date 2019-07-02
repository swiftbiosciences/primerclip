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
                        "primerclip -- Swift Biosciences Accel-Amplicon targeted panel primer trimming tool v0.3.9")
    args <- execParser opts
    runstats <- selectRunmode args
    {--
    runstats <- case (sereads args) of
                    True  -> runPrimerTrimmingSE args
                    False -> runPrimerTrimmingPE args
    --}
    putStrLn "primer trimming complete."
    writeRunStats (outfilename args) runstats -- 180226
-- end main

-- 190702
selectRunmode :: Opts -> IO RunStats
selectRunmode args
    | (not isSE) && (not isFQ) = runPrimerTrimmingPE args
    | (not isSE) && isFQ = runPrimerTrimmingPE_Fastq args
    | isSE && (not isFQ) = runPrimerTrimmingSE args
    | isSE && isFQ = runPrimerTrimmingSE_Fastq args
    | otherwise = error "[!] Invalid command line arguments: please check your command"
        where isSE = sereads args
              isFQ = fqout args

-- 190627 TODO: implement runPrimerTrimmingPE which outputs R1 and R2 FASTQ
-- files (properly handle reads which are no longer paired)

-- 180329 parse and trim as PairedAln sets
runPrimerTrimmingPE :: Opts -> IO RunStats
runPrimerTrimmingPE args = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              P..| CA.conduitParserEither parseSAMtoPairedAlns -- parsePairedAlnsOrHdr
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

-- 190701 output as primer-trimmed FASTQ files
runPrimerTrimmingPE_Fastq :: Opts -> IO RunStats
runPrimerTrimmingPE_Fastq args = do
    (fmp, rmp) <- createprimerbedmaps args
    let outfnameprefix = reverse $ drop 4 $ reverse $ outfilename args
    runstats <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              P..| CA.conduitParserEither parseSAMtoPairedAlns -- parsePairedAlnsOrHdr
              P..| P.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimerPairsE fmp rmp)
              P..| P.mapC flattenPairedAlnFastq
              P..| P.concatC
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.getZipSink
                       (P.ZipSink
                            (printAlnStreamToFastqs outfnameprefix)
                                *> calcRunStats)
    return runstats

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

-- 190702 print R1 FASTQ file as output
runPrimerTrimmingSE_Fastq :: Opts -> IO RunStats
runPrimerTrimmingSE_Fastq args = do
    (fmp, rmp) <- createprimerbedmaps args
    let outfnameprefix = reverse $ drop 4 $ reverse $ outfilename args
    runstats <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              P..| CA.conduitParserEither parseSingleAlnsOrHdr
              P..| P.mapC rightOrDefaultSingle -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimersE fmp rmp)
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.getZipSink
                       (P.ZipSink
                            (printAlnStreamToFastq1 outfnameprefix)
                                *> calcRunStats) -- 180226 --}
    return runstats
