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
                        "primerclip -- Swift Biosciences Accel-Amplicon™ targeted panel primer trimming tool v0.3")
    args <- execParser opts
    runstats <- runPrimerTrimming args
    putStrLn "primer trimming complete."
    -- writeRunStats (outfilename args) runstats -- 180226
-- end main

-- {--
-- 180329 parse and trim as PairedAln sets
runPrimerTrimming :: Opts -> IO [PairedAln]
runPrimerTrimming args = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              -- P..| CB.lines
              P..| CA.conduitParserEither parsePairedAlnsOrHdr
              -- P..| P.mapC (A.parseOnly (hdralnparser <|> alnparser))
              P..| P.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              P..| P.mapCE (trimprimerPairsE fmp rmp)
              -- P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              {-- P..| P.getZipSink
                       (P.ZipSink (printAlnStreamToFile (outfilename args))
                    *> calcRunStats) -- 180226 --}
              -- P..| P.mapC (B.pack . show)
              P..| P.concatC
              P..| P.sinkList
    return runstats
--}

-- 180206 adapt to optional primer coords input file formats
runPrimerTrimming2 :: Opts -> IO RunStats
runPrimerTrimming2 args = do
    (fmp, rmp) <- createprimerbedmaps args
    runstats <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              P..| CB.lines
              P..| P.mapC (A.parseOnly (hdralnparser <|> alnparser))
              P..| P.mapC rightOrDefault -- convert parse fails to defaultAlignment
              P..| P.mapC (trimprimersE fmp rmp)
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.getZipSink
                       (P.ZipSink (printAlnStreamToFile (outfilename args))
                    *> calcRunStats) -- 180226
    return runstats

