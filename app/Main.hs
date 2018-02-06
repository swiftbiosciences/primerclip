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

-- main
main :: IO ()
main = do
    let opts = info (helper <*> optargs)
            (fullDesc <> progDesc
                        "Trim PCR primer sequences from aligned reads"
                      <> header
                        "primerclip -- Swift Biosciences primer trimming tool v0.2")
    args <- execParser opts
    trimstats <- runPrimerTrimming args
    putStrLn $ (show trimstats) ++ " alignments with >=1 primer bases trimmed."
    putStrLn "primers trimming complete."
-- end main

-- record to store command line arguments
data Cmd = Cmd { samfile :: String
               , masterfile :: String
               , outsamfile :: String
               } deriving (Show, Eq)

-- parse command line arguments
cargs :: Parser Cmd
cargs = Cmd
    <$> argument str (metavar "SAM_INFILE")
    <*> argument str (metavar "PANEL_MASTER_INFILE")
    <*> argument str (metavar "OUTPUT_SAM_FILENAME")

data CoordFmt = MASTER | BEDPE deriving (Show, Eq, Read)

-- trimming function (returns number of alignments clipped by >=1 base [DEBUG])
runPrimerTrimming' :: Cmd -> IO Integer
runPrimerTrimming' args = do
    m <- getMasterFile $ masterfile args
    let fmp = makechrbedmap $ masterToFPrimerBED m
        rmp = makechrbedmap $ masterToRPrimerBED m
    trimstats <- P.runConduitRes
               $ P.sourceFile (samfile args)
               P..| CB.lines
               P..| P.mapC (A.parseOnly (hdralnparser <|> alnparser))
               P..| P.mapC rightOrDefault -- convert parse fails to defaultAlignment
               P..| P.mapC (trimprimersE fmp rmp)
               P..| P.filterC checknonzeroCigMatch
               P..| P.mapC checkCigarSeqlen
               P..| P.filterC (\x -> (qname x) /= "NONE") -- remove malformed alignments
               P..| P.getZipSink
                        (P.ZipSink (printAlnStreamToFile (outsamfile args))
                     *>  P.ZipSink calculateTrimStats)
    return trimstats

-- 180206 adapt to optional primer coords input file formats
runPrimerTrimming :: Opts -> IO Integer
runPrimerTrimming args = do
    (fmp, rmp) <- createprimerbedmaps args
    -- m <- getMasterFile $ masterfile args
    -- let fmp = makechrbedmap $ masterToFPrimerBED m
    --     rmp = makechrbedmap $ masterToRPrimerBED m
    trimstats <- P.runConduitRes
               $ P.sourceFile (insamfile args)
               P..| CB.lines
               P..| P.mapC (A.parseOnly (hdralnparser <|> alnparser))
               P..| P.mapC rightOrDefault -- convert parse fails to defaultAlignment
               P..| P.mapC (trimprimersE fmp rmp)
               P..| P.filterC checknonzeroCigMatch
               P..| P.mapC checkCigarSeqlen
               P..| P.filterC (\x -> (qname x) /= "NONE") -- remove malformed alignments
               P..| P.getZipSink
                        (P.ZipSink (printAlnStreamToFile (outfilename args))
                     *>  P.ZipSink calculateTrimStats)
    return trimstats

