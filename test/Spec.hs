module Spec where

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
import qualified Data.Conduit.Binary as CB
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Data.Conduit.Attoparsec as CA

main :: IO ()
main = do
    putStrLn "test suite under construction."


-- test functions for running in main

-- 180329 parse and trim as PairedAln sets
runPrimerTrimmingPETest :: Opts -> IO [AlignedRead]
runPrimerTrimmingPETest args = do
    (fmp, rmp) <- createprimerbedmaps args
    trimdalns <- C.runConduitRes
              $ C.sourceFile (insamfile args)
              C..| CA.conduitParserEither parsePairedAlnsOrHdr
              C..| C.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              C..| C.concatC
              C..| C.mapC (trimprimerPairsE fmp rmp)
              C..| C.mapC flattenPairedAln
              C..| C.concatC
              C..| C.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              C..| C.sinkList
    return trimdalns

-- 181125 parse and trim single-end read alignments
runPrimerTrimmingSEtest :: Opts -> IO [AlignedRead]
runPrimerTrimmingSEtest args = do
    (fmp, rmp) <- createprimerbedmaps args
    trimdalns <- C.runConduitRes
              $ C.sourceFile (insamfile args)
              C..| CA.conduitParserEither parseSingleAlnsOrHdr
              C..| C.mapC rightOrDefaultSingle -- convert parse fails to defaultAlignment
              C..| C.concatC
              C..| C.mapC (trimprimersE fmp rmp)
              C..| C.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              C..| C.sinkList
    return trimdalns
