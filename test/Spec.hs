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
--import Data.Either (isRight, rights, fromRight, partitionEithers)
import Data.Either
import qualified Conduit as P
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
    -- determine input type (stdin or filename)
    let fnames = filenames args
        insource
            | null fnames          = P.stdinC
            | (head fnames) == "-" = P.stdinC
            | otherwise            = P.sourceFile insamfile
                where insamfile = head fnames -- should be safe
    (fmp, rmp) <- createprimerbedmaps args
    trimdalns <- P.runConduitRes
              $ insource
              P..| CA.conduitParserEither parseSAMtoPairedAlns -- parsePairedAlnsOrHdr
              P..| P.mapC rightOrDefaultPaird -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimerPairsE fmp rmp)
              P..| P.mapC flattenPairedAln
              P..| P.concatC
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.sinkList
    return trimdalns

-- 181125 parse and trim single-end read alignments
runPrimerTrimmingSEtest :: Opts -> IO [AlignedRead]
runPrimerTrimmingSEtest args = do
    let fnames = filenames args
        insource
            | null fnames          = P.stdinC
            | (head fnames) == "-" = P.stdinC
            | otherwise            = P.sourceFile insamfile
                where insamfile = head fnames -- should be safe
    (fmp, rmp) <- createprimerbedmaps args
    trimdalns <- P.runConduitRes
              $ insource
              P..| CA.conduitParserEither parseSingleAlnsOrHdr
              P..| P.mapC rightOrDefaultSingle -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimersE fmp rmp)
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.sinkList
    return trimdalns

-- 181230
-- readSAMTest :: Opts -> IO [AlignedRead]
readSAMTest args = do
    let fnames = filenames args
        insource
            | null fnames          = P.stdinC
            | (head fnames) == "-" = P.stdinC
            | otherwise            = P.sourceFile insamfile
                where insamfile = head fnames -- should be safe
    (fmp, rmp) <- createprimerbedmaps args
    alns <- P.runConduitRes
              $ insource
              P..| CA.conduitParserEither parseSAMtoPairedAlns -- parsePairedAlnsOrHdr
              P..| P.mapC rightOrDefaultPaird -- P..| P.concatC
              P..| P.sinkList
    return alns
