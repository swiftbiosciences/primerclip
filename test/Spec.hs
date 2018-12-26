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
import qualified Conduit as P
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

main :: IO ()
main = do
    putStrLn "test suite under construction."


-- test functions for running in main

-- 180329 parse and trim as PairedAln sets
runPrimerTrimmingPETest :: Opts -> IO [AlignedRead]
runPrimerTrimmingPETest args = do
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
    return trimdalns

-- 181125 parse and trim single-end read alignments
runPrimerTrimmingSEtest :: Opts -> IO [AlignedRead]
runPrimerTrimmingSEtest args = do
    (fmp, rmp) <- createprimerbedmaps args
    trimdalns <- P.runConduitRes
              $ P.sourceFile (insamfile args)
              P..| CA.conduitParserEither parseSingleAlnsOrHdr
              P..| P.mapC rightOrDefaultSingle -- convert parse fails to defaultAlignment
              P..| P.concatC
              P..| P.mapC (trimprimersE fmp rmp)
              P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
              P..| P.sinkList
    return trimdalns

-- 181205 test parallel execution using a TBMQueue and async
-- runPrimerTrimmingPEpar :: Opts -> IO ()
runPrimerTrimmingSEpar args = do
    (fmp, rmp) <- createprimerbedmaps args
    q <- atomically $ newTBMQueue 20
    _ <- async $ P.runResourceT $ do
        _ <- register $ atomically $ closeTBMQueue q
        P.runConduitRes $ P.sourceFile (insamfile args) P..| sinkTBMQueue q
    r <- replicateConcurrently 16 ( P.runConduitRes
                                  $ sourceTBMQueue q
                               P..| CA.conduitParserEither parseSingleAlnsOrHdr
                               P..| P.mapC rightOrDefaultSingle
                               P..| P.concatC
                               {--
                               P..| P.mapC (trimprimerPairsE fmp rmp)
                               P..| P.mapC flattenPairedAln
                               P..| P.concatC
                               P..| P.filterC (\x -> (qname x) /= "NONE") -- remove dummy alignments
                               --}
                               P..| P.mapC printAlignmentOrHdr
                               P..| P.unlinesAsciiC
                               P..| P.stdoutC )
    return r
