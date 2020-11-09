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
import Test.Hspec

main :: IO ()
main = hspec $ do
    describe "Lib" $ do
        describe "uchrparser" $ do
            context "with valid input" $ do
                it "can parse chrN chromosome names" $ do
                    A.parseOnly uchrparser "chr1"
                    `shouldBe`
                    (Right $ ChrAlt "chr1" :: Either String UChr)
                it "can parse N chromosome names" $ do
                    A.parseOnly uchrparser "1"
                    `shouldBe`
                    (Right $ ChrAlt "1" :: Either String UChr)
                it "can parse non-humna FASTA sequence names" $ do
                    A.parseOnly uchrparser "NC_045512.2 DON'T_PARSE_ME"
                    `shouldBe`
                    (Right $ ChrAlt "NC_045512.2" :: Either String UChr)
            {--
            context "with invalid input" $ do
                it "returns parse error" $ do
                    A.parseOnly uchrparser ""
                    `shouldBe`
                    (Left "parse error" :: Either String UChr)
            --}
        describe "bedrecparser" $ do
            context "with valid input" $ do
                it "can parse a valid BED record" $ do
                    A.parseOnly
                        bedrecparser
                        "NC_045512.2\t791\t894\tcovid19genome_200-29703_s7490_U_81"
                        `shouldBe`
                        ( Right
                        $ BedRecord  (ChrAlt "NC_045512.2")
                                      791
                                      894
                                      "covid19genome_200-29703_s7490_U_81"
                          :: Either String BedRecord)
        describe "bedPEparser" $ do
            context "with valid input" $ do
                it "can parse a valid BED record" $ do
                    A.parseOnly
                        bedPEparser
                        "NC_045512.2\t772\t791\tNC_045512.2\t894\t915\tcovid19genome_200-29703_s7490_U_81"
                        `shouldBe`
                        ( Right
                        $ BEDPE  (ChrAlt "NC_045512.2")
                                  772
                                  791
                                 (ChrAlt "NC_045512.2")
                                  894
                                  915
                                  "covid19genome_200-29703_s7490_U_81"
                          :: Either String BEDPE)
        describe "masterparser" $ do
            context "with valid input" $ do
                it "can parse a valid masterfile record" $ do
                    A.parseOnly
                        masterparser
                        "NC_045512.2\t791\t894\tcovid19genome_200-29703_s7490_U_81\t772\t791\tcovid19genome_200-29703_s7490_U_81F\t894\t915\tcovid19genome_200-29703_s7490_U_81R\tACCCGTGAACTCATGCGTG\tTCGGACAAAGTGCATGAAGCT"
                        `shouldBe`
                        ( Right
                        $ MasterRecord  (ChrAlt "NC_045512.2")
                                        791
                                        894
                                        "covid19genome_200-29703_s7490_U_81"
                                        772
                                        791
                                        894
                                        915
                                        "covid19genome_200-29703_s7490_U_81F"
                                        "covid19genome_200-29703_s7490_U_81R"
                                        "ACCCGTGAACTCATGCGTG"
                                        "TCGGACAAAGTGCATGAAGCT"
                          :: Either String MasterRecord)


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
