module Lib
    ( module Lib
    ) where

import System.Process
import Control.Applicative
import Options.Applicative hiding (flag)
import Data.Semigroup ((<>))
import Control.Monad
import Control.Monad.IO.Class (liftIO)
import Control.Exception (evaluate)
import qualified Data.Either as U
import Data.Ord (comparing)
import qualified Data.Char as C
import Data.Foldable (toList)
import qualified Data.Map.Strict as M
import Data.List
import qualified Data.ByteString.Char8 as B
import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Conduit as P
import qualified Data.Vector as V
import Data.Bits
import Data.Digits
import Data.Ord (comparing, Down)
import Data.Monoid hiding ((<>))
import qualified Data.IntMap.Strict as I
import Data.Maybe
import qualified Data.Set as S
import GHC.Generics (Generic)

{--
    2016-10-14 Jonathan Irish
    Post-alignment primer trimming tool v0.2

    NOTE: reverted to last commit before CIGAR trimming code was updated, which
    resulted in some edge case errors which require debugging.

--}

------------------------------------------------------------------------------
--------------------------------  Data Types  --------------------------------
------------------------------------------------------------------------------

-- alignment record containing alignment data for each read
data AlignedRead = AlignedRead { qname :: B.ByteString
                               , flag :: Int
                               , rname :: UChr
                               , pos :: Integer
                               , endpos :: Integer
                               , mapqual :: Integer
                               , cigar :: B.ByteString
                               , trimdcigar :: B.ByteString
                               , cigmap :: CigarMap
                               , trimdcigmap :: CigarMap
                               , rnext :: B.ByteString
                               , pnext :: Integer
                               , tlen :: Integer
                               , refseq :: B.ByteString
                               , basequal :: B.ByteString
                               , optfields :: B.ByteString
                               , strand :: B.ByteString
                               , paired :: Bool
                               , mapped :: Bool
                               , fint :: [BedRecord]
                               , rint :: [BedRecord]
                               , pintflag :: Bool
                               , trimdflag :: Bool
                               , trimdpos :: Integer
                               , trimdendpos :: Integer
                               , headerstrings :: Header
                               , isheader :: Bool
                               , mid :: B.ByteString
                               , xmid :: B.ByteString
                               , tbed :: BedRecord
                               } deriving (Eq, Show, Generic)

instance Ord AlignedRead where
    compare =  comparing rname
            <> comparing pos
            <> comparing mapqual

defaultAlignment = AlignedRead { qname = "NONE"
                               , flag = 0
                               , rname = NONE
                               , pos = 0
                               , endpos = 0
                               , mapqual = 0
                               , cigar = ""
                               , trimdcigar = "*"
                               , cigmap = [(0,"*")]
                               , trimdcigmap = [(0,"*")]
                               , rnext = "*"
                               , pnext = 0
                               , tlen = 0
                               , refseq = ""
                               , basequal = ""
                               , optfields = ""
                               , strand = ""
                               , paired = False
                               , mapped = False
                               , fint = []
                               , rint = []
                               , pintflag = False
                               , trimdflag = False
                               , trimdpos = 0
                               , trimdendpos = 0
                               , headerstrings = []
                               , isheader = False
                               , mid = ""
                               , xmid = ""
                               , tbed = defaultBed
                               }

-- 180206
data BEDPE = BEDPE { chr1 :: UChr
                   , start1 :: Integer
                   , end1 :: Integer
                   , chr2 :: UChr
                   , start2 :: Integer
                   , end2 :: Integer
                   , bedpename :: B.ByteString
                   } deriving (Show, Eq, Generic)

instance Ord BEDPE where
    compare =  comparing chr1
            <> comparing start1

defaultBEDPE = BEDPE NONE 0 0 NONE 0 0 "DEFAULTBEDPE"


data BedRecord = BedRecord { bedchr :: UChr
                           , bedstart :: Integer
                           , bedend :: Integer
                           , bedname :: B.ByteString
                           } deriving (Show, Eq, Generic)

instance Ord BedRecord where
    compare =  comparing bedchr
            <> comparing bedstart

defaultBed = BedRecord { bedchr = NONE
                       , bedstart = 0
                       , bedend = 0
                       , bedname = "DEFAULT"
                       }

-- 170926 use master file as source of primer and target coords
data MasterRecord = MasterRecord { mchrom :: UChr
                                 , mtargetstart :: Integer
                                 , mtargetend :: Integer
                                 , mampname :: B.ByteString
                                 , mforstart :: Integer
                                 , mforend :: Integer
                                 , mrevstart :: Integer
                                 , mrevend :: Integer
                                 , mforname :: B.ByteString
                                 , mrevname :: B.ByteString
                                 , mforseq :: B.ByteString
                                 , mrevseq :: B.ByteString
                                 } deriving (Show, Eq)

instance Ord MasterRecord where
    compare = comparing mchrom
           <> comparing mtargetstart
           <> comparing mtargetend

defaultMasterRec = MasterRecord { mchrom = NONE
                                , mtargetstart = 0
                                , mtargetend = 0
                                , mampname = "EMPTY"
                                , mforstart = 0
                                , mforend = 0
                                , mrevstart = 0
                                , mrevend = 0
                                , mforname = "EMPTY"
                                , mrevname = "EMPTY"
                                , mforseq = "EMPTY"
                                , mrevseq = "EMPTY"
                                }

-- 170919 add MidFamily record type to group [AlignedRead] by mID
data MidFamily = MidFamily { chrom :: UChr
                           , alnstart :: Integer
                           , alnend :: Integer
                           , midseq :: B.ByteString
                           , alns :: [AlignedRead]
                           , alncount :: Integer
                           , txposedcount :: Integer
                           , postxposdcnt :: Integer
                           , negtxposdcnt :: Integer
                           , targetbed :: BedRecord
                           } deriving (Show, Eq)

instance Ord MidFamily where
    compare = comparing chrom
           <> comparing alnstart
           <> comparing alnend
           <> comparing alncount

defaultMidFam = MidFamily { chrom = NONE
                          , alnstart = 0
                          , alnend = 0
                          , midseq = "EMPTY"
                          , alns = []
                          , alncount = 0
                          , txposedcount = 0
                          , postxposdcnt = 0
                          , negtxposdcnt = 0
                          , targetbed = defaultBed
                          }

type Header = [B.ByteString]
type Alignments = [AlignedRead]
type BED = V.Vector BedRecord
type SAM = (Header, Alignments)
type CMap = M.Map UChr BedMap
type BedMap = I.IntMap BedRecord
type CigarMap = [(Integer, B.ByteString)]


-- Internal sum type for all chromosome names in both UCSC and ensembl naming
-- schemes.
-- NOTE: only GRCh37 non-canonical chromosomes recognized this version
data UChr = C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 | C13
          | C14 | C15 | C16 | C17 | C18 | C19 | C20 | C21 | C22 | CX | CY | CMT
          | Chr1 | Chr2 | Chr3 | Chr4 | Chr5 | Chr6 | Chr7 | Chr8 | Chr9
          | Chr10 | Chr11 | Chr12 | Chr13 | Chr14 | Chr15 | Chr16 | Chr17
          | Chr18 | Chr19 | Chr20 | Chr21 | Chr22 | ChrX | ChrY | ChrM
          | GL000207P1
          | GL000226P1
          | GL000229P1
          | GL000231P1
          | GL000210P1
          | GL000239P1
          | GL000235P1
          | GL000201P1
          | GL000247P1
          | GL000245P1
          | GL000197P1
          | GL000203P1
          | GL000246P1
          | GL000249P1
          | GL000196P1
          | GL000248P1
          | GL000244P1
          | GL000238P1
          | GL000202P1
          | GL000234P1
          | GL000232P1
          | GL000206P1
          | GL000240P1
          | GL000236P1
          | GL000241P1
          | GL000243P1
          | GL000242P1
          | GL000230P1
          | GL000237P1
          | GL000233P1
          | GL000204P1
          | GL000198P1
          | GL000208P1
          | GL000191P1
          | GL000227P1
          | GL000228P1
          | GL000214P1
          | GL000221P1
          | GL000209P1
          | GL000218P1
          | GL000220P1
          | GL000213P1
          | GL000211P1
          | GL000199P1
          | GL000217P1
          | GL000216P1
          | GL000215P1
          | GL000205P1
          | GL000219P1
          | GL000224P1
          | GL000223P1
          | GL000195P1
          | GL000212P1
          | GL000222P1
          | GL000200P1
          | GL000193P1
          | GL000194P1
          | GL000225P1
          | GL000192P1
          | NC_007605
          | NONE
            deriving (Eq, Ord, Bounded, Enum, Generic)

-- 170508 type for tracking chromosome name type
data ChromNameFormat = GRC | UCSC
                        deriving (Show, Eq, Ord, Bounded, Enum, Read)

instance Show UChr where
    show C1 = "1"
    show C2 = "2"
    show C3 = "3"
    show C4 = "4"
    show C5 = "5"
    show C6 = "6"
    show C7 = "7"
    show C8 = "8"
    show C9 = "9"
    show C10 = "10"
    show C11 = "11"
    show C12 = "12"
    show C13 = "13"
    show C14 = "14"
    show C15 = "15"
    show C16 = "16"
    show C17 = "17"
    show C18 = "18"
    show C19 = "19"
    show C20 = "20"
    show C21 = "21"
    show C22 = "22"
    show CX = "X"
    show CY = "Y"
    show CMT = "MT"
    show GL000207P1 = "GL000207.1"
    show GL000226P1 = "GL000226.1"
    show GL000229P1 = "GL000229.1"
    show GL000231P1 = "GL000231.1"
    show GL000210P1 = "GL000210.1"
    show GL000239P1 = "GL000239.1"
    show GL000235P1 = "GL000235.1"
    show GL000201P1 = "GL000201.1"
    show GL000247P1 = "GL000247.1"
    show GL000245P1 = "GL000245.1"
    show GL000197P1 = "GL000197.1"
    show GL000203P1 = "GL000203.1"
    show GL000246P1 = "GL000246.1"
    show GL000249P1 = "GL000249.1"
    show GL000196P1 = "GL000196.1"
    show GL000248P1 = "GL000248.1"
    show GL000244P1 = "GL000244.1"
    show GL000238P1 = "GL000238.1"
    show GL000202P1 = "GL000202.1"
    show GL000234P1 = "GL000234.1"
    show GL000232P1 = "GL000232.1"
    show GL000206P1 = "GL000206.1"
    show GL000240P1 = "GL000240.1"
    show GL000236P1 = "GL000236.1"
    show GL000241P1 = "GL000241.1"
    show GL000243P1 = "GL000243.1"
    show GL000242P1 = "GL000242.1"
    show GL000230P1 = "GL000230.1"
    show GL000237P1 = "GL000237.1"
    show GL000233P1 = "GL000233.1"
    show GL000204P1 = "GL000204.1"
    show GL000198P1 = "GL000198.1"
    show GL000208P1 = "GL000208.1"
    show GL000191P1 = "GL000191.1"
    show GL000227P1 = "GL000227.1"
    show GL000228P1 = "GL000228.1"
    show GL000214P1 = "GL000214.1"
    show GL000221P1 = "GL000221.1"
    show GL000209P1 = "GL000209.1"
    show GL000218P1 = "GL000218.1"
    show GL000220P1 = "GL000220.1"
    show GL000213P1 = "GL000213.1"
    show GL000211P1 = "GL000211.1"
    show GL000199P1 = "GL000199.1"
    show GL000217P1 = "GL000217.1"
    show GL000216P1 = "GL000216.1"
    show GL000215P1 = "GL000215.1"
    show GL000205P1 = "GL000205.1"
    show GL000219P1 = "GL000219.1"
    show GL000224P1 = "GL000224.1"
    show GL000223P1 = "GL000223.1"
    show GL000195P1 = "GL000195.1"
    show GL000212P1 = "GL000212.1"
    show GL000222P1 = "GL000222.1"
    show GL000200P1 = "GL000200.1"
    show GL000193P1 = "GL000193.1"
    show GL000194P1 = "GL000194.1"
    show GL000225P1 = "GL000225.1"
    show GL000192P1 = "GL000192.1"
    show NC_007605 = "NC_007605"
    show Chr1 = "chr1"
    show Chr2 = "chr2"
    show Chr3 = "chr3"
    show Chr4 = "chr4"
    show Chr5 = "chr5"
    show Chr6 = "chr6"
    show Chr7 = "chr7"
    show Chr8 = "chr8"
    show Chr9 = "chr9"
    show Chr10 = "chr10"
    show Chr11 = "chr11"
    show Chr12 = "chr12"
    show Chr13 = "chr13"
    show Chr14 = "chr14"
    show Chr15 = "chr15"
    show Chr16 = "chr16"
    show Chr17 = "chr17"
    show Chr18 = "chr18"
    show Chr19 = "chr19"
    show Chr20 = "chr20"
    show Chr21 = "chr21"
    show Chr22 = "chr22"
    show ChrX = "chrX"
    show ChrY = "chrY"
    show ChrM = "chrM"
    show NONE = "*"

-- 170508 add alternate "show" function to print chromosome names based on
-- config flag setting
showChrom :: ChromNameFormat -> UChr -> String
showChrom fmt chr = case (fmt, chr) of
    (GRC,Chr1)   -> "1"
    (GRC,Chr2)   -> "2"
    (GRC,Chr3)   -> "3"
    (GRC,Chr4)   -> "4"
    (GRC,Chr5)   -> "5"
    (GRC,Chr6)   -> "6"
    (GRC,Chr7)   -> "7"
    (GRC,Chr8)   -> "8"
    (GRC,Chr9)   -> "9"
    (GRC,Chr10)  -> "10"
    (GRC,Chr11)  -> "11"
    (GRC,Chr12)  -> "12"
    (GRC,Chr13)  -> "13"
    (GRC,Chr14)  -> "14"
    (GRC,Chr15)  -> "15"
    (GRC,Chr16)  -> "16"
    (GRC,Chr17)  -> "17"
    (GRC,Chr18)  -> "18"
    (GRC,Chr19)  -> "19"
    (GRC,Chr20)  -> "20"
    (GRC,Chr21)  -> "21"
    (GRC,Chr22)  -> "22"
    (GRC,ChrX)   -> "X"
    (GRC,ChrY)   -> "Y"
    (GRC,ChrM)   -> "MT"
    (UCSC,Chr1)  -> "chr1"
    (UCSC,Chr2)  -> "chr2"
    (UCSC,Chr3)  -> "chr3"
    (UCSC,Chr4)  -> "chr4"
    (UCSC,Chr5)  -> "chr5"
    (UCSC,Chr6)  -> "chr6"
    (UCSC,Chr7)  -> "chr7"
    (UCSC,Chr8)  -> "chr8"
    (UCSC,Chr9)  -> "chr9"
    (UCSC,Chr10) -> "chr10"
    (UCSC,Chr11) -> "chr11"
    (UCSC,Chr12) -> "chr12"
    (UCSC,Chr13) -> "chr13"
    (UCSC,Chr14) -> "chr14"
    (UCSC,Chr15) -> "chr15"
    (UCSC,Chr16) -> "chr16"
    (UCSC,Chr17) -> "chr17"
    (UCSC,Chr18) -> "chr18"
    (UCSC,Chr19) -> "chr19"
    (UCSC,Chr20) -> "chr20"
    (UCSC,Chr21) -> "chr21"
    (UCSC,Chr22) -> "chr22"
    (UCSC,ChrX)  -> "chrX"
    (UCSC,ChrY)  -> "chrY"
    (UCSC,ChrM)  -> "chrM"
    otherwise    -> "*"

-- 170510 try parsing either GRC or UCSC chromosome names
uchrparser :: A.Parser UChr
uchrparser =      (A.string "chr10" >> return Chr10)
        <|>  (A.string "chr11" >> return Chr11)
        <|>  (A.string "chr12" >> return Chr12)
        <|>  (A.string "chr13" >> return Chr13)
        <|>  (A.string "chr14" >> return Chr14)
        <|>  (A.string "chr15" >> return Chr15)
        <|>  (A.string "chr16" >> return Chr16)
        <|>  (A.string "chr17" >> return Chr17)
        <|>  (A.string "chr18" >> return Chr18)
        <|>  (A.string "chr19" >> return Chr19)
        <|>  (A.string "chr20" >> return Chr20)
        <|>  (A.string "chr21" >> return Chr21)
        <|>  (A.string "chr22" >> return Chr22)
        <|>  (A.string "chr1" >> return Chr1)
        <|>  (A.string "chr2" >> return Chr2)
        <|>  (A.string "chr3" >> return Chr3)
        <|>  (A.string "chr4" >> return Chr4)
        <|>  (A.string "chr5" >> return Chr5)
        <|>  (A.string "chr6" >> return Chr6)
        <|>  (A.string "chr7" >> return Chr7)
        <|>  (A.string "chr8" >> return Chr8)
        <|>  (A.string "chr9" >> return Chr9)
        <|>  (A.string "chrX" >> return ChrX)
        <|>  (A.string "chrY" >> return ChrY)
        <|>  (A.string "chrM" >> return ChrM)
        <|>  (A.string "10" >> return C10)
        <|>  (A.string "11" >> return C11)
        <|>  (A.string "12" >> return C12)
        <|>  (A.string "13" >> return C13)
        <|>  (A.string "14" >> return C14)
        <|>  (A.string "15" >> return C15)
        <|>  (A.string "16" >> return C16)
        <|>  (A.string "17" >> return C17)
        <|>  (A.string "18" >> return C18)
        <|>  (A.string "19" >> return C19)
        <|>  (A.string "20" >> return C20)
        <|>  (A.string "21" >> return C21)
        <|>  (A.string "22" >> return C22)
        <|>  (A.string "1" >> return C1)
        <|>  (A.string "2" >> return C2)
        <|>  (A.string "3" >> return C3)
        <|>  (A.string "4" >> return C4)
        <|>  (A.string "5" >> return C5)
        <|>  (A.string "6" >> return C6)
        <|>  (A.string "7" >> return C7)
        <|>  (A.string "8" >> return C8)
        <|>  (A.string "9" >> return C9)
        <|>  (A.string "X" >> return CX)
        <|>  (A.string "Y" >> return CY)
        <|>  (A.string "MT" >> return CMT)
        <|>  (A.string "GL000207.1" >> return GL000207P1)
        <|>  (A.string "GL000226.1" >> return GL000226P1)
        <|>  (A.string "GL000229.1" >> return GL000229P1)
        <|>  (A.string "GL000231.1" >> return GL000231P1)
        <|>  (A.string "GL000210.1" >> return GL000210P1)
        <|>  (A.string "GL000239.1" >> return GL000239P1)
        <|>  (A.string "GL000235.1" >> return GL000235P1)
        <|>  (A.string "GL000201.1" >> return GL000201P1)
        <|>  (A.string "GL000247.1" >> return GL000247P1)
        <|>  (A.string "GL000245.1" >> return GL000245P1)
        <|>  (A.string "GL000197.1" >> return GL000197P1)
        <|>  (A.string "GL000203.1" >> return GL000203P1)
        <|>  (A.string "GL000246.1" >> return GL000246P1)
        <|>  (A.string "GL000249.1" >> return GL000249P1)
        <|>  (A.string "GL000196.1" >> return GL000196P1)
        <|>  (A.string "GL000248.1" >> return GL000248P1)
        <|>  (A.string "GL000244.1" >> return GL000244P1)
        <|>  (A.string "GL000238.1" >> return GL000238P1)
        <|>  (A.string "GL000202.1" >> return GL000202P1)
        <|>  (A.string "GL000234.1" >> return GL000234P1)
        <|>  (A.string "GL000232.1" >> return GL000232P1)
        <|>  (A.string "GL000206.1" >> return GL000206P1)
        <|>  (A.string "GL000240.1" >> return GL000240P1)
        <|>  (A.string "GL000236.1" >> return GL000236P1)
        <|>  (A.string "GL000241.1" >> return GL000241P1)
        <|>  (A.string "GL000243.1" >> return GL000243P1)
        <|>  (A.string "GL000242.1" >> return GL000242P1)
        <|>  (A.string "GL000230.1" >> return GL000230P1)
        <|>  (A.string "GL000237.1" >> return GL000237P1)
        <|>  (A.string "GL000233.1" >> return GL000233P1)
        <|>  (A.string "GL000204.1" >> return GL000204P1)
        <|>  (A.string "GL000198.1" >> return GL000198P1)
        <|>  (A.string "GL000208.1" >> return GL000208P1)
        <|>  (A.string "GL000191.1" >> return GL000191P1)
        <|>  (A.string "GL000227.1" >> return GL000227P1)
        <|>  (A.string "GL000228.1" >> return GL000228P1)
        <|>  (A.string "GL000214.1" >> return GL000214P1)
        <|>  (A.string "GL000221.1" >> return GL000221P1)
        <|>  (A.string "GL000209.1" >> return GL000209P1)
        <|>  (A.string "GL000218.1" >> return GL000218P1)
        <|>  (A.string "GL000220.1" >> return GL000220P1)
        <|>  (A.string "GL000213.1" >> return GL000213P1)
        <|>  (A.string "GL000211.1" >> return GL000211P1)
        <|>  (A.string "GL000199.1" >> return GL000199P1)
        <|>  (A.string "GL000217.1" >> return GL000217P1)
        <|>  (A.string "GL000216.1" >> return GL000216P1)
        <|>  (A.string "GL000215.1" >> return GL000215P1)
        <|>  (A.string "GL000205.1" >> return GL000205P1)
        <|>  (A.string "GL000219.1" >> return GL000219P1)
        <|>  (A.string "GL000224.1" >> return GL000224P1)
        <|>  (A.string "GL000223.1" >> return GL000223P1)
        <|>  (A.string "GL000195.1" >> return GL000195P1)
        <|>  (A.string "GL000212.1" >> return GL000212P1)
        <|>  (A.string "GL000222.1" >> return GL000222P1)
        <|>  (A.string "GL000200.1" >> return GL000200P1)
        <|>  (A.string "GL000193.1" >> return GL000193P1)
        <|>  (A.string "GL000194.1" >> return GL000194P1)
        <|>  (A.string "GL000225.1" >> return GL000225P1)
        <|>  (A.string "GL000192.1" >> return GL000192P1)
        <|>  (A.string "NC_007605" >> return NC_007605)
        <|>  (A.string "*" >> return NONE)


-- 180206 include option for providing primer coords in BED or BEDPE format
-- in addition to master file
optargs :: Parser Opts
optargs = Opts
    <$> switch
        ( short 'b'
       <> long "bedpe"
       <> help "add this switch to use BEDPE coordinate input format (default format is master file)"
        )
    <*> argument str (metavar "PRIMER_COORDS_INFILE")
    <*> argument str (metavar "SAM_INFILE")
    <*> argument str (metavar "OUTPUT_SAM_FILENAME")

-- record to store command line arguments
data Opts = Opts { coordfileformat :: Bool
                 , incoordsfile :: String
                 , insamfile :: String
                 , outfilename :: String
                 } deriving (Show, Eq)



-- 170927 parse master file for primer and target intervals
masterparser :: A.Parser MasterRecord
masterparser = do
    tbed <- bedrecparser
    A.skipSpace
    fpstart <- A.decimal
    A.skipSpace
    fpend <- A.decimal
    A.skipSpace
    fpname <- txtfieldp
    A.skipSpace
    rpstart <- A.decimal
    A.skipSpace
    rpend <- A.decimal
    A.skipSpace
    rpname <- txtfieldp
    A.skipSpace
    fseq <- seqp
    A.skipSpace
    rseq <- seqp
    return  MasterRecord { mchrom = bedchr tbed
                         , mtargetstart = bedstart tbed
                         , mtargetend = bedend tbed
                         , mampname = bedname tbed
                         , mforstart = fpstart
                         , mforend = fpend
                         , mrevstart = rpstart
                         , mrevend = rpend
                         , mforname = fpname
                         , mrevname = rpname
                         , mforseq = fseq
                         , mrevseq = rseq
                         }

bedrecparser :: A.Parser BedRecord
bedrecparser = do
    chr <- uchrparser
    A.skipSpace
    bstart <- A.decimal
    A.skipSpace
    bend <- A.decimal
    A.skipSpace
    bname <- txtfieldp
    let bedrec = defaultBed { bedchr = chr
                            , bedstart = bstart
                            , bedend = bend
                            , bedname = bname
                            }
    return bedrec

-- 20180206
bedPEparser :: A.Parser BEDPE
bedPEparser = do
    c1 <- uchrparser
    A.skipSpace
    s1 <- A.decimal
    A.skipSpace
    e1 <- A.decimal
    A.skipSpace
    c2 <- uchrparser
    A.skipSpace
    s2 <- A.decimal
    A.skipSpace
    e2 <- A.decimal
    bname <- txtfieldp
    -- ignore any additional columns present
    return $ BEDPE c1 s1 e1 c2 s2 e2 bname

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- read SAM file and parse header and alignments, returning both in tuple
-- NOTE: as of October 2017 these functions are deprecated due to use of the
-- Conduit library for lower memory usage. They remain in the library in the
-- event they are needed for future use.
------------------------------------------------------------------------------
------------------------------------------------------------------------------
getSAM2 :: FilePath -> IO SAM
getSAM2 fp = parseSAM <$> B.lines <$> B.readFile fp

parseSAM :: [B.ByteString] -> SAM
parseSAM lns =
    let (hdr, bdytxt) = partition (\x -> (B.head x) == '@') lns
        alignments = parseAlns bdytxt
    in (hdr, alignments)

-- NOTE: changed to make use of U.rights to leave out parsing failures
parseAlns :: [B.ByteString] -> Alignments
parseAlns as = U.rights $ (A.parseOnly alnparser <$> as)

parseAln as = A.parseOnly alnparser as

-- 170316 parse with logging of parse failures
-- BUG: mis-assignment of rname (not sure if this is the source) -- 170510 FIXED??
-- 170510 report lines which failed SAM parsing into file containing actual lines (not line numbers)
getSAM :: FilePath -> IO SAM
getSAM fp = do
    flines <- B.lines <$> B.readFile fp
    let (hdr, bdy) = partition (\x -> (B.head x) == '@') flines
        nbdy = length bdy
        as = parseAln <$> bdy
        lks = [0..] :: [Int]
        amap = M.fromList $ zip lks as
        (succm, failm) = M.partition U.isRight amap
        succs = U.rights $ snd <$> M.toAscList succm
        faildlineixs = M.keys failm
        acnt = length succs
        failcnt = length faildlineixs
        parsestatus = parsechkSAM nbdy acnt faildlineixs
        parsefaildLines = B.unlines $ (bdy !!) <$> faildlineixs
    putStrLn parsestatus
    B.writeFile "samparsefails.log" parsefaildLines
    writeFile "sam_parsing.log" parsestatus
    return (hdr, succs)

-- read BED file into list of BedRecord (primer and target intervals)
--getBED :: FilePath -> IO BED
getBED fp = do
    flines <- B.lines <$> B.readFile fp
    let nr = length flines
        bs = parseBED <$> flines
        lks = [0..] :: [Int]
        bmap = M.fromList $ zip lks bs
        (succm, failm) = M.partition U.isRight bmap
        succs = V.fromList $ U.rights $ snd <$> M.toAscList succm
        faildlineixs = M.keys failm
        bcnt = V.length succs
        failcnt = length faildlineixs
        parsestatus = parsechkBED nr bcnt faildlineixs
        parsefaildLines = B.unlines $ (flines !!) <$> faildlineixs
    putStrLn parsestatus
    B.writeFile "bedparsefails.log" parsefaildLines
    writeFile "primer_parsing.log" parsestatus
    return succs

getBEDPE :: FilePath -> IO [BEDPE]
getBEDPE fp = do
    flines <- B.lines <$> B.readFile fp
    let nr = length flines
        bs = parseBEDPE <$> flines -- [BEDPE]
        lks = [0..] :: [Int]
        bmap = M.fromList $ zip lks bs
        (succm, failm) = M.partition U.isRight bmap
        succs = U.rights $ snd <$> M.toAscList succm
        faildlineixs = M.keys failm
        bcnt = length succs
        failcnt = length faildlineixs
        parsestatus = parsechkBED nr bcnt faildlineixs
        parsefaildLines = B.unlines $ (flines !!) <$> faildlineixs
    putStrLn parsestatus
    B.writeFile "bedPEparsefails.log" parsefaildLines
    writeFile "primer_BEDPE_parsing.log" parsestatus
    return succs

parseBED bs = A.parseOnly bedrecparser bs
parseBEDPE bs = A.parseOnly bedPEparser bs

-- 170510 write parse fail lines to log file (don't print line nums to stdout)
parsechkSAM :: Int -> Int -> [Int] -> String
parsechkSAM numrec numsucc failedlines
    | numrec == numsucc = reportOK
    | otherwise = reportFails
        where reportOK = "all " ++ (show numrec)
                                ++ " SAM alignments parsed successfully.\n"
              parsediff = length failedlines
              --flinestrs = show <$> failedlines
              --fstring = intercalate "\n" flinestrs
              reportFails = "WARNING: "
                             ++ (show parsediff)
                             ++ " of " ++ (show numrec)
                             ++ " SAM alignments failed to parse (see samparsefails.log)"
                             ++ "\n"
                             -- ++ fstring ++ "\n"

parsechkBED :: Int -> Int -> [Int] -> String
parsechkBED numrec numsucc failedlines
    | numrec == numsucc = reportOK
    | otherwise = reportFails
        where reportOK = "all " ++ (show numrec)
                                ++ " primer BEDPE records parsed successfully.\n"
              parsediff = length failedlines
              -- flinestrs = show <$> failedlines
              -- fstring = intercalate "\n" flinestrs
              reportFails = "WARNING: "
                             ++ (show parsediff)
                             ++ " of " ++ (show numrec)
                             ++ " primer BEDPE records failed to parse (see bedparsefails.log)"
                             ++ "\n"


-- 170926 read and parse master file
getMasterFile :: FilePath -> IO [MasterRecord]
getMasterFile fp = do
    mlines <- B.lines <$> B.readFile fp
    let nr = length mlines
        mrecs = (A.parseOnly masterparser) <$> mlines
        lks = [0..] :: [Int]
        mrecmap = M.fromList $ zip lks mrecs
        (succm, failm) = M.partition U.isRight mrecmap
        succs = U.rights $ snd <$> M.toAscList succm
        failedlinenums = M.keys failm
        mcnt = length succs
        failcnt = length failedlinenums
        parseStatus = parsechkMaster nr mcnt failedlinenums
        parsefailrecs = B.unlines $ (mlines !!) <$> failedlinenums
    putStrLn parseStatus
    B.writeFile "masterparsefails.log" $ parsefailrecs
    return succs

parsechkMaster :: Int -> Int -> [Int] -> String
parsechkMaster numrec numsucc failedlines
    | numrec == numsucc = reportOK
    | otherwise = reportFails
        where reportOK = "all " ++ (show numrec)
                                ++ " master records parsed successfully.\n"
              parsediff = length failedlines
              reportFails = "WARNING: "
                             ++ (show parsediff)
                             ++ " of " ++ (show numrec)
                             ++ " master records failed to parse (see masterparsefails.log)"
                             ++ "\n"

masterToFPrimerBED :: [MasterRecord] -> BED
masterToFPrimerBED mrecs = V.fromList
                         $ getFPrimerBEDfromMaster <$> (sort mrecs)

masterToRPrimerBED :: [MasterRecord] -> BED
masterToRPrimerBED mrecs = V.fromList
                         $ getRPrimerBEDfromMaster <$> (sort mrecs)

masterToTargetBED :: [MasterRecord] -> BED
masterToTargetBED mrecs = V.fromList
                         $ getTargetBEDfromMaster <$> (sort mrecs)

masterRecToTargetBED :: MasterRecord -> BedRecord
masterRecToTargetBED m = defaultBed { bedchr = mchrom m
                                    , bedstart = mtargetstart m
                                    , bedend = mtargetend m
                                    , bedname = mampname m
                                    }

getFPrimerBEDfromMaster :: MasterRecord -> BedRecord
getFPrimerBEDfromMaster mr =
    let c = mchrom mr
        fpstart = mforstart mr
        fpend = mforend mr
        fpname = mforname mr
    in defaultBed { bedchr = c
                  , bedstart = fpstart
                  , bedend = fpend
                  , bedname = fpname
                  }

getRPrimerBEDfromMaster :: MasterRecord -> BedRecord
getRPrimerBEDfromMaster mr =
    let c = mchrom mr
        rpstart = mrevstart mr
        rpend = mrevend mr
        rpname = mrevname mr
    in defaultBed { bedchr = c
                  , bedstart = rpstart
                  , bedend = rpend
                  , bedname = rpname
                  }

getTargetBEDfromMaster :: MasterRecord -> BedRecord
getTargetBEDfromMaster mr =
    let c = mchrom mr
        tstart = mtargetstart mr
        tend = mtargetend mr
        tname = mampname mr
    in defaultBed { bedchr = c
                  , bedstart = tstart
                  , bedend = tend
                  , bedname = tname
                  }

-- new alignment parser which allows parse failures to be dropped
-- from Alignments
alnparser :: A.Parser AlignedRead
alnparser = do
    qn <- txtfieldp
    A.space
    f <- A.decimal
    A.space
    chr <- uchrparser
    A.space
    p <- A.decimal
    A.space
    mpscore <- A.decimal
    A.space
    cig <- txtfieldp
    A.space
    nchr <- txtfieldp
    A.space
    pn <- A.decimal
    A.space
    tl <- A.double
    A.space
    seq <- txtfieldp
    A.space
    qual <- txtfieldp
    A.space
    optfs <- optfieldsp -- optfieldpEOL : only for use with conduitParserEither
    let flag = f
        strand = case testBit flag 4 of
            True  -> B.pack "-"
            False -> B.pack "+"
        prd = testBit flag 0
        mapd = not $ testBit flag 2
        tln = floor $ tl
        cigm = exrights $ parseCigar cig
        end = (sumMatches cigm) + p - 1
        optfieldstr = B.intercalate "\t" optfs -- B.pack optfs
        midstr = parsemIDstring optfieldstr
        a = defaultAlignment { qname = qn
                             , flag = f
                             , rname = chr
                             , pos = p - 1 -- modify as necessary
                             , trimdpos = p - 1
                             , endpos = end
                             , trimdendpos = end
                             , mapqual = mpscore
                             , cigar = cig
                             , cigmap = cigm
                             , rnext = nchr
                             , pnext = pn - 1 -- modify as necessary
                             , tlen = tln
                             , refseq = seq
                             , basequal = qual
                             , optfields = optfieldstr
                             , strand = strand
                             , paired = prd
                             , mapped = mapd
                             , mid = midstr
                             }
    return a

-- 171017 header parser more suited for conduit and parsing header into AlignedRead
hdralnparser :: A.Parser AlignedRead
hdralnparser = do
    hlines <- A.many1 samhdrparser
    return $ defaultAlignment { headerstrings = hlines
                              , isheader = True
                              , qname = "HEADERLINE"
                              }

-- SAM header chromosome names parser (parses chromosome names ONLY for now)
-- NOTE: for simplicity, only canonical chromosome names are parsed at this time.
hdrchromp :: A.Parser UChr
hdrchromp = do
    A.string "@SQ"
    A.skipSpace
    A.string "SN:"
    c <- uchrparser
    A.many1 A.anyChar
    return c

samhdrparser :: A.Parser B.ByteString
samhdrparser = do
    atprefix <- A.char '@'
    restofhdrline <- A.takeTill (A.inClass "\n\r")
    let hdrln = B.append "@" restofhdrline
    return hdrln

-- 171204 parser for header line with newline char present
samhdrparserEOL :: A.Parser B.ByteString
samhdrparserEOL = do
    atprefix <- A.char '@'
    fields <- A.takeTill (A.inClass "\n\r")
    A.endOfLine
    let hdrln = B.append "@" fields
    return hdrln

-- 171204 Attoparsec version of getRight (below)
-- use of Conduit to read input alignment file precludes use of "rights" in Data.Either
-- unless/until we come up with a more idiomatic implementation of the parsing step.
rightOrDefault :: Either String AlignedRead -> AlignedRead
rightOrDefault e = case e of
    Left _ -> defaultAlignment
    Right a -> a

-- 171017 extract parse success from Either (as singleton, for conduit)
getRight :: Either t (a, AlignedRead) -> AlignedRead
getRight e = case e of
    Left _  -> defaultAlignment
    Right x -> snd x

-- print parse failures to log file (try reporting attoparsec output to start)
-- getRightM :: (P.MonadIO m, Show t) => Either t (a, AlignedRead) -> IO AlignedRead
getRightM ev = case ev of
    Left e  -> do
        -- let epos = errorPosition e
        liftIO $ writeFile "samparsefails.log" $ show e
        return defaultAlignment
    Right a -> return $ snd a

-- cigar parsers
starcigarP = A.string "*" >> return [(0, star)]
    where star = "*" :: B.ByteString

cigarP = A.many1 ((,) <$> A.decimal <*> (A.takeWhile (A.inClass "MIDNSHP=X")))

parseCigar = A.parseOnly $ cigarP <|> starcigarP

-- used on single cigar ByteString
getAlignedLength cigcol = sumMatches ciglist
    where ciglist = head $ U.rights $ cig
          cig = (parseCigar cigcol) : []

sumMatches :: [(Integer, B.ByteString)] -> Integer
sumMatches cigs = sum [ x | (x, y) <- cigs, y == "M" || y == "I" ]

-- 161023 check if cigar "length" equal to read length (refseq)
checkcigseqlen :: AlignedRead -> Bool
checkcigseqlen a
    | cigmatchlen == refseqlen = True
    | otherwise = False
        where cigmatchlen = sum [ x | (x, y) <- tcmap, y == "M"
                                                    || y == "I"
                                                    || y == "S" ]
              tcmap = mapcig $ trimdcigar a
              refseqlen = genericLength $ B.unpack $ refseq a

-- 170206 filter nonmatching and zero-length CIGAR/sequence lengths
-- however: keep "*"
checkcigseqlen2 :: AlignedRead -> Bool
checkcigseqlen2 a
    | cigar a == "*" = True
    | (cigmatchlen == refseqlen) && (matchcnt > 0) = True
    | otherwise = False
        where matchcnt = sumMatches tcmap
              cigmatchlen = sum [ x | (x, y) <- tcmap, y == "M"
                                                    || y == "I"
                                                    || y == "S" ]
              tcmap = mapcig $ trimdcigar a
              refseqlen = genericLength $ B.unpack $ refseq a

-- 171018 new CIGAR sequence length check which only checks non-header AlignedRead recs
checkCigarSeqlen :: AlignedRead -> AlignedRead
checkCigarSeqlen a
    | (isheader a) = a
    | trimdcigarmap == Nothing = defaultAlignment
    | (cigseqlenHdrPassTest a) == True = a
    | otherwise = defaultAlignment
        where trimdcigarmap = safemapcig $ trimdcigar a

cigseqlenHdrPassTest :: AlignedRead -> Bool
cigseqlenHdrPassTest a
    -- | (isheader a) = True
    | cigar a == "*" = True
    | (cigmatchlen == refseqlen) && (matchcnt > 0) = True
    | otherwise = False
        where matchcnt = sumMatches tcmap
              cigmatchlen = sum [ x | (x, y) <- tcmap, y == "M"
                                                    || y == "I"
                                                    || y == "S" ]
              tcmap = mapcig $ trimdcigar a
              refseqlen = genericLength $ B.unpack $ refseq a

-- translate CIGAR string into aligned interval with CIGAR op labels
-- 171018 safer version due to head on empty list error
safemapcig :: B.ByteString -> Maybe CigarMap
safemapcig cigstr
    | (length cigtups) > 0 = Just $ head cigtups
    | otherwise = Nothing
        where cigtups = U.rights $ [parseCigar cigstr] -- [(Integer, B.ByteString)]

-- translate CIGAR string into aligned interval with CIGAR op labels
mapcig :: B.ByteString -> CigarMap
mapcig cigstr =
    let cigtups = exrights $ parseCigar cigstr -- [(Integer, B.ByteString)]
    in cigtups

-- tuple map helper functions
mapfst :: Num a => (a -> c) -> (a, b) -> (c, b)
mapfst f (x, y) = (f x, y)

-- parsers and helper functions
readint = read :: String -> Integer
rdint = read :: String -> Int

parseint = A.parseOnly A.decimal

parsedbl = A.parseOnly A.double

parsesignedint i = floor $ exrights $ parsedbl i :: Integer

spaces = A.many1 A.skipSpace


--txtfieldp = A.takeTill A.isHorizontalSpace
txtfieldp = A.takeTill A.isSpace

optfieldsp = A.sepBy' txtfieldp A.space -- parses correctly

-- 171017 prevent entire SAM file being parsed into optional fields field of
-- first AlignedRead
optfieldpEOL = A.manyTill (A.satisfy (A.notInClass "\r\n")) A.endOfLine

optfieldsp2 = A.sepBy' optfieldpEOL A.space

optfieldsp3 = A.sepBy' txtfieldp A.space

optfieldp :: A.Parser B.ByteString
optfieldp = optnmp <|> optasp <|> optxsp

-- number of mismatches in alignment
optnmp :: A.Parser B.ByteString
optnmp = do
    id <- optidp
    skipcolon
    dtype <- itypep
    skipcolon
    cnt <- A.decimal
    let cntstr = B.pack $ show cnt
        nmoptstr = B.intercalate ":" [id, dtype, cntstr]
    return nmoptstr

-- alignment score of this alignment
optasp :: A.Parser B.ByteString
optasp = do
    id <- optidp
    skipcolon
    dtype <- itypep
    skipcolon
    score <- A.decimal
    let scorestr = B.pack $ show score
        asoptstr = B.intercalate ":" [id, dtype, scorestr]
    return asoptstr

-- alignment score of best alternate alignment (if AS is >= XS
-- then this aln is primary)
optxsp :: A.Parser B.ByteString
optxsp = do
    id <- optidp
    skipcolon
    dtype <- itypep
    skipcolon
    score <- A.decimal
    let scorestr = B.pack $ show score
        xsoptstr = B.intercalate ":" [id, dtype, scorestr]
    return xsoptstr

skipcolon = A.skipWhile (\x -> x == ':')

optidp :: A.Parser B.ByteString
optidp =    ("NM" >> return "NM")
        <|> ("AS" >> return "AS")
        <|> ("XS" >> return "XS")

itypep :: A.Parser B.ByteString
itypep = "i" >> return "i"

-- 171017 mID parsing (additional mID processing to be added to master as needed)
mIDp :: A.Parser B.ByteString
mIDp = do
    tag <- A.string "RX:Z"
    mid <- A.takeWhile (A.inClass "ACGTN")
    return mid

seqp :: A.Parser B.ByteString
seqp = A.takeWhile (A.inClass "ACTGatcgN")

parsemIDstring :: B.ByteString -> B.ByteString
parsemIDstring bs
    | length midvals <= 0 = "NOMIDSTRING"
    | otherwise = head midvals
        where midvals = U.rights $ (A.parseOnly mIDp) <$> (B.words bs)


intgr2int :: Integer -> Int
intgr2int n = fromIntegral n

-- convert decimal flag to binary representation for parsing flag
toBin :: Int -> String
toBin n = concatMap show $ digits 2 n

toBinInt :: String -> Int
toBinInt n = read n

-- possible data structure for SAM flag
data SAMFlag = SAMFlag { pairedRead :: Bool
                       , mateMapped :: Bool
                       , notMapped  :: Bool
                       , pairNotMapped :: Bool
                       , queryPlus :: Bool
                       , mateNeg :: Bool
                       , fstInPair :: Bool
                       , sndInPair :: Bool
                       , notPrimary :: Bool
                       , failedQC :: Bool
                       , dupRead :: Bool
                       , intflag :: Int
                       } deriving (Show, Eq)

readSAMFlag :: Int -> SAMFlag
readSAMFlag flag =
       SAMFlag { pairedRead = testBit flag 0
               , mateMapped = testBit flag 1
               , notMapped = testBit flag 2
               , pairNotMapped = testBit flag 3
               , queryPlus = testBit flag 4
               , mateNeg =  testBit flag 5
               , fstInPair = testBit flag 6
               , sndInPair = testBit flag 7
               , notPrimary = testBit flag 8
               , failedQC = testBit flag 9
               , dupRead = testBit flag 10
               , intflag = flag
               }

-- select element of nested vector
getcol :: Int -> V.Vector (V.Vector a) -> V.Vector a
getcol n txt = fmap (V.! n) txt

-- short notation for selecting a vector element by index
ix :: Int -> V.Vector a -> a
ix i v = (V.! i) v

-- parsing hack
exrights x = head $ U.rights $ [x]

-- safer exrights
exrights2 :: [a] -> Maybe a
exrights2 xs
    | length xs == 0 = Nothing
    | otherwise = Just (head xs)

getlengths seqs = fmap B.length seqs


------------------------------------------------------------------------------
----------------------  HIGH-LEVEL FUNCTIONS IN MAIN -------------------------
------------------------------------------------------------------------------

-- element (single alignment)-wise primer trimming (for use with conduit)
trimprimersE :: CMap -> CMap -> AlignedRead -> AlignedRead
trimprimersE fmap rmap a =
    let intaln = addprimerints fmap rmap a
        trimd = trimalignment intaln
    in trimd

-- 171017 function to convert conduit ZipSink stream to output and write
-- output to file in SAM format.
printAlnStreamToFile :: P.MonadResource m => FilePath -> P.ConduitM AlignedRead c m ()
printAlnStreamToFile outfile = P.mapC printAlignmentOrHdr
                          P..| P.unlinesAsciiC
                          P..| P.sinkFile outfile

-- 171017 calculate trim stats and print to stdout (TODO: print full stats to file)
calculateTrimStats :: P.ConduitM AlignedRead c (P.ResourceT IO) Integer
calculateTrimStats = P.filterC (\x -> trimdflag x) P..| P.lengthC

-- 170926 calculate and populate amplicon target BED field
addprimerints :: CMap -> CMap -> AlignedRead -> AlignedRead
addprimerints fpmap rpmap aln =
    let achr = rname aln
        start = pos aln
        end = endpos aln
        fbedmap = justchrmaps $ [M.lookup achr fpmap] -- possibly []
        rbedmap = justchrmaps $ [M.lookup achr rpmap] -- possibly []
        starthits = justbedmaps $ (bedmaplookup start) <$> fbedmap
        endhits = justbedmaps $ (bedmaplookup end) <$> rbedmap
        allhits = sort (starthits ++ endhits) -- [BedRecord]
        flagstatus = setpintflag allhits -- set True if alignment intersects one or more primers
    in aln { fint = starthits
           , rint = endhits
           , pintflag = flagstatus
           }

setpintflag hits
    | (length hits) > 0 = True
    | otherwise = False


---------- Functions to trim AlignedRead if it intersects >=1 primer ---------
------------------------------------------------------------------------------

-- 161017 modify AlignedRead to trim primer sequence by softclips
trimalignment :: AlignedRead -> AlignedRead
trimalignment a
    | (fint a /= []) && (rint a /= []) = btrimdaln
    | (fint a /= []) && (rint a == []) = ftrimdaln
    | (fint a == []) && (rint a /= []) = rtrimdaln
    | (fint a == []) && (rint a == []) = a { trimdcigar = cigar a }
    | otherwise = a { trimdcigar = cigar a }
         where ftrimdaln = trimfwd a
               rtrimdaln = trimrev a
               btrimdaln = trimboth a

-- for alignments intersecting a "forward" primer only ( ref (+)-strand orientation)
trimfwd :: AlignedRead -> AlignedRead
trimfwd a =
    let newpos = bedend $ head $ fint a
        as = pos a
        fdiff = newpos - as -- handle alignment ending "inside" primer
        oldcigar = cigar a
        newcig = updateCigF fdiff oldcigar
        newcigmap = exrights $ parseCigar newcig
        tpos = (if (fdiff < 0) then as else newpos)
    in a { trimdpos = tpos
         , trimdcigar = newcig
         , trimdcigmap = newcigmap
         , trimdflag = if (as /= tpos) then True else False
         }

-- for alignments intersecting a "reverse" primer only ( ref (+)-strand orientation)
trimrev :: AlignedRead -> AlignedRead
trimrev a =
    let newendpos = bedstart $ head $ rint a
        ae = endpos a
        rdiff = (ae - newendpos) -- handle alignment ending "inside" primer
        oldcigar = cigar a
        newcig = updateCigR rdiff oldcigar
        newcigmap = exrights $ parseCigar newcig
        tendpos = (if (rdiff < 0) then ae else newendpos)
    in a { trimdendpos = tendpos
         , trimdcigar = newcig
         , trimdcigmap = newcigmap
         , trimdflag = if (ae /= tendpos) then True else False
         }

-- for alignments with primer intersections at both ends
trimboth :: AlignedRead -> AlignedRead
trimboth a =
    let newpos = bedend $ head $ fint a
        newendpos = bedstart $ head $ rint a
        as = pos a
        ae = endpos a
        fdiff = (newpos - as)
        rdiff = (ae - newendpos)
        oldcigar = cigar a
        newcig = updateCigB fdiff rdiff oldcigar
        newcigmap = exrights $ parseCigar newcig
        tpos = (if (fdiff < 0) then as else newpos)
        tendpos = (if (rdiff < 0) then ae else newendpos)
    in a { trimdpos = tpos
         , trimdendpos = tendpos
         , trimdcigar = newcig
         , trimdcigmap = newcigmap
         , trimdflag = if ((as /= tpos) || (ae /= tendpos))
                       then True
                       else False
         }


-- UPDATE 17-02-01 simplify CIGAR arithmetic by adjusting read coords to match
--                 reference coords (D CIGAR chars don't increment coords, and
--                 I chars double-increment coords)
updateCigF :: Integer -> B.ByteString -> B.ByteString
updateCigF fdiff cigar
    | snd (head cmap) == "*" = "*"
    | fdiffi <= 0 = cigar
    | (nopadlen - fdiffi) == 0 = cigar
    | ((nopadlen - fdiffi) > 0) = newcig
    | otherwise = "*"
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              fdiffi = intgr2int fdiff
              cignoclip = filter noclip cmap
              cigexp = expandcigar2 cignoclip -- [(Int, B.ByteString)] assoc. list
              nopadlen = length cigexp
              adjcig = adjustcrds cigexp
              ftrim = taketrim fdiff adjcig
              trimDs = countDs ftrim
              trimDcorr = intgr2int trimDs
              adjix = (genericLength ftrim) + trimDs
              remcigraw = drop trimDcorr $ trimrem adjix cigexp
              newss = B.replicate (length ftrim) 'S'
              newcigarcore = B.append newss $ B.concat $ snd <$> remcigraw
              newcigar = B.append (B.append fSs newcigarcore) rSs
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar

updateCigR :: Integer -> B.ByteString -> B.ByteString
updateCigR rdiff cigar
    | snd (head cmap) == "*" = "*"
    | rdiffi <= 0 = cigar
    | (nopadlen - rdiffi) == 0 = cigar
    | ((nopadlen - rdiffi) > 0) = newcig
    | otherwise = "*"
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              rdiffi = intgr2int rdiff
              cignoclip = filter noclip cmap
              cigexpR = expandcigar2 $ reverse cignoclip
              nopadlen = length cigexpR
              adjcigR = adjustcrds cigexpR
              rtrim = taketrim rdiff adjcigR
              trimDs = countDs rtrim
              trimDcorr = intgr2int trimDs
              adjix = (genericLength rtrim) + trimDs
              remcigraw = drop trimDcorr $ trimrem adjix cigexpR
              remDs = intgr2int $ countDs remcigraw
              remcigDadj = reverse $ drop remDs remcigraw
              newss = B.replicate ((length rtrim)
                                  + trimDcorr
                                  + remDs) 'S'
              newcigarcore = B.append (B.concat $ snd <$> remcigDadj) newss
              newcigar = B.append fSs (B.append newcigarcore rSs)
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar

updateCigB :: Integer -> Integer -> B.ByteString -> B.ByteString
updateCigB fdiff rdiff cigar
    | snd (head cmap) == "*" = "*"
    | fdiffi <= 0 = updateCigR rdiff cigar
    | rdiffi <= 0 = updateCigF fdiff cigar
    | (nopadlen - fdiffi) == 0 = updateCigR rdiff cigar
    | (nopadlen - rdiffi) == 0 = updateCigF fdiff cigar
    | ((nopadlen - fdiffi - rdiffi) > 0) = newcig
    | otherwise = "*"
        where cmap = mapcig cigar
              grps = B.group $ expandcigar cmap
              nohardgrps = B.group $ expandcigar $ filter nohardclip cmap
              fHs = B.filter (== 'H') $ head grps -- maybe ""
              rHs = B.filter (== 'H') $ last grps -- maybe ""
              fSs = B.filter (== 'S') $ head nohardgrps
              rSs = B.filter (== 'S') $ last nohardgrps
              fdiffi = intgr2int fdiff
              rdiffi = intgr2int rdiff
              -- 5p trim
              cignoclipf = filter noclip cmap
              cigexpf = expandcigar2 cignoclipf
              nopadlen = length cigexpf
              adjcigf = adjustcrds cigexpf
              ftrim = taketrim fdiff adjcigf
              ftrimDs = countDs ftrim
              ftrimDcorr = intgr2int ftrimDs
              adjixf = (genericLength ftrim) + ftrimDs
              remcigrawf = drop ftrimDcorr $ trimrem adjixf cigexpf
              newfss = B.replicate (length ftrim) 'S'
              -- 3p trim3p
              cigexpr = zipWith (\x (i, j) -> (x, j))
                                [1..]
                                (reverse remcigrawf)
              adjcigr = adjustcrds cigexpr
              rtrim = taketrim rdiff adjcigr
              rtrimDs = countDs rtrim
              rtrimDcorr = intgr2int rtrimDs
              adjixr = (genericLength rtrim) + rtrimDs
              remcigrawr = drop rtrimDcorr $ trimrem adjixr cigexpr
              remDsr = intgr2int $ countDs remcigrawr
              remcigrDadj = reverse $ drop remDsr remcigrawr
              newrss = B.replicate
                        ((length rtrim) + rtrimDcorr + remDsr + ftrimDcorr) 'S'
              totfss = B.append fSs newfss
              totrss = B.append rSs newrss
              newcigarcore = B.concat $ snd <$> remcigrDadj
              newcigar = B.append totfss (B.append newcigarcore totrss)
              newfullcigar = B.append fHs (B.append newcigar rHs)
              newcig = contractcigar newfullcigar

showcigar :: (Integer, B.ByteString) -> B.ByteString
showcigar cm =
    let cnt = B.pack $ show $ fst cm
    in B.append cnt (snd cm)

expandcigar :: CigarMap -> B.ByteString
expandcigar cmap = B.concat
                $ [ B.replicate (intgr2int n) $ B.head c
                  | (n, c) <- cmap ]

expandcigar2 :: CigarMap -> [(Integer, B.ByteString)]
expandcigar2 cmap =
    let ops = B.singleton <$> (B.unpack $ expandcigar cmap)
        opslen = genericLength ops
        ixs = [1..opslen] :: [Integer]
    in zipWith (,) ixs ops

-- convert expanded cigar string to standard CIGAR format (after trim update)
contractcigar :: B.ByteString -> B.ByteString
contractcigar longcig
    | longcig == "*" = "*"
    | otherwise = cigarstring
        where grpd = B.group longcig -- maintains order of groups
              shorten x = B.append (B.pack $ show $ B.length x)
                                   (B.singleton $ B.head x)
              cigarstring = B.concat $ shorten <$> grpd

-- 161023 sum all "S", "I", and "M" chars in trimmed sequence (drop all "D")
-- to get accurate CIGAR string
filtpadassoc :: [(Int, B.ByteString)] -> [(Int, B.ByteString)]
filtpadassoc asclist = filter nopadassoc asclist
    where nopadassoc = (\x -> (snd x) /= "D")

-- 161107 remove "D"s before making associative array ([(Int, B.ByteString)])
filtpadcmap :: CigarMap -> CigarMap
filtpadcmap cmap = filter nopad cmap
    where nopad = (\x -> (snd x) /= "D")

-- 161023 remove "H" and "D" from CigarMap (for trimming by soft clipping)
nopadding :: (Integer, B.ByteString) -> Bool
nopadding cp
    | (op /= "H") && (op /= "D") = True
    | otherwise = False
        where op = snd cp

-- 170201
removeDs :: [(Integer, B.ByteString)] -> [(Integer, B.ByteString)]
removeDs cs = filter (\x -> (snd x) /= "D") cs

countDs :: [(Integer, B.ByteString)] -> Integer
countDs cs = genericLength $ filter (\x -> (snd x) == "D") cs

nohardclip :: (Integer, B.ByteString) -> Bool
nohardclip op
    | (snd op) /= "H" = True
    | otherwise = False

noclip :: (Integer, B.ByteString) -> Bool
noclip op
    | ((snd op) /= "H") && ((snd op) /= "S") = True
    | otherwise = False

-- collect cigar string group (usually 'H' or 'S') into ordered list
getNclips :: [B.ByteString] -> Char -> [B.ByteString]
getNclips ciggrps cigchar = filter (\x -> B.length x > 0)
                          $ B.filter (== cigchar) <$> ciggrps

-- 171018 conduit-compatible filters for cigar string match count > 0
-- this function is designed to filter out alignments which have zero matches to
-- the reference after primer trimming
checknonzeroCigMatch :: AlignedRead -> Bool
checknonzeroCigMatch a
    | (isheader a) = True
    | (not $ mapped a) = True
    | (B.any (== 'M') (trimdcigar a)) = True
    | otherwise = False

-- 161017 create I.IntMap Int B.ByteString from (Int, BedRecord)
-- NOTE: this is an attempt at fast lookup of primer intervals
--       key represent chromosome coordinate (IntMap for each chromosome)
makechrbedmap :: BED -> M.Map UChr (I.IntMap BedRecord)
makechrbedmap bs =
    let blist = V.toList bs
        bbychr = groupBy (\x y -> bedchr x == bedchr y) blist -- [[BedRecord]]
        bchrs = bedchr <$> (head <$> bbychr) -- [UChr]
        chrintmaps = makebedmap <$> bbychr -- [I.IntMap BedRecord]
        tups = zipWith (,) bchrs chrintmaps
    in M.fromList tups

makebedmap :: [BedRecord] -> I.IntMap BedRecord
makebedmap bs =
    let tups = join $ makeprimertups <$> bs
    in I.fromList tups

-- expand primer BED interval
-- TODO: make target BED interval padding configurable by command line option
makeprimertups :: BedRecord -> [(Int, BedRecord)]
makeprimertups b =
    let ps = (bedstart b) - 2 -- pad border with arbitrary 2 bases for now
        pe = (bedend b) + 2 -- pad border for now
        rng = intgr2int $ (pe - ps) + 1
        recs = replicate rng $ b
        ints = intgr2int <$> [ps..pe]
    in zipWith (,) ints recs

-- IntMap lookup with default result when key not found (represents no primer
-- at that position)
bedmaplookup :: Integer -> I.IntMap BedRecord -> Maybe BedRecord
bedmaplookup pos bmap = I.lookup (intgr2int pos) bmap

justbedmaps :: [Maybe BedRecord] -> [BedRecord]
justbedmaps ms = catMaybes ms

justchrmaps :: [Maybe (I.IntMap BedRecord)] -> [I.IntMap BedRecord]
justchrmaps mcmaps = catMaybes mcmaps

createprimerbedmaps :: Opts -> IO ( M.Map UChr (I.IntMap BedRecord)
                                  , M.Map UChr (I.IntMap BedRecord) )
createprimerbedmaps args = case (coordfileformat args) of
    False -> do
        m <- getMasterFile $ incoordsfile args
        let fm = makechrbedmap $ masterToFPrimerBED m
            rm = makechrbedmap $ masterToRPrimerBED m
        return (fm, rm)
    True  -> do
        bedpe <- getBEDPE $ incoordsfile args
        let fb = bedpeToFbed <$> bedpe
            rb = bedpeToRbed <$> bedpe
            fmpFromBedPE = makechrbedmap $ V.fromList fb
            rmpFromBedPE = makechrbedmap $ V.fromList rb
        return (fmpFromBedPE, rmpFromBedPE)
    -- _        -> error "(!) Incorrect primer coord input string"

-- 180206 split BEDPE to forward and reverse BedRecords
bedpeToFbed :: BEDPE -> BedRecord
bedpeToFbed b = BedRecord (chr1 b) (start1 b) (end1 b)
                          (B.append (bedpename b) "_F")

bedpeToRbed :: BEDPE -> BedRecord
bedpeToRbed b = BedRecord (chr2 b) (start2 b) (end2 b)
                          (B.append (bedpename b) "_R")

-- 171017 convert AlignedRead to output string format, handling either header
-- information stored in an AlignedRead record or an actual alignment
-- NOTE: this simplifies parsing the input SAM file using Conduit
printAlignmentOrHdr :: AlignedRead -> B.ByteString
printAlignmentOrHdr a
    | (isheader a) = headerstr
    | otherwise = alnstr
        where rn = B.pack $ show $ rname a
              tpos = if (rn == "*") then 0 else (trimdpos a) + 1 -- SAM format
              q = qname a
              f = B.pack $ show $ flag a
              p = B.pack $ show $ tpos
              mq = B.pack $ show $ mapqual a
              c = trimdcigar a -- 170202 prints trimmed alignment CIGAR
              rnxt = rnext a
              pnxt = B.pack $ show $ (pnext a) + 1
              tl = B.pack $ show $ tlen a
              sq = refseq a
              bq = basequal a
              optfs = optfields a
              alnstr = B.intercalate "\t" [ q, f, rn, p, mq, c, rnxt,
                                            pnxt, tl, sq, bq, optfs ]
              headerstr = B.reverse $ B.drop 1 $ B.reverse
                        $ B.unlines $ headerstrings a

printAlignment :: AlignedRead -> B.ByteString
printAlignment a =
    let rn = B.pack $ show $ rname a
        tpos = if (rn == "*") then 0 else (trimdpos a) + 1 -- SAM format rules
        q = qname a
        f = B.pack $ show $ flag a
        p = B.pack $ show $ tpos
        mq = B.pack $ show $ mapqual a
        c = trimdcigar a -- 170202 prints trimmed alignment CIGAR
        rnxt = rnext a
        pnxt = B.pack $ show $ (pnext a) + 1
        tl = B.pack $ show $ tlen a
        sq = refseq a
        bq = basequal a
        optfs = optfields a
        alnstr = B.intercalate "\t" [ q, f, rn, p, mq, c, rnxt,
                                      pnxt, tl, sq, bq, optfs ]
    in alnstr

-- max(0,i)
checkpos :: Integer -> Integer
checkpos i
    | i < 0 = 0
    | otherwise = i

-- 170201 new trimming arithmetic adjusts read coords to reference
adjustcrds :: [(Integer, B.ByteString)] -> [(Integer, B.ByteString)]
adjustcrds cigs = scanl1 shiftcrds cigs

-- 'H' and 'S' should be removed from CIGAR string before calling this function
-- this function adjusts read positions to align with reference (relative)
-- coordinates, allowing correct calculation of the primer-trimmed CIGAR string
shiftcrds :: (Integer, B.ByteString) -> (Integer, B.ByteString)
        -> (Integer, B.ByteString)
shiftcrds a b
    | thiscig == "I" = adjustposI
    | thiscig == "D" = adjustposD
    | otherwise = keeppos
    where lastpos = fst a
          thispos = fst b
          thiscig = snd b
          adjustposI = (lastpos, thiscig)
          adjustposD = (lastpos + 2, thiscig)
          keeppos = (lastpos + 1, thiscig)

taketrim :: Integer -> [(Integer, B.ByteString)] -> [(Integer, B.ByteString)]
taketrim cnt cs = takeWhile (\x -> (fst x) <= cnt) cs

-- 170203 try to fix edge case where D is first base in remaining CIGAR string
trimrem :: Integer -> [(Integer, B.ByteString)] -> [(Integer, B.ByteString)]
trimrem cnt cs
    | length rem == 0 = []
    | (snd $ head rem) == "D" = tail rem
    | otherwise = rem
        where rem = dropWhile (\x -> (fst x) <= cnt) cs

-- 170510 print warning if chromosome name formats from SAM and master input
-- files do not match
checkChromNameMatchStatus :: Header -> BED -> IO ()
checkChromNameMatchStatus hdr bed = do
    let bedchrs = S.toList $ S.fromList $ bedchr <$> (V.toList bed) -- set of BED chrs
        hdrchrs = S.toList $ S.fromList $ U.rights $ A.parseOnly hdrchromp
                <$> ((hdr !!) <$> [1..25]) -- set of SAM hdr chrs
        matches = filter (\x -> x /= NONE) $ intersect bedchrs hdrchrs
    if (length matches) >= 1
        then putStrLn "SAM and BED chromosome name formats match."
        else error "ERROR: different chromosome name formats in SAM and BED!"

-- end of Library
